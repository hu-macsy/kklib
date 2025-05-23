#pragma once

#include <kklib/constants.hpp>
#include <kklib/graph_data_structure.hpp>
#include <kklib/mpi_helper.hpp>
#include <kklib/type.hpp>
#include <kklib/util.hpp>

#include <omp.h>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <type_traits>
#include <thread>
#include <random>
#include <functional>
#include <unistd.h>
#include <sys/mman.h>
#include <limits.h>

template<typename edge_data_t>
class GraphEngine
{
protected:
    VertexID v_num;
    edge_id_t e_num;
    int worker_num;
    edge_id_t local_e_num;

    VertexID* vertex_partition_begin;
    VertexID* vertex_partition_end;

    partition_id_t local_partition_id;
    partition_id_t partition_num;

    MessageBuffer **thread_local_msg_buffer; 
    MessageBuffer **msg_send_buffer;
    MessageBuffer **msg_recv_buffer;
    std::mutex *send_locks;
    std::mutex *recv_locks;

    DistributedExecutionCtx dist_exec_ctx;
public:
    VertexID* vertex_in_degree;
    VertexID* vertex_out_degree;
    partition_id_t *vertex_partition_id;

    EdgeContainer<edge_data_t> *csr;

protected:
    void set_graph_engine_concurrency(int worker_num_param)
    {
        this->worker_num = worker_num_param;
        omp_set_dynamic(0);
        omp_set_num_threads(worker_num);
        //message buffer depends on worker number
        free_msg_buffer();
    }

public:
    inline bool is_local_vertex(VertexID v_id)
    {
        return v_id >= vertex_partition_begin[local_partition_id]
            && v_id < vertex_partition_end[local_partition_id];
    }
    inline bool is_valid_edge(Edge<edge_data_t> e)
    {
        return e.src < v_num && e.dst < v_num;
    }
    inline VertexID get_vertex_num() { return v_num; }
    inline edge_id_t get_edge_num()
    {
        return e_num;
    }
    inline int get_worker_num()
    {
        return worker_num;
    }
    inline VertexID get_local_vertex_begin() { return vertex_partition_begin[local_partition_id]; }
    inline VertexID get_local_vertex_end() { return vertex_partition_end[local_partition_id]; }
    inline VertexID get_vertex_begin(partition_id_t p) { return vertex_partition_begin[p]; }
    inline VertexID get_vertex_end(partition_id_t p) { return vertex_partition_end[p]; }

public:
    // deallocate a vertex array
    template<typename T>
    void dealloc_vertex_array(T * array)
    {
        dealloc_array(array, v_num);
    }

    template<typename T>
    void dealloc_array(T * array, size_t num)
    {
        munmap(array, sizeof(T) * num);
    }

    // allocate a vertex array
    template<typename T>
    T * alloc_vertex_array()
    {
        return alloc_array<T>(v_num);
    }

    template<typename T>
    T * alloc_array(size_t num)
    {
        T* array = (T*) mmap(NULL, sizeof(T) * num, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        assert(array != nullptr);
        return array;
    }

    GraphEngine()
    {
        vertex_partition_begin = nullptr;
        vertex_partition_end = nullptr;

        thread_local_msg_buffer = nullptr;
        msg_send_buffer = nullptr;
        msg_recv_buffer = nullptr;

        send_locks = nullptr;
        recv_locks = nullptr;

        vertex_in_degree = nullptr;
        vertex_out_degree = nullptr;
        vertex_partition_id = nullptr;

        csr = nullptr;

        this->worker_num = std::max(1, ((int)std::thread::hardware_concurrency()) - 1);
        omp_set_dynamic(0);
        omp_set_num_threads(worker_num);
    }

    GraphEngine(VertexID _v_num,
                Edge<edge_data_t>* read_edges,
                edge_id_t read_e_num,
                partition_id_t const _partition_num = get_mpi_size(),
                partition_id_t const _local_partition_id = get_mpi_rank())
    : GraphEngine()
    {
        v_num = _v_num;
        partition_num = _partition_num;
        local_partition_id = _local_partition_id;
        local_e_num = 0;

        vertex_in_degree = alloc_vertex_array<VertexID>();
        vertex_out_degree = alloc_vertex_array<VertexID>();

        std::vector<VertexID> local_vertex_degree(v_num, 0);
        for (edge_id_t e_i = 0; e_i < read_e_num; e_i++)
        {
            local_vertex_degree[read_edges[e_i].src]++;
        }

        MPI_Allreduce(local_vertex_degree.data(), vertex_out_degree, v_num, kklib::deduce_mpi_data_type<VertexID>(),
                      MPI_SUM, MPI_COMM_WORLD);

        std::fill(local_vertex_degree.begin(), local_vertex_degree.end(), 0);
        for (edge_id_t e_i = 0; e_i < read_e_num; e_i++)
        {
            local_vertex_degree[read_edges[e_i].dst]++;
        }
        MPI_Allreduce(local_vertex_degree.data(), vertex_in_degree, v_num, kklib::deduce_mpi_data_type<VertexID>(),
                      MPI_SUM, MPI_COMM_WORLD);

        vertex_partition_begin = new VertexID[partition_num];
        vertex_partition_end = new VertexID[partition_num];

        edge_id_t total_workload = 0;
        e_num = 0;
        for (VertexID v_i = 0; v_i < v_num; v_i++)
        {
            total_workload += 5 + vertex_out_degree[v_i];
            e_num += vertex_out_degree[v_i];
        }
        edge_id_t workload_per_node = (total_workload + partition_num - 1) / partition_num;
        for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
        {
            if (p_i == 0)
            {
                vertex_partition_begin[p_i] = 0;
            }
            else
            {
                vertex_partition_begin[p_i] = vertex_partition_end[p_i - 1];
            }
            vertex_partition_end[p_i] = vertex_partition_begin[p_i];
            edge_id_t workload = 0;
            for (VertexID v_i = vertex_partition_begin[p_i]; v_i < v_num && workload < workload_per_node; v_i++)
            {
                workload += 5 + vertex_out_degree[v_i];
                vertex_partition_end[p_i]++;
            }
#ifdef PERF_PROF
            if (local_partition_id == 0)
            {
                printf("partition %d: %u %u (%zu %zu)\n", p_i, vertex_partition_begin[p_i], vertex_partition_end[p_i],
                       workload, workload_per_node);
            }
#endif
        }
        assert(vertex_partition_end[partition_num - 1] == v_num);

        vertex_partition_id = alloc_vertex_array<partition_id_t>();
        for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
        {
            for (VertexID v_i = vertex_partition_begin[p_i]; v_i < vertex_partition_end[p_i]; v_i++)
            {
                vertex_partition_id[v_i] = p_i;
            }
        }

        for (VertexID v_i = vertex_partition_begin[local_partition_id]; v_i < vertex_partition_end[local_partition_id]; v_i++)
        {
            local_e_num += vertex_out_degree[v_i];
        }
        Edge<edge_data_t>* local_edges = new Edge<edge_data_t>[local_e_num];

        //! Shuffle edges to remove any pattern when reading in file.
        shuffle_edges(read_edges, read_e_num, local_edges, local_e_num);

        csr = build_edge_container(local_edges, local_e_num, vertex_out_degree, v_num, local_partition_id,
                                   vertex_partition_begin, vertex_partition_end);

        send_locks = new std::mutex[partition_num];
        recv_locks = new std::mutex[partition_num];

        dist_exec_ctx.progress = new size_t*[worker_num];
        for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
        {
            dist_exec_ctx.progress[t_i] = new size_t[partition_num];
        }

#ifdef PERF_PROF
        printf("finish build graph, time %.3lfs\n", timer.duration());
#endif
        //! Handle this deletion differently.
        delete[] read_edges;
        delete[] local_edges;
    }

    virtual ~GraphEngine()
    {
        if (vertex_partition_begin != nullptr)
        {
            delete []vertex_partition_begin;
        }
        if (vertex_partition_end != nullptr)
        {
            delete []vertex_partition_end;
        }

        if (send_locks != nullptr)
        {
            delete []send_locks;
        }
        if (recv_locks != nullptr)
        {
            delete []recv_locks;
        }

        if (vertex_in_degree != nullptr)
        {
            dealloc_vertex_array<VertexID>(vertex_in_degree);
        }
        if (vertex_out_degree != nullptr)
        {
            dealloc_vertex_array<VertexID>(vertex_out_degree);
        }
        if (vertex_partition_id != nullptr)
        {
            dealloc_vertex_array<partition_id_t>(vertex_partition_id);
        }

        if (csr != nullptr)
        {
            delete csr;
        }

        if (dist_exec_ctx.progress != nullptr)
        {
            for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
            {
                delete []dist_exec_ctx.progress[t_i];
            }
            delete []dist_exec_ctx.progress;
        }

        free_msg_buffer();
    }

    void shuffle_edges(Edge<edge_data_t> *misc_edges, edge_id_t misc_e_num, Edge<edge_data_t> *local_edges, edge_id_t local_e_num)
    {
        std::vector<edge_id_t> e_count(partition_num, 0);
        for (edge_id_t e_i = 0; e_i < misc_e_num; e_i++)
        {
            e_count[vertex_partition_id[misc_edges[e_i].src]]++;
        }
        Edge<edge_data_t> *tmp_es  = new Edge<edge_data_t>[misc_e_num];
        std::vector<edge_id_t> e_p(partition_num, 0);
        for (partition_id_t p_i = 1; p_i < partition_num; p_i++)
        {
            e_p[p_i] = e_p[p_i - 1] + e_count[p_i - 1];
        }
        auto e_begin = e_p;
        for (edge_id_t e_i = 0; e_i < misc_e_num; e_i++)
        {
            auto pt = vertex_partition_id[misc_edges[e_i].src];
            tmp_es[e_p[pt] ++] = misc_edges[e_i];
        }
        edge_id_t local_edge_p = 0;
        std::thread send_thread([&](){
            for (partition_id_t step = 0; step < partition_num; step++)
            {
                partition_id_t dst = (local_partition_id + step) % partition_num;
                size_t tot_send_sz = e_count[dst] * sizeof(Edge<edge_data_t>);
#ifdef UNIT_TEST
                const int max_single_send_sz = (1 << 8) / sizeof(Edge<edge_data_t>) * sizeof(Edge<edge_data_t>);
#else
                const int max_single_send_sz = (1 << 28) / sizeof(Edge<edge_data_t>) * sizeof(Edge<edge_data_t>);
#endif
                void* send_data = tmp_es + e_begin[dst];
                while (true)
                {
                    MPI_Send(&tot_send_sz, 1, kklib::deduce_mpi_data_type<size_t>(), dst, Tag_ShuffleGraph, MPI_COMM_WORLD);
                    if (tot_send_sz == 0)
                    {
                        break;
                    }
                    int send_sz = std::min((size_t)max_single_send_sz, tot_send_sz);
                    tot_send_sz -= send_sz;
                    MPI_Send(send_data, send_sz, MPI_CHAR, dst, Tag_ShuffleGraph, MPI_COMM_WORLD);
                    send_data = (char*)send_data + send_sz;
                }
                usleep(100000);
            }
        });
        std::thread recv_thread([&](){
            for (partition_id_t step = 0; step < partition_num; step++)
            {
                partition_id_t src = (local_partition_id + partition_num - step) % partition_num;
                while (true)
                {
                    size_t remained_sz;
                    MPI_Recv(&remained_sz, 1, kklib::deduce_mpi_data_type<size_t>(), src, Tag_ShuffleGraph,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    if (remained_sz == 0)
                    {
                        break;
                    }

                    MPI_Status recv_status;
                    MPI_Probe(src, Tag_ShuffleGraph, MPI_COMM_WORLD, &recv_status);
                    int sz;
                    MPI_Get_count(&recv_status, MPI_CHAR, &sz);

                    MPI_Recv(local_edges + local_edge_p, sz, MPI_CHAR, src, Tag_ShuffleGraph, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    local_edge_p += sz / sizeof(Edge<edge_data_t>);
                }
                usleep(100000);
            }
        });
        send_thread.join();
        recv_thread.join();
        delete []tmp_es;
        assert(local_e_num == local_edge_p);
    }

    void set_msg_buffer(size_t max_msg_num, size_t max_msg_size)
    {
        if (thread_local_msg_buffer == nullptr)
        {
            thread_local_msg_buffer = new MessageBuffer*[worker_num];
            #pragma omp parallel
            {
                int worker_id = omp_get_thread_num();
                thread_local_msg_buffer[worker_id] = new MessageBuffer();
            }
        }
        if (msg_send_buffer == nullptr)
        {
            msg_send_buffer = new MessageBuffer*[partition_num];
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                msg_send_buffer[p_i] = new MessageBuffer();
            }
        }
        if (msg_recv_buffer == nullptr)
        {
            msg_recv_buffer = new MessageBuffer*[partition_num];
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                msg_recv_buffer[p_i] = new MessageBuffer();
            }
        }

        size_t local_buf_size = max_msg_size * kklib::thread_local_buf_capacity;
#pragma omp parallel
        {
            int worker_id = omp_get_thread_num();
            if (thread_local_msg_buffer[worker_id]->sz < local_buf_size)
            {
                thread_local_msg_buffer[worker_id]->alloc(local_buf_size);
            }
        }
        size_t comm_buf_size = max_msg_size * max_msg_num;
        for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
        {
            if (msg_send_buffer[p_i]->sz < comm_buf_size)
            {
                msg_send_buffer[p_i]->alloc(comm_buf_size);
            }
            if (msg_recv_buffer[p_i]->sz < comm_buf_size)
            {
                msg_recv_buffer[p_i]->alloc(comm_buf_size);
            }
        }
    }

    void free_msg_buffer()
    {
        if (thread_local_msg_buffer != nullptr)
        {
            for (partition_id_t t_i = 0; t_i < worker_num; t_i ++)
            {
                delete thread_local_msg_buffer[t_i];
            }
            delete []thread_local_msg_buffer;
        }
        if (msg_send_buffer != nullptr)
        {
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                delete msg_send_buffer[p_i];
            }
            delete []msg_send_buffer;
        }
        if (msg_recv_buffer != nullptr)
        {
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                delete msg_recv_buffer[p_i];
            }
            delete []msg_recv_buffer;
        }
    }

    template <typename msg_data_t> void emit(VertexID dst_id, msg_data_t data, int worker_id)
    {
        typedef Message<msg_data_t> msg_t;
        msg_t* buf_data = (msg_t*)thread_local_msg_buffer[worker_id]->data;
        auto &count = thread_local_msg_buffer[worker_id]->count;
        buf_data[count].dst_vertex_id = dst_id;
        buf_data[count].data = data;
        count++;
#ifdef UNIT_TEST
        thread_local_msg_buffer[worker_id]->self_check<msg_t>();
        assert(dst_id < v_num);
#endif
        if (count == kklib::thread_local_buf_capacity)
        {
            flush_thread_local_msg_buffer<msg_t>(worker_id);
        }
    }


    template <typename msg_data_t> void emit(VertexID dst_id, msg_data_t data)
    {
        emit(dst_id, data, omp_get_thread_num());
    }

    template<typename msg_t>
    void flush_thread_local_msg_buffer(partition_id_t worker_id)
    {
        auto local_buf = thread_local_msg_buffer[worker_id];
        msg_t *local_data = (msg_t*)local_buf->data;
        auto &local_msg_count = local_buf->count;
        if (local_msg_count != 0)
        {
            VertexID dst_count[partition_num];
            std::fill(dst_count, dst_count + partition_num, 0);
            for (VertexID m_i = 0; m_i < local_msg_count; m_i++)
            {
                dst_count[vertex_partition_id[local_data[m_i].dst_vertex_id]] ++;
            }
            msg_t *dst_data_pos[partition_num];
            size_t end_data_pos[partition_num];
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                VertexID start_pos = __sync_fetch_and_add(&msg_send_buffer[p_i]->count, dst_count[p_i]);
#ifdef UNIT_TEST
                msg_send_buffer[p_i]->self_check<msg_t>();
#endif
                dst_data_pos[p_i] = (msg_t*)(msg_send_buffer[p_i]->data) + start_pos;
                end_data_pos[p_i] = start_pos + dst_count[p_i];
            }
            for (VertexID m_i = 0; m_i < local_msg_count; m_i++)
            {
                *(dst_data_pos[vertex_partition_id[local_data[m_i].dst_vertex_id]]++) =  local_data[m_i];
            }
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                dist_exec_ctx.progress[worker_id][p_i] = end_data_pos[p_i];
            }
            local_msg_count = 0;
        }
    }

    void notify_progress(VertexID progress_begin, VertexID progress_end, VertexID workload, bool phased_exec)
    {
        int phase_num = phased_exec ? kklib::distributed_execution_ctx_phasenum : 1;
        if (phase_num > 1)
        {
            VertexID work_per_phase = workload / phase_num + 1;
            int phase_begin = 0;
            while (progress_begin >= work_per_phase)
            {
                phase_begin ++;
                progress_begin -= work_per_phase;
            }
            int phase_end = 0;
            while (progress_end >= work_per_phase)
            {
                phase_end ++;
                progress_end -= work_per_phase;
            }
            if (phase_end == phase_num)
            {
                phase_end --;
            }
            for (int phase_i = phase_begin; phase_i < phase_end; phase_i++)
            {
                dist_exec_ctx.phase_locks[phase_i].unlock();
                __sync_fetch_and_add(&dist_exec_ctx.unlocked_phase, 1);
            }
        }
    }

    template<typename msg_data_t>
    size_t distributed_execute(
        std::function<void(void)> msg_producer,
        std::function<void(Message<msg_data_t> *, Message<msg_data_t> *)> msg_consumer,
        Message<msg_data_t> *zero_copy_data = nullptr,
        bool phased_exec = false
    )
    {
        typedef Message<msg_data_t> msg_t;
        Timer timer;

        int phase_num = phased_exec ? kklib::distributed_execution_ctx_phasenum : 1;
        for (int phase_i = 0; phase_i < phase_num; phase_i++)
        {
            dist_exec_ctx.phase_locks[phase_i].lock();
        }
        for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
        {
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                dist_exec_ctx.progress[t_i][p_i] = 0;
            }
        }
        dist_exec_ctx.unlocked_phase = 0;
        for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
        {
            recv_locks[p_i].lock();
        }
        volatile size_t zero_copy_recv_count = 0;
        std::thread recv_thread([&](){
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                msg_recv_buffer[p_i]->count = 0;
            }
            std::vector<MPI_Request*> requests[partition_num];
            auto recv_func = [&] (partition_id_t src)
            {
                MPI_Status prob_status;
                MPI_Probe(src, Tag_Msg, MPI_COMM_WORLD, &prob_status);
                int sz;
                MPI_Get_count(&prob_status, MPI_CHAR, &sz);
                //printf("recv %u <- %u: %zu\n", local_partition_id, src, sz / sizeof(msg_t));
                MPI_Request *recv_req = new MPI_Request();
                requests[src].push_back(recv_req);
                if (zero_copy_data == nullptr)
                {
                    MPI_Irecv(((msg_t*)msg_recv_buffer[src]->data) + msg_recv_buffer[src]->count, sz, MPI_CHAR, src,
                              Tag_Msg, MPI_COMM_WORLD, recv_req);
                    msg_recv_buffer[src]->count += sz / sizeof(msg_t);
                    msg_recv_buffer[src]->template self_check<msg_t>();
                } else
                {
                    MPI_Irecv(zero_copy_data + zero_copy_recv_count, sz, MPI_CHAR, src, Tag_Msg, MPI_COMM_WORLD, recv_req);
                    zero_copy_recv_count += sz / sizeof(msg_t);
                }
            };
            for (int phase_i = 0; phase_i < phase_num; phase_i ++)
            {
                if (phase_i + 1 == phase_num)
                {
                    for (partition_id_t step = 0; step < partition_num; step++)
                    {
                        partition_id_t src = (partition_num + local_partition_id - step) % partition_num;
                        recv_func(src);
                    }
                } else
                {
                    partition_id_t src = (partition_num + local_partition_id - phase_i % partition_num) % partition_num;
                    recv_func(src);
                }
            }
            for (partition_id_t step = 0; step < partition_num; step++)
            {
                partition_id_t src = (partition_num + local_partition_id - step) % partition_num;
                for (auto req : requests[src])
                {
                    MPI_Status status;
                    MPI_Wait(req, &status);
                    delete req;
                }
                recv_locks[src].unlock();
            }
        });

        std::thread send_thread([&](){
            size_t send_progress[partition_num];
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                send_progress[p_i] = 0;
            }
            std::vector<MPI_Request*> requests;
            auto send_func = [&] (partition_id_t dst, size_t diff)
            {
                msg_send_buffer[dst]->template self_check<msg_t>();
                MPI_Request* req = new MPI_Request();
                requests.push_back(req);
                MPI_Isend(((msg_t*)msg_send_buffer[dst]->data) + send_progress[dst], diff * sizeof(msg_t), MPI_CHAR,
                          dst, Tag_Msg, MPI_COMM_WORLD, req);
#ifdef PERF_PROF
                if (local_partition_id == 0)
                {
                    printf("end send %u -> %u: %zu time %lf\n", local_partition_id, dst, diff, timer.duration());
                }
#endif
                send_progress[dst] += diff;
            };
            for (int phase_i = 0; phase_i < phase_num; phase_i++)
            {
                dist_exec_ctx.phase_locks[phase_i].lock();
                if (phase_i + 1 == phase_num)
                {
                    for (partition_id_t step = 0; step < partition_num; step++)
                    {
                        partition_id_t dst = (local_partition_id + step) % partition_num;
                        size_t max_progress = 0;
                        for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
                        {
                            volatile size_t temp_val = dist_exec_ctx.progress[t_i][dst];
                            if (temp_val > max_progress)
                            {
                                max_progress = temp_val;
                            }
                        }
                        size_t diff = max_progress - send_progress[dst];
                        send_func(dst, diff);
                    }
                } else
                {
                    partition_id_t dst = (local_partition_id + phase_i) % partition_num;
                    size_t min_progress = UINT_MAX;
                    for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
                    {
                        volatile size_t temp_val = dist_exec_ctx.progress[t_i][dst];
                        if (temp_val < min_progress)
                        {
                            min_progress = temp_val;
                        }
                    }
                    size_t diff = min_progress - send_progress[dst];
                    send_func(dst, diff);
                }
                dist_exec_ctx.phase_locks[phase_i].unlock();
            }
            for (partition_id_t p_i = 0; p_i < partition_num; p_i++)
            {
                msg_send_buffer[p_i]->count = 0;
            }
            for (auto req : requests)
            {
                MPI_Status status;
                MPI_Wait(req, &status);
                delete req;
            }
        });

        msg_producer();

        size_t flush_workload = 0;
        for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
        {
            flush_workload += thread_local_msg_buffer[t_i]->count;
        }
#pragma omp parallel for if (flush_workload * 2 >= kklib::omp_parallel_threshold)
        for (partition_id_t t_i = 0; t_i < worker_num; t_i++)
        {
            flush_thread_local_msg_buffer<msg_t>(t_i);
        }

        for (int phase_i = dist_exec_ctx.unlocked_phase; phase_i < phase_num; phase_i++)
        {
            dist_exec_ctx.phase_locks[phase_i].unlock();
        }
#ifdef PERF_PROF
        if (local_partition_id == 0)
        {
            printf("%u: finish msg_producer in %lfs\n", local_partition_id, timer.duration());
        }
#endif

        size_t msg_num = 0;
        for (int step = 0; step < partition_num; step++)
        {
            partition_id_t src_partition_id = (partition_num + local_partition_id - step) % partition_num;
            recv_locks[src_partition_id].lock();
            if (zero_copy_data == nullptr)
            {
                size_t data_amount = msg_recv_buffer[src_partition_id]->count;
                msg_num += data_amount;
                msg_t* data_begin = (msg_t*)(msg_recv_buffer[src_partition_id]->data);
                msg_t* data_end = data_begin + data_amount;
                msg_consumer(data_begin, data_end);
            }
            recv_locks[src_partition_id].unlock();
        }
        if (zero_copy_data != nullptr)
        {
            msg_consumer(zero_copy_data, zero_copy_data + zero_copy_recv_count);
            msg_num = zero_copy_recv_count;
        }
#ifdef PERF_PROF
        if (local_partition_id == 0)
        {
            printf("%u: finish msg_consumer in %lfs\n", local_partition_id, timer.duration());
        }
#endif
        recv_thread.join();
        send_thread.join();
#ifdef PERF_PROF
        if (local_partition_id == 0)
        {
            printf("%u: finish transmission in %lfs\n", local_partition_id, timer.duration());
        }
#endif
        size_t glb_msg_num;
        MPI_Allreduce(&msg_num, &glb_msg_num, 1, kklib::deduce_mpi_data_type<size_t>(), MPI_SUM, MPI_COMM_WORLD);
        return glb_msg_num;
    }

    template <typename reducer_data_t> reducer_data_t process_vertices(std::function<reducer_data_t(VertexID)> process)
    {
        VertexID progress = vertex_partition_begin[local_partition_id];
        VertexID step_length = kklib::parallel_chunk_size;
        reducer_data_t reducer = 0;
#pragma omp parallel reduction(+:reducer)
        {
            VertexID work_begin, work_end;
            while ((work_begin = __sync_fetch_and_add(&progress, step_length)) < vertex_partition_end[local_partition_id])
            {
                VertexID work_end = std::min(work_begin + step_length, vertex_partition_end[local_partition_id]);
                for (VertexID v_i = work_begin; v_i != work_end; v_i++)
                {
                    reducer += process(v_i);
                }
            }
        }
        reducer_data_t glb_reducer;
        MPI_Allreduce(&reducer, &glb_reducer, 1, kklib::deduce_mpi_data_type<reducer_data_t>(), MPI_SUM, MPI_COMM_WORLD);
        return glb_reducer;
    }
};
