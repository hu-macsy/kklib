
#include <gtest/gtest.h>

#include "test.hpp"
#include "test_walk.hpp"

#include <kklib/graph.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/storage.hpp>
#include <kklib/util.hpp>
#include <kklib/walk.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <type_traits>
#include <utility>
#include <vector>

typedef uint64_t hash_t;

struct HashWalkState
{
    hash_t hash;
    VertexID previous_vertex;
};

struct HashWalkConf
{
    walker_id_t walker_num;
    step_t walk_length;
    const real_t upper_bound = 3;
    const real_t lower_bound = 1;
    const hash_t magic_num = 10000007;
    hash_t get_walker_init_hash(walker_id_t walker, VertexID vertex) { return walker * magic_num + vertex; }
    hash_t get_walker_new_hash(hash_t old_hash, walker_id_t walker, VertexID next_vertex)
    {
        return old_hash * magic_num + walker ^ next_vertex;
    }
    real_t get_dynamic_comp(hash_t hash, VertexID current, VertexID next)
    {
        real_t temp = hash % 3 + (current ^ next);
        return lower_bound + (real_t) fmod(temp, upper_bound - lower_bound);
    }
};

template<typename edge_data_t>
void hashwalk(WalkEngine<edge_data_t, HashWalkState>* graph, HashWalkConf conf, int order, std::vector<hash_t> &walker_hash)
{
    struct WHItem
    {
        walker_id_t walker;
        hash_t hash;
        step_t step;
        WHItem() : walker(0), hash(0), step(0) {}
        WHItem(walker_id_t _walker, hash_t _hash, step_t _step) : walker(_walker), hash(_hash), step(_step) {}
    };
    std::vector<std::vector<WHItem> > wh_collector(graph->get_worker_num());
    WalkerConfig<edge_data_t, HashWalkState> walker_conf(
        conf.walker_num,
        [&](Walker<HashWalkState>& walker, VertexID start_vertex)
        {
            walker.data.hash = conf.get_walker_init_hash(walker.id, start_vertex);
            wh_collector[omp_get_thread_num()].push_back(WHItem(walker.id, walker.data.hash, walker.step));
        },
        [&](Walker<HashWalkState>& walker, VertexID current_v, AdjUnit<edge_data_t>* edge)
        {
            walker.data.hash = conf.get_walker_new_hash(walker.data.hash, walker.id, edge->neighbour);
            walker.data.previous_vertex = current_v;
            wh_collector[omp_get_thread_num()].push_back(WHItem(walker.id, walker.data.hash, walker.step));
        });
    auto extension_comp = [&](Walker<HashWalkState>& walker, VertexID current_v)
    { return walker.step >= conf.walk_length ? 0.0 : 1.0; };
    auto static_comp = get_trivial_static_comp(graph);
    auto upper_bound_func = [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return conf.upper_bound; };
    auto lower_bound_func = [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return conf.lower_bound; };
    if (order == 1) 
    {
        TransitionConfig<edge_data_t, HashWalkState> tr_conf(
            extension_comp, static_comp,
            [&](Walker<HashWalkState>& walker, VertexID vertex, AdjUnit<edge_data_t>* edge)
            {
                if (walker.step == 0)
                {
                    return conf.upper_bound;
                }
                else
                {
                    return conf.get_dynamic_comp(walker.data.hash, walker.data.previous_vertex, edge->neighbour);
                }
            },
            upper_bound_func, lower_bound_func);
        graph->random_walk(&walker_conf, &tr_conf);
    } else
    {
        SecondOrderTransitionConfig<edge_data_t, HashWalkState, EmptyData, VertexID> tr_conf(
            extension_comp, static_comp,
            [&](Walker<HashWalkState>& walker, walker_id_t walker_idx, VertexID current_v, AdjUnit<edge_data_t>* edge)
            {
                if (walker.step != 0)
                {
                    stateQuery<EmptyData> query;
                    query.src_v = current_v;
                    query.walker_idx = walker_idx;
                    graph->emit(walker.data.previous_vertex, query);
                }
            },
            [&](VertexID vtx, stateQuery<EmptyData> query, AdjList<edge_data_t>* adj_list)
            {
                stateResponse<VertexID> response;
                response.walker_idx = query.walker_idx;
                response.data = vtx;
                graph->emit(query.src_v, response);
            },
            [&](Walker<HashWalkState>& walker, stateResponse<VertexID>& response, VertexID current_v, AdjUnit<edge_data_t>* edge)
            {
                if (walker.step == 0)
                {
                    return conf.upper_bound;
                }
                else
                {
                    return conf.get_dynamic_comp(walker.data.hash, response.data, edge->neighbour);
                }
            },
            upper_bound_func, lower_bound_func);
        graph->random_walk(&walker_conf, &tr_conf);
    }

	std::thread send_thread([&]() {
		for (int t_i = 0; t_i < graph->get_worker_num(); t_i++)
		{
            MPI_Send(wh_collector[t_i].data(), wh_collector[t_i].size() * sizeof(WHItem), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        }
	});
    if (get_mpi_rank() == 0)
    {
        std::vector<WHItem> wh_items;
        for (partition_id_t p_i = 0; p_i < get_mpi_size(); p_i++)
        {
            for (int t_i = 0; t_i < graph->get_worker_num(); t_i++)
            {
                int recv_size = 0;
                MPI_Status recv_status;
                MPI_Probe(p_i, 0, MPI_COMM_WORLD, &recv_status);
                MPI_Get_count(&recv_status, MPI_CHAR, &recv_size);
                int recv_n = recv_size / sizeof(WHItem);
                size_t old_size = wh_items.size();
                wh_items.resize(old_size + recv_n);
                MPI_Recv(wh_items.data() + old_size, recv_size, MPI_CHAR, p_i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        walker_hash.clear();
        walker_hash.resize(conf.walker_num, 0);
        std::vector<step_t> max_step(conf.walker_num, 0);
        for (auto wh : wh_items)
        {
            if (wh.step >= max_step[wh.walker])
            {
                walker_hash[wh.walker] = wh.hash;
                max_step[wh.walker] = wh.step;
            }
        }
    }
    send_thread.join();
}

template <typename edge_data_t>
void check_hashwalk_random_walk(VertexID v_num,
                                Edge<edge_data_t>* edges,
                                edge_id_t e_num,
                                HashWalkConf conf,
                                std::vector<hash_t> real_walker_hash,
                                std::vector<std::vector<VertexID>>& seq)
{
    std::vector<hash_t> std_walker_hash;
    for (walker_id_t w_i = 0; w_i < seq.size(); w_i++)
    {
        hash_t hash = conf.get_walker_init_hash(w_i, seq[w_i][0]); 
        for (size_t s_i = 1; s_i < seq[w_i].size(); s_i++)
        {
            hash = conf.get_walker_new_hash(hash, w_i, seq[w_i][s_i]); 
        }
        std_walker_hash.push_back(hash);
    }
    std::sort(real_walker_hash.begin(), real_walker_hash.end());
    std::sort(std_walker_hash.begin(), std_walker_hash.end());
    ASSERT_EQ(real_walker_hash.size(), std_walker_hash.size());
    for (size_t h_i = 0; h_i < std_walker_hash.size(); h_i++)
    {
        EXPECT_EQ(real_walker_hash[h_i], std_walker_hash[h_i]);
    }
}

template <typename edge_data_t> void test_walker(VertexID v_num, int worker_number, int order)
{
    kklib::LoadedGraphData graph_data = kklib::load_graph<edge_data_t>(v_num, test_data_file);
    WalkEngine<edge_data_t, HashWalkState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
    graph.set_concurrency(worker_number);

    HashWalkConf conf;
    conf.walk_length = 20 + rand() % 20;
    conf.walker_num = graph.get_vertex_num() * 100 + graph.get_edge_num() * 100 + rand() % 100;
    MPI_Bcast(&conf, sizeof(conf), MPI_CHAR, 0, MPI_COMM_WORLD);

    std::vector<hash_t> walker_hash;
    hashwalk(&graph, conf, order, walker_hash);

    std::vector<std::vector<VertexID>> rw_sequences;
    graph.collect_walk_sequence(rw_sequences, conf.walker_num);

    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        check_hashwalk_random_walk(v_num, std_edges, std_edge_num, conf, walker_hash, rw_sequences);
    }
}

template<typename edge_data_t>
void test_walker(int order)
{
    edge_id_t e_nums_arr[] = {200, 400, 600, 800, 1000, 1200, 1400, 1600};
    VertexID v_num = 100 + rand() % 50;
    std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 8);
    /*
    size_t e_nums_arr[] = {30};
    VertexID v_num = 10;
    std::vector<size_t> e_nums(e_nums_arr, e_nums_arr + 1);
    */

    MPI_Bcast(&v_num, 1, kklib::deduce_mpi_data_type<VertexID>(), 0, MPI_COMM_WORLD);

    for (auto &e_num : e_nums_arr)
    {
        if (get_mpi_rank() == 0)
        {
            gen_undirected_graph_file<edge_data_t>(v_num, e_num);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int worker_number = rand() % 8 + 1;
        MPI_Bcast(&worker_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
        test_walker<edge_data_t>(v_num, worker_number, order);
    }
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(Walker, UnbiasedFirstOrder)
{
    test_walker<EmptyData>(1);
}

TEST(Walker, BiasedFirstOrder)
{
    test_walker<real_t>(1);
}

TEST(Walker, UnbiasedSecondOrder)
{
    test_walker<EmptyData>(2);
}

TEST(Walker, BiasedSecondOrder)
{
    test_walker<EmptyData>(2);
}


GTEST_API_ int main(int argc, char *argv[])
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}
