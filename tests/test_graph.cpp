#include "test.hpp"

#include <kklib/graph.hpp>
#include <kklib/next_storage.hpp>
#include <kklib/util.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

template<typename edge_data_t>
class GraphTester: public GraphEngine<edge_data_t>
{
public:
    GraphTester(VertexID _v_num,
                Edge<edge_data_t>* read_edges,
                edge_id_t read_e_num,
                partition_id_t const _partition_num = get_mpi_size(),
                partition_id_t const _local_partition_id = get_mpi_rank())
    : GraphEngine<edge_data_t>(_v_num, read_edges, read_e_num, _partition_num, _local_partition_id)
    {
    }

    void get_edges(EdgeContainer<edge_data_t> *ec, std::vector<Edge<edge_data_t>> &ret)
    {
        ret.clear();
        for (VertexID v_i = this->vertex_partition_begin[this->local_partition_id];
             v_i < this->vertex_partition_end[this->local_partition_id]; v_i++)
        {
            for (auto p = ec->adj_lists[v_i].begin; p != ec->adj_lists[v_i].end; p++)
            {
                Edge<edge_data_t> e;
                e.src = v_i;
                e.dst = p->neighbour;
				if (!std::is_same<edge_data_t, EmptyData>::value)
				{
					e.data = p->data;
				}
                ret.push_back(e);
            }
        }
    }

    void set_concurrency(int worker_number)
    {
        this->set_graph_engine_concurrency(worker_number);
    }
};

template <typename edge_data_t> void check_edges(GraphTester<edge_data_t>& graph, bool load_as_undirected = false)
{
    std::vector<Edge<edge_data_t> > local_graph_edges;
    graph.get_edges(graph.csr, local_graph_edges);
    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        if (load_as_undirected)
        {
            std::vector<Edge<edge_data_t> > temp;
            for (edge_id_t e_i = 0; e_i < std_edge_num; e_i++)
            {
                temp.push_back(std_edges[e_i]);
                std::swap(std_edges[e_i].src, std_edges[e_i].dst);
                temp.push_back(std_edges[e_i]);
            }
            delete []std_edges;
            std_edge_num *= 2;
            std_edges = new Edge<edge_data_t>[std_edge_num];
            memcpy(std_edges, temp.data(), sizeof(Edge<edge_data_t>) * std_edge_num);
        }
        auto graph_edges = local_graph_edges;
        for (partition_id_t p_i = 1; p_i < get_mpi_size(); p_i++)
        {
            int recv_size = 0;
            MPI_Status recv_status;
            MPI_Probe(p_i, Tag_ShuffleGraph, MPI_COMM_WORLD, &recv_status);
            MPI_Get_count(&recv_status, MPI_CHAR, &recv_size);
            std::vector<Edge<edge_data_t> > remote_edges(recv_size / sizeof(Edge<edge_data_t>));
            MPI_Recv(remote_edges.data(), recv_size, MPI_CHAR, p_i, Tag_ShuffleGraph, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (auto e : remote_edges)
            {
                graph_edges.push_back(e);
            }
        }
        cmp_edges(graph_edges.data(), graph_edges.size(), std_edges, std_edge_num); 
        delete []std_edges;
    } else
    {
        MPI_Send(local_graph_edges.data(), local_graph_edges.size() * sizeof(Edge<edge_data_t>), MPI_CHAR, 0,
                 Tag_ShuffleGraph, MPI_COMM_WORLD);
    }
}

template <typename edge_data_t> void test_static_edge(VertexID v_num, bool load_as_undirected = false)
{
    kklib::LoadedGraphData graph_data = kklib::load_graph<edge_data_t>(v_num, test_data_file);
    GraphTester<edge_data_t> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };

    int worker_number = rand() % 8 + 1;
    graph.set_concurrency(worker_number);

    check_edges(graph, load_as_undirected);
}

template <typename edge_data_t> void test_edges(bool const load_as_undirected = false)
{
    // edge_id_t e_nums_arr[] = {0, 2, 6, 16, 8888, 10000, 20000, 100000};
    // VertexID v_num = 1000 + rand() % 1000;
    // std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 8);

    std::array<size_t, 1> e_nums_arr = { 20 };
    VertexID v_num = 20;


    for (auto e_num : e_nums_arr)
    {
        if (get_mpi_rank() == 0)
        {
            if (load_as_undirected)
            {
                gen_directed_graph_file<edge_data_t>(v_num, e_num);
            } else
            {
                gen_undirected_graph_file<edge_data_t>(v_num, e_num);
            }
        }

        MPI_Bcast(&v_num, 1, kklib::deduce_mpi_data_type<VertexID>(), 0, MPI_COMM_WORLD);
        test_static_edge<edge_data_t>(v_num, load_as_undirected);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(GraphEngine, DefaultLoad)
{
    test_edges<EmptyData>();
    test_edges<real_t>();
}

TEST(GraphEngine, LoadAsUndirected)
{
    test_edges<EmptyData>(true);
    test_edges<real_t>(true);
}

GTEST_API_ int main(int argc, char *argv[])
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}

