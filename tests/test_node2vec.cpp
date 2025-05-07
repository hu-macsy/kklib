#include <gtest/gtest.h>

#include "test.hpp"
#include "test_walk.hpp"

#include <kklib/graph.hpp>
#include <kklib/node2vec.hpp>
#include <kklib/storage.hpp>
#include <kklib/util.hpp>
#include <kklib/walk.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

template <typename edge_data_t>
void get_node2vec_trans_matrix(
    VertexID v_num, Edge<edge_data_t>* edges, edge_id_t e_num, double p, double q, std::vector<std::vector<double>>& trans_mat)
{
    std::vector<std::vector<Edge<edge_data_t> > > graph(v_num);
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        graph[edges[e_i].src].push_back(edges[e_i]);
    }
    for (VertexID v_i = 0; v_i < v_num; v_i++)
    {
        std::sort(graph[v_i].begin(), graph[v_i].end(), [](const Edge<edge_data_t> a, const Edge<edge_data_t> b){return a.dst < b.dst;});
    }
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        VertexID src = edges[e_i].src;
        VertexID dst = edges[e_i].dst;
        assert(src != dst);
        //must be undirected graph
        assert(graph[dst].size() != 0);
        for (auto e : graph[dst])
        {
            if (e.dst == src)
            {
                trans_mat[e_i][e.dst] += 1 / p * get_edge_trans_weight(e);
            } else if (std::binary_search(graph[src].begin(), graph[src].end(), e, [](const Edge<edge_data_t> a, const Edge<edge_data_t> b){return a.dst < b.dst;}))
            {
                trans_mat[e_i][e.dst] += 1 * get_edge_trans_weight(e);
            } else
            {
                trans_mat[e_i][e.dst] += 1 / q * get_edge_trans_weight(e);
            }
        }
    }
    mat_normalization(trans_mat);
}

template <typename edge_data_t>
void check_node2vec_random_walk(
    VertexID v_num, Edge<edge_data_t>* edges, edge_id_t e_num, double p, double q, std::vector<std::vector<VertexID>> rw_sequences)
{
    std::vector<std::vector<double> > trans_mat(e_num);
    for (auto &vec : trans_mat)
    {
        vec.resize(v_num, 0.0);
    }
    get_node2vec_trans_matrix(v_num, edges, e_num, p, q, trans_mat);

    //check if sequences are legal
    std::vector<VertexID> out_degree(v_num, 0);
    std::vector<std::vector<bool> > adj_mat(v_num);
    for (auto &vec : adj_mat)
    {
        vec.resize(v_num, false);
    }
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        adj_mat[edges[e_i].src][edges[e_i].dst] = true;
        out_degree[edges[e_i].src]++;
    }
    for (auto &s : rw_sequences)
    {
        if (out_degree[s[0]] == 0)
        {
            for (auto v : s)
            {
                ASSERT_EQ(v, s[0]);
            }
        } else
        {
            for (size_t v_i = 0; v_i + 1 < s.size(); v_i++)
            {
                if (adj_mat[s[v_i]][s[v_i + 1]] == false)
                {
                    printf("fault %u %u\n", s[v_i], s[v_i + 1]);
                }
                ASSERT_TRUE(adj_mat[s[v_i]][s[v_i + 1]]);
            }
        }
    }

    std::map<std::pair<VertexID, VertexID>, edge_id_t> dict;
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        std::pair<VertexID, VertexID> key = std::pair<VertexID, VertexID>(edges[e_i].src, edges[e_i].dst);
        assert(dict.find(key) == dict.end());
        dict[key] = e_i;
    }

    std::vector<std::vector<double> > real_trans_mat(e_num);
    for (auto &vec : real_trans_mat)
    {
        vec.resize(v_num, 0.0);
    }
    for (auto &s : rw_sequences)
    {
        if (out_degree[s[0]] != 0)
        {
            for (size_t v_i = 0; v_i + 2 < s.size(); v_i++)
            {
                real_trans_mat[dict[std::pair<VertexID, VertexID>(s[v_i], s[v_i + 1])]][s[v_i + 2]] += 1;
            }
        }
    }
    mat_normalization(real_trans_mat);
    cmp_trans_matrix(real_trans_mat, trans_mat, 10.0);
}

template <typename edge_data_t> void test_node2vec(VertexID v_num, int worker_number)
{
    kklib::LoadedGraphData graph_data = kklib::load_graph<edge_data_t>(v_num, test_data_file);
    WalkEngine<edge_data_t, Node2vecState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
    graph.set_concurrency(worker_number);

    Node2vecConf n2v_conf;
    n2v_conf.walk_length = 80 + rand() % 20;
    n2v_conf.walker_num = graph.get_vertex_num() * 500 + graph.get_edge_num() * 100 + rand() % 100;
    n2v_conf.p = rand() % 4 + 1;
    n2v_conf.q = rand() % 4 + 1;
    if (rand() % 2 == 0)
    {
        n2v_conf.p = 1.0 / n2v_conf.p;
    }
    if (rand() % 2 == 0)
    {
        n2v_conf.q = 1.0 / n2v_conf.q;
    }
    MPI_Bcast(&n2v_conf, sizeof(n2v_conf), MPI_CHAR, 0, MPI_COMM_WORLD);

    node2vec(&graph, n2v_conf);
    std::vector<std::vector<VertexID>> rw_sequences;
    graph.collect_walk_sequence(rw_sequences, n2v_conf.walker_num);

    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        check_node2vec_random_walk(v_num, std_edges, std_edge_num, n2v_conf.p, n2v_conf.q, rw_sequences);
    }
}

template<typename edge_data_t>
void test_node2vec()
{
    edge_id_t e_nums_arr[] = {100, 200, 400, 556, 888, 1000, 1200, 1500};
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
        test_node2vec<edge_data_t>(v_num, worker_number);
    }
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(Node2vec, Unbiased)
{
    test_node2vec<EmptyData>();
}

TEST(Node2vec, Biased)
{
    Node2vecConf n2v_conf;
    test_node2vec<real_t>();
}


GTEST_API_ int main(int argc, char *argv[])
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}
