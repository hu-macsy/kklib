#include "test.hpp"
#include "test_walk.hpp"

#include <kklib/graph.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/storage.hpp>
#include <kklib/util.hpp>
#include <kklib/walk.hpp>

#include <gtest/gtest.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>


template <typename edge_data_t>
void deepwalk(WalkEngine<edge_data_t, EmptyData>* graph, walker_id_t walker_num, step_t walk_length, WalkConfig walk_conf = WalkConfig{})
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    WalkerConfig<edge_data_t, EmptyData> walker_conf(walker_num);
    auto extension_comp = [&](Walker<EmptyData>& walker, VertexID current_v)
    { return walker.step >= walk_length ? 0.0 : 1.0; };
    auto static_comp = get_trivial_static_comp(graph);
    TransitionConfig<edge_data_t, EmptyData> tr_conf(extension_comp, static_comp);
    graph->random_walk(&walker_conf, &tr_conf, walk_conf);
    printf("total time %lfs\n", timer.duration());
}

template <typename edge_data_t> void test_deepwalk(VertexID v_num, int worker_number)
{
    kklib::LoadedGraphData graph_data = kklib::load_graph<edge_data_t>(v_num, test_data_file);
    WalkEngine<edge_data_t, EmptyData> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };

    graph.set_concurrency(worker_number);

    step_t walk_length = 60 + rand() % 20;
    MPI_Bcast(&walk_length, 1, kklib::deduce_mpi_data_type<step_t>(), 0, MPI_COMM_WORLD);
    walker_id_t walker_num = graph.get_vertex_num() * 50 + graph.get_edge_num() * 10 + rand() % 100;
    MPI_Bcast(&walker_num, 1, kklib::deduce_mpi_data_type<walker_id_t>(), 0, MPI_COMM_WORLD);

    deepwalk(&graph, walker_num, walk_length);
    std::vector<std::vector<VertexID>> rw_sequences;
    graph.collect_walk_sequence(rw_sequences, walker_num);

    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        check_static_random_walk(v_num, std_edges, std_edge_num, rw_sequences);
    }
}

template<typename edge_data_t>
void test_deepwalk()
{
    size_t e_nums_arr[] = {1000, 2000, 4000, 5556, 8888, 10000, 20000, 100000};
    VertexID v_num = 1000 + rand() % 500;
    std::vector<size_t> e_nums(e_nums_arr, e_nums_arr + 8);
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
        test_deepwalk<edge_data_t>(v_num, worker_number);
    }
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(DeepWalk, UnbiasedDeepWalk)
{
    test_deepwalk<EmptyData>();
}

TEST(DeepWalk, BiasedDeepWalk)
{
    test_deepwalk<real_t>();
}

GTEST_API_ int main(int argc, char *argv[])
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}
