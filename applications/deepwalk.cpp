#include "option_helper.hpp"

#include <kklib/static_comp.hpp>
#include <kklib/storage.hpp>
#include <kklib/walk.hpp>

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

#ifndef UNIT_TEST
    printf("total time %lfs\n", timer.duration());
#endif
}

template <typename edge_data_t>
void run(WalkEngine<edge_data_t, EmptyData>* graph, STruncatedRandomWalkOptionHelper* opt)
{
    WalkConfig walk_conf;
    if (!opt->output_path.empty())
    {
        walk_conf.set_output_file(opt->output_path.c_str());
    }
    if (opt->set_rate)
    {
        walk_conf.set_walk_rate(opt->rate);
    }
    deepwalk(graph, opt->walker_num, opt->walk_length, walk_conf);
}

int main(int argc, char** argv)
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);

    STruncatedRandomWalkOptionHelper opt;
    opt.parse(argc, argv);

    if (opt.static_comp.compare("weighted") == 0)
    {
        kklib::LoadedGraphData graph_data = kklib::load_graph<real_t>(opt.v_num, opt.graph_path.c_str());
        WalkEngine<real_t, EmptyData> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
        run(&graph, &opt);
    }
    else if (opt.static_comp.compare("unweighted") == 0)
    {
        kklib::LoadedGraphData graph_data = kklib::load_graph<EmptyData>(opt.v_num, opt.graph_path.c_str());
        WalkEngine<EmptyData, EmptyData> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
        run(&graph, &opt);
    }
    else
    {
        exit(1);
    }
    return 0;
}
