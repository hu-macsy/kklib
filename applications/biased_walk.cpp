#include "option_helper.hpp"

#include <kklib/mpi_helper.hpp>
#include <kklib/storage.hpp>
#include <kklib/walk.hpp>

//! Quick use:
//! ./biased_walk -g ./karate.data -v 34 -w 1 -o ./out/walks.txt -l 5
int main(int argc, char** argv)
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);

    TruncatedRandomWalkOptionHelper opt;
    opt.parse(argc, argv);
    kklib::LoadedGraphData graph_data = kklib::load_graph<real_t>(opt.v_num, opt.graph_path.c_str());

    WalkEngine<real_t, EmptyData> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };

    WalkConfig walk_conf;
    if (!opt.output_path.empty())
    {
        walk_conf.set_output_file(opt.output_path.c_str());
    }
    if (opt.set_rate)
    {
        walk_conf.set_walk_rate(opt.rate);
    }
    WalkerConfig<real_t, EmptyData> walker_conf(opt.walker_num);
    auto extension_comp = [&](Walker<EmptyData>& walker, VertexID current_v)
    {
        return walker.step >= opt.walk_length ? 0.0 : 1.0; /*walk opt.walk_length steps then terminate*/
    };
    auto static_comp = [&](VertexID v, AdjUnit<real_t>* edge)
    {
        return edge->data; /*edge->data is a real number denoting edge weight*/
    };
    TransitionConfig<real_t, EmptyData> tr_conf(extension_comp, static_comp);
    graph.random_walk(&walker_conf, &tr_conf, walk_conf);
    return 0;
}
