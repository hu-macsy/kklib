#include "option_helper.hpp"

#include <kklib/storage.hpp>
#include <kklib/walk.hpp>

struct WalkState
{
    VertexID last_vertex;
};

int main(int argc, char** argv)
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);

    TruncatedRandomWalkOptionHelper opt;
    opt.parse(argc, argv);

    kklib::LoadedGraphData graph_data = kklib::load_graph<real_t>(opt.v_num, opt.graph_path.c_str());
    WalkEngine<real_t, WalkState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };

    WalkConfig walk_conf;
    if (!opt.output_path.empty())
    {
        walk_conf.set_output_file(opt.output_path.c_str());
    }
    if (opt.set_rate)
    {
        walk_conf.set_walk_rate(opt.rate);
    }

    auto init_walker_func = [&](Walker<WalkState>& walker, VertexID start_vertex)
    {
        /*At first, the last vertex is not defined*/
        walker.data.last_vertex = UINT_MAX;
    };
    auto update_walker_func = [&](Walker<WalkState>& walker, VertexID current_v, AdjUnit<real_t>* edge)
    { walker.data.last_vertex = current_v; };
    WalkerConfig<real_t, WalkState> walker_conf(34, init_walker_func, update_walker_func);

    auto extension_comp = [&](Walker<WalkState>& walker, VertexID current_v)
    {
        return walker.step >= opt.walk_length ? 0.0 : 1.0; /*walk walk_length steps then terminate*/
    };
    auto static_comp = [&](VertexID v, AdjUnit<real_t>* edge)
    {
        return edge->data; /*edge->data is a real number denoting edge weight*/
    };
    auto dynamic_comp = [&](Walker<WalkState>& walker, VertexID current_v, AdjUnit<real_t>* edge)
    {
        if (walker.step == 0)
        {
            /*No return edge for the first step*/
            return 1.0;
        }
        else if (edge->neighbour == walker.data.last_vertex)
        {
            /*if return edge, double the un-normalized transition probability*/
            return 2.0;
        }
        else
        {
            /*if not return edge*/
            return 1.0;
        }
    };
    auto dynamic_comp_upperbound = [&](VertexID v_id, AdjList<real_t>* adj_lists) { return 2.0; };

    TransitionConfig<real_t, WalkState> tr_conf(extension_comp, static_comp, dynamic_comp, dynamic_comp_upperbound);
    graph.random_walk(&walker_conf, &tr_conf, walk_conf);

    return 0;
}
