#include "option_helper.hpp"

#include <kklib/metapath.hpp>
#include <kklib/metascheme.hpp>
#include <kklib/storage.hpp>
#include <kklib/walk.hpp>

class MetapathOptionHelper : public STruncatedRandomWalkOptionHelper
{
private:
    args::ValueFlag<std::string> schemes_path_flag;

public:
    std::string schemes_path;
    MetapathOptionHelper()
    : schemes_path_flag(parser, "schemes", "schemes file path", { "schemes" })
    {
    }
    virtual void parse(int argc, char** argv)
    {
        STruncatedRandomWalkOptionHelper::parse(argc, argv);

        assert(schemes_path_flag);
        schemes_path = args::get(schemes_path_flag);
    }
};

template <typename edge_data_t> void run(WalkEngine<edge_data_t, MetapathState>* graph, MetapathOptionHelper* opt)
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
    metapath(graph, read_metapath_schemes(opt->schemes_path.c_str()), opt->walker_num, opt->walk_length, walk_conf);
}

int main(int argc, char** argv)
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);

    MetapathOptionHelper opt;
    opt.parse(argc, argv);

    if (opt.static_comp.compare("weighted") == 0)
    {
        kklib::LoadedGraphData graph_data = kklib::load_graph<WeightedMetaData>(opt.v_num, opt.graph_path.c_str());
        WalkEngine<WeightedMetaData, MetapathState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };

        run(&graph, &opt);
    }
    else if (opt.static_comp.compare("unweighted") == 0)
    {
        kklib::LoadedGraphData graph_data =
            kklib::load_graph<UnweightedMetaData>(opt.v_num, opt.graph_path.c_str(), opt.make_undirected);
        WalkEngine<UnweightedMetaData, MetapathState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };

        run(&graph, &opt);
    }
    else
    {
        exit(1);
    }
    return 0;
}
