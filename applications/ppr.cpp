#include "option_helper.hpp"

#include <kklib/ppr.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/storage.hpp>
#include <kklib/walk.hpp>

class PPROptionHelper : public SRandomWalkOptionHelper
{
private:
    args::ValueFlag<real_t> terminate_prob_flag;

public:
    real_t terminate_prob;
    PPROptionHelper()
    : terminate_prob_flag(parser, "terminate", "terminate probabiility", { 't' })
    {
    }
    virtual void parse(int argc, char** argv)
    {
        SRandomWalkOptionHelper::parse(argc, argv);

        assert(terminate_prob_flag);
        terminate_prob = args::get(terminate_prob_flag);
    }
};

template <typename edge_data_t> void run(WalkEngine<edge_data_t, EmptyData>* graph, PPROptionHelper* opt)
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
    ppr(graph, opt->walker_num, opt->terminate_prob, nullptr, walk_conf);
}

int main(int argc, char** argv)
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);

    PPROptionHelper opt;
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
