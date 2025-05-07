#include "option_helper.hpp"

#include <kklib/node2vec.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/storage.hpp>
#include <kklib/walk.hpp>

class Node2vecOptionHelper : public STruncatedRandomWalkOptionHelper
{
private:
    args::ValueFlag<float> p_flag;
    args::ValueFlag<float> q_flag;

public:
    float p;
    float q;
    Node2vecOptionHelper()
    : p_flag(parser, "p", "hyperparameter p", { 'p' })
    , q_flag(parser, "q", "hyperparameter q", { 'q' })
    {
    }
    virtual void parse(int argc, char** argv)
    {
        STruncatedRandomWalkOptionHelper::parse(argc, argv);

        assert(p_flag);
        p = args::get(p_flag);

        assert(q_flag);
        q = args::get(q_flag);
    }
    virtual Node2vecConf get_n2v_conf()
    {
        Node2vecConf n2v_conf;
        n2v_conf.p = p;
        n2v_conf.q = q;
        n2v_conf.walker_num = walker_num;
        n2v_conf.walk_length = walk_length;
        return n2v_conf;
    }
};

template <typename edge_data_t> void run(WalkEngine<edge_data_t, Node2vecState>* graph, Node2vecOptionHelper* opt)
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
    Node2vecConf n2v_conf = opt->get_n2v_conf();
    node2vec(graph, n2v_conf, walk_conf);
}

int main(int argc, char** argv)
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);

    Node2vecOptionHelper opt;
    opt.parse(argc, argv);

    if (opt.static_comp.compare("weighted") == 0)
    {
        kklib::LoadedGraphData graph_data = kklib::load_graph<real_t>(opt.v_num, opt.graph_path.c_str());
        WalkEngine<real_t, Node2vecState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
        run(&graph, &opt);
    }
    else if (opt.static_comp.compare("unweighted") == 0)
    {
        kklib::LoadedGraphData graph_data = kklib::load_graph<EmptyData>(opt.v_num, opt.graph_path.c_str());
        WalkEngine<EmptyData, Node2vecState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
        run(&graph, &opt);
    }
    else
    {
        exit(1);
    }
    return 0;
}
