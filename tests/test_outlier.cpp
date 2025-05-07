#include <gtest/gtest.h>

#include "test.hpp"
#include "test_walk.hpp"

#include <kklib/graph.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/storage.hpp>
#include <kklib/util.hpp>
#include <kklib/walk.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <queue>
#include <type_traits>
#include <utility>
#include <vector>

typedef int tag_t;
const tag_t tag_num = 4;

struct TagWalkState
{
    tag_t tag;
    union
    {
        VertexID previous_vertex;
        tag_t previous_vertex_tag;
    };
};

struct WeightedTag 
{
    tag_t tag;
    real_t weight;
    real_t get_weight()
    {
        return weight;
    }
};

struct UnweightedTag
{
    tag_t tag;
    real_t get_weight()
    {
        return 1.0;
    }
};

struct TagWalkConf
{
    walker_id_t walker_num;
    step_t walk_length;
    real_t outlier_amplifier;
    tag_t tag_num;
};

template<typename edge_data_t>
void tagwalk(WalkEngine<edge_data_t, TagWalkState>* graph, TagWalkConf conf, tag_t* vertex_tag, int order)
{
    struct TagSet
    {
        AdjUnit<edge_data_t> *begin[tag_num];
        AdjUnit<edge_data_t> *end[tag_num];
    };
    TagSet* adj_tags = graph->template alloc_vertex_array<TagSet>();
    graph->template process_vertices<VertexID>(
        [&](VertexID v_i)
        {
            std::sort(graph->csr->adj_lists[v_i].begin, graph->csr->adj_lists[v_i].end,
                      [](const AdjUnit<edge_data_t> a, const AdjUnit<edge_data_t> b) { return a.data.tag < b.data.tag; });
            for (tag_t c_i = 0; c_i < tag_num; c_i++)
            {
                adj_tags[v_i].begin[c_i] = nullptr;
                adj_tags[v_i].end[c_i] = nullptr;
            }
            for (auto p = graph->csr->adj_lists[v_i].begin; p < graph->csr->adj_lists[v_i].end; p++)
            {
                tag_t c = p->data.tag;
                if (adj_tags[v_i].begin[c] == nullptr)
                {
                    adj_tags[v_i].begin[c] = p;
                }
                adj_tags[v_i].end[c] = p + 1;
            }
            return 0;
        });

    WalkerConfig<edge_data_t, TagWalkState> walker_conf(
        conf.walker_num,
        [&](Walker<TagWalkState>& walker, VertexID start_vertex) { walker.data.tag = walker.id % tag_num; },
        [&](Walker<TagWalkState>& walker, VertexID current_v, AdjUnit<edge_data_t>* edge)
        {
            walker.data.tag = (walker.data.tag + 1) % tag_num;
            walker.data.previous_vertex = current_v;
        });
    auto extension_comp = [&](Walker<TagWalkState>& walker, VertexID current_v)
    { return walker.step >= conf.walk_length ? 0.0 : 1.0; };
    auto static_comp = get_trivial_static_comp<edge_data_t>();
    auto upper_bound_func = [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return (real_t)tag_num; };
    auto lower_bound_func = [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return (real_t)1; };
    auto outlier_upperbound_func = [&](Walker<TagWalkState>& walker, VertexID vertex, AdjList<edge_data_t>* adj_list,
                                       real_t& prob_upperbound, VertexID& num_upperbound)
    {
        if (walker.step != 0)
        {
            auto begin = adj_tags[vertex].begin[walker.data.tag];
            auto end = adj_tags[vertex].end[walker.data.tag];
            prob_upperbound = 0;
            num_upperbound = end - begin;
            for (auto p = begin; p < end; p++)
            {
                prob_upperbound = std::max(prob_upperbound, tag_num * (conf.outlier_amplifier - 1) * p->data.get_weight());
            }
        }
        else
        {
            prob_upperbound = 0;
            num_upperbound = 0;
        }
    };
    auto outlier_search_func = [&](Walker<TagWalkState>& walker, VertexID vertex, AdjList<edge_data_t>* adj_list, VertexID outlier_idx)
    { return adj_tags[vertex].begin[walker.data.tag] + outlier_idx; };
    if (order == 1) 
    {
        TransitionConfig<edge_data_t, TagWalkState> tr_conf(
            extension_comp, static_comp,
            [&](Walker<TagWalkState>& walker, VertexID vertex, AdjUnit<edge_data_t>* edge)
            {
                if (walker.step == 0)
                {
                    return (real_t)tag_num;
                }
                else
                {
                    real_t temp = (real_t)1 + (vertex_tag[walker.data.previous_vertex] + walker.data.tag + edge->data.tag) % tag_num;
                    if (walker.data.tag == edge->data.tag)
                    {
                        temp *= conf.outlier_amplifier;
                    }
                    return temp;
                }
            },
            upper_bound_func, lower_bound_func, outlier_upperbound_func, outlier_search_func);
        graph->random_walk(&walker_conf, &tr_conf);
    } else
    {
        SecondOrderTransitionConfig<edge_data_t, TagWalkState, EmptyData, tag_t> tr_conf(
            extension_comp, static_comp,
            [&](Walker<TagWalkState>& walker, walker_id_t walker_idx, VertexID current_v, AdjUnit<edge_data_t>* edge)
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
                stateResponse<tag_t> response;
                response.walker_idx = query.walker_idx;
                response.data = vertex_tag[vtx];
                graph->emit(query.src_v, response);
            },
            [&](Walker<TagWalkState>& walker, stateResponse<tag_t>& response, VertexID current_v, AdjUnit<edge_data_t>* edge)
            {
                if (walker.step == 0)
                {
                    return (real_t)tag_num;
                }
                else
                {
                    real_t temp = (real_t)1 + (response.data + walker.data.tag + edge->data.tag) % tag_num;
                    if (walker.data.tag == edge->data.tag)
                    {
                        temp *= conf.outlier_amplifier;
                    }
                    return temp;
                }
            },
            upper_bound_func, lower_bound_func, outlier_upperbound_func, outlier_search_func);
        graph->random_walk(&walker_conf, &tr_conf);
    }

    graph->dealloc_vertex_array(adj_tags);
}

template <typename edge_data_t>
void check_tagwalk_random_walk(VertexID v_num,
                               Edge<edge_data_t>* edges,
                               edge_id_t e_num,
                               TagWalkConf conf,
                               tag_t* vertex_tag,
                               std::vector<std::vector<VertexID>>& seq)
{
    const size_t max_state_num = tag_num * tag_num;
    auto get_state_id = [&] (tag_t walker_tag, tag_t previous_vertex_tag)
    {
        return (size_t)walker_tag * tag_num + previous_vertex_tag;
    };
    std::vector<std::vector<Edge<edge_data_t> > > graph(v_num);
    for (edge_id_t e_i = 0; e_i < e_num; e_i++)
    {
        graph[edges[e_i].src].push_back(edges[e_i]);
    }
    for (auto &adj : graph)
    {
        std::sort(adj.begin(), adj.end(), [](const Edge<edge_data_t> a, const Edge<edge_data_t> b){return a.dst < b.dst;});
    }
    auto get_edge_idx = [&](VertexID src, VertexID dst)
    {
        Edge<edge_data_t> target;
        target.dst = dst;
        auto p = std::lower_bound(graph[src].begin(), graph[src].end(), target, [](const Edge<edge_data_t> &a, const Edge<edge_data_t> &b) { return a.dst < b.dst; });
        assert(p != graph[src].end());
        return p - graph[src].begin();
    };
    for (auto &s : seq)
    {
        VertexID start = s[0];
        assert((graph[start].size() == 0 && s.size() == 1) || (s.size() == conf.walk_length + 1));
    }
    std::vector<std::vector<bool> > vis(v_num);
    for (auto &x : vis)
    {
        x.resize(max_state_num, false);
    }
    struct QueueItem
    {
        tag_t wk_tag;
        VertexID vertex;
    };
    for (walker_id_t w_i = 0; w_i < seq.size(); w_i ++)
    {
        std::queue<QueueItem> q;
        auto expand_func = [&] (QueueItem current)
        {
            for (auto edge : graph[current.vertex])
            {
                QueueItem next;
                next.wk_tag = (current.wk_tag + 1) % tag_num;
                next.vertex = edge.dst;
                tag_t pv_tag = vertex_tag[current.vertex];
                size_t state_id = get_state_id(next.wk_tag, pv_tag);
                if (!vis[next.vertex][state_id])
                {
                    vis[next.vertex][state_id] = true;
                    q.push(next);
                }
            }

        };
        QueueItem start;
        start.vertex = seq[w_i][0];
        start.wk_tag = w_i % tag_num;
        expand_func(start);
        while (q.empty() == false)
        {
            QueueItem current = q.front();
            q.pop();
            expand_func(current);
        }
    }
    std::vector<std::vector<std::vector<double> > > std_trans_mat(v_num);
    for (VertexID v_i = 0; v_i < v_num; v_i++)
    {
        std_trans_mat[v_i].resize(max_state_num);
        for (tag_t walker_tag = 0; walker_tag < tag_num; walker_tag++)
        {
            for (tag_t pv_tag = 0; pv_tag < tag_num; pv_tag++)
            {
                size_t s_i = get_state_id(walker_tag, pv_tag);
                auto &dist = std_trans_mat[v_i][s_i];
                dist.resize(graph[v_i].size(), 0.0);
                if (!vis[v_i][s_i])
                {
                    continue;
                }
                for (size_t e_i = 0; e_i < graph[v_i].size(); e_i++)
                {
                    auto &edge = graph[v_i][e_i];
                    double val = (real_t)1 + (pv_tag + walker_tag + edge.data.tag) % tag_num;
                    if (walker_tag == edge.data.tag)
                    {
                        val *= conf.outlier_amplifier;
                    }
                    val *= edge.data.get_weight();
                    dist[e_i] = val;
                }
            }
        }
    }
    std::vector<std::vector<std::vector<double> > > real_trans_mat(v_num);
    for (VertexID v_i = 0; v_i < v_num; v_i++)
    {
        real_trans_mat[v_i].resize(max_state_num);
        for (size_t s_i = 0; s_i < max_state_num; s_i++)
        {
            real_trans_mat[v_i][s_i].resize(graph[v_i].size(), 0.0);
        }
    }
    for (walker_id_t w_i = 0; w_i < seq.size(); w_i++)
    {
        for (step_t s_i = 1; s_i + 1 < seq[w_i].size(); s_i++)
        {
            VertexID current_v = seq[w_i][s_i];
            size_t state_id = get_state_id((w_i + s_i) % tag_num, vertex_tag[seq[w_i][s_i - 1]]);
            size_t edge_idx = get_edge_idx(seq[w_i][s_i], seq[w_i][s_i + 1]);
            real_trans_mat[current_v][state_id][edge_idx] += 1.0;
        }
    }
    auto get_flat_mat = [] (std::vector<std::vector<std::vector<double> > > &three_d_mat)
    {
        std::vector<std::vector<double> > two_d_mat;
        for (auto &x : three_d_mat)
        {
            for (auto &y : x)
            {
                two_d_mat.push_back(y);
            }
        }
        return two_d_mat;
    };
    auto std_flat_mat = get_flat_mat(std_trans_mat);
    auto real_flat_mat = get_flat_mat(real_trans_mat);
    mat_normalization(std_flat_mat);
    mat_normalization(real_flat_mat);
    cmp_trans_matrix(real_flat_mat, std_flat_mat);
}

template <typename edge_data_t> void test_outlier(VertexID v_num, int worker_number, int order)
{
    kklib::LoadedGraphData graph_data = kklib::load_graph<edge_data_t>(v_num, test_data_file);
    WalkEngine<edge_data_t, TagWalkState> graph{ graph_data.v_num_param, graph_data.read_edges, graph_data.read_e_num };
    graph.set_concurrency(worker_number);

    TagWalkConf conf;
    conf.walk_length = 80 + rand() % 20;
    conf.walker_num = graph.get_vertex_num() * 500 + graph.get_edge_num() * 500 + rand() % 100;
    conf.tag_num = 3 + rand() % 5;
    conf.outlier_amplifier = 1.0 + 5.0 / (rand() % 10 + 1);
    MPI_Bcast(&conf, sizeof(conf), MPI_CHAR, 0, MPI_COMM_WORLD);

    tag_t* vertex_tag = graph.template alloc_vertex_array<tag_t>();
    if (get_mpi_rank() == 0)
    {
        kklib::RandomEngine<tag_t> random_engine;
        for (VertexID v_i = 0; v_i < v_num; v_i++)
        {
            vertex_tag[v_i] = random_engine(0, tag_num - 1);
        }
    }
    MPI_Bcast(vertex_tag, v_num, kklib::deduce_mpi_data_type<tag_t>(), 0, MPI_COMM_WORLD);

    tagwalk(&graph, conf, vertex_tag, order);

    std::vector<std::vector<VertexID>> rw_sequences;
    graph.collect_walk_sequence(rw_sequences, conf.walker_num);

    if (get_mpi_rank() == 0)
    {
        Edge<edge_data_t> *std_edges;
        edge_id_t std_edge_num;
        read_graph(test_data_file, 0, 1, std_edges, std_edge_num);
        check_tagwalk_random_walk(v_num, std_edges, std_edge_num, conf, vertex_tag, rw_sequences);
    }
    graph.dealloc_vertex_array(vertex_tag);
}

template<typename T>
std::function<void(T&)> get_tag_gen_func()
{
    printf("[error] Undefined edge data generator\n");
    exit(1);
}

template<>
std::function<void(UnweightedTag&)> get_tag_gen_func<UnweightedTag>()
{
    auto func = [] (UnweightedTag &data)
    {
        data.tag = rand() % tag_num;
    };
    return func;
}

template<>
std::function<void(WeightedTag&)> get_tag_gen_func<WeightedTag>()
{
    auto func = [] (WeightedTag &data)
    {
        data.tag = rand() % tag_num;
        kklib::RandomEngine<real_t> random_engine;
        constexpr uint32_t range = 5;
        data.weight = random_engine(0., range);
    };
    return func;
}

template<typename edge_data_t>
void test_outlier(int order)
{
    edge_id_t e_nums_arr[] = {100, 200, 300, 400, 500, 600};
    VertexID v_num = 50 + rand() % 25;
    std::vector<edge_id_t> e_nums(e_nums_arr, e_nums_arr + 6);
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
            gen_undirected_graph_file<edge_data_t>(v_num, e_num, get_tag_gen_func<edge_data_t>());
        }
        MPI_Barrier(MPI_COMM_WORLD);
        int worker_number = rand() % 8 + 1;
        MPI_Bcast(&worker_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
        test_outlier<edge_data_t>(v_num, worker_number, order);
    }
    if (get_mpi_rank() == 0)
    {
        rm_test_graph_temp_file();
    }
}

TEST(Outlier, UnbiasedFirstOrder)
{
    test_outlier<UnweightedTag>(1);
}

TEST(Outlier, BiasedFirstOrder)
{
    test_outlier<WeightedTag>(1);
}

TEST(Outlier, UnbiasedSecondOrder)
{
    test_outlier<UnweightedTag>(2);
}

TEST(Outlier, BiasedSecondOrder)
{
    test_outlier<WeightedTag>(2);
}


GTEST_API_ int main(int argc, char *argv[])
{
    kklib::MPI_Instance mpi_instance(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    mute_nonroot_gtest_events();
    int result = RUN_ALL_TESTS();
    return result;
}
