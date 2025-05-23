#pragma once

#include <kklib/graph_data_structure.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/type.hpp>
#include <kklib/walk.hpp>

#include <functional>

struct Node2vecState
{
    VertexID previous_vertex;
};

struct Node2vecConf
{
    real_t p;
    real_t q;
    walker_id_t walker_num;
    step_t walk_length;
};

// For the outlier upperbound function and outlier
// search function, assume no duplicate edge exists.
// If there are duplicate edges, just use binary search to find
// how many they are or just giva an upper bound number.
struct Node2vecOutlierFuncWrapper
{
    static real_t overflow_prob;
    static std::function<void(Walker<Node2vecState>&, VertexID, AdjList<EmptyData>*, real_t&, VertexID&)>
    get_node2vec_outlier_upperbound_func(WalkEngine<EmptyData, Node2vecState>* graph)
    {
        // printf("overflow prob %f\n", overflow_prob);
        if (overflow_prob > 0)
        {
            auto outlier_upperbound_func = [&](Walker<Node2vecState>& walker, VertexID vertex, AdjList<EmptyData>* adj_list,
                                               real_t& prob_upperbound, VertexID& num_upperbound)
            {
                prob_upperbound = overflow_prob;
                num_upperbound = 1;
            };
            return outlier_upperbound_func;
        }
        else
        {
            return nullptr;
        }
    }

    static std::function<void(Walker<Node2vecState>&, VertexID, AdjList<real_t>*, real_t&, VertexID&)>
    get_node2vec_outlier_upperbound_func(WalkEngine<real_t, Node2vecState>* graph)
    {
        if (overflow_prob > 0)
        {
            auto outlier_upperbound_func = [&](Walker<Node2vecState>& walker, VertexID vertex, AdjList<real_t>* adj_list,
                                               real_t& prob_upperbound, VertexID& num_upperbound)
            {
                AdjUnit<real_t> target;
                target.neighbour = walker.data.previous_vertex;
                auto return_edge = std::lower_bound(adj_list->begin, adj_list->end, target,
                                                    [](const AdjUnit<real_t>& a, const AdjUnit<real_t>& b)
                                                    { return a.neighbour < b.neighbour; });
                prob_upperbound = return_edge->data * overflow_prob;
                num_upperbound = 1;
            };
            return outlier_upperbound_func;
        }
        else
        {
            return nullptr;
        }
    }

    template <typename edge_data_t>
    static std::function<AdjUnit<edge_data_t>*(Walker<Node2vecState>&, VertexID, AdjList<edge_data_t>*, VertexID)>
    get_node2vec_outlier_search_func(WalkEngine<edge_data_t, Node2vecState>* graph)
    {
        if (overflow_prob > 0)
        {
            auto outlier_search_func = [&](Walker<Node2vecState>& walker, VertexID vertex,
                                           AdjList<edge_data_t>* adj_list, VertexID outlier_idx)
            {
                AdjUnit<edge_data_t> target;
                target.neighbour = walker.data.previous_vertex;
                auto return_edge = std::lower_bound(adj_list->begin, adj_list->end, target,
                                                    [](const AdjUnit<edge_data_t>& a, const AdjUnit<edge_data_t>& b)
                                                    { return a.neighbour < b.neighbour; });
                return return_edge;
            };
            return outlier_search_func;
        }
        else
        {
            return nullptr;
        }
    }
};

real_t Node2vecOutlierFuncWrapper::overflow_prob = 0.0;

template <typename edge_data_t>
void node2vec(WalkEngine<edge_data_t, Node2vecState>* graph, Node2vecConf conf, WalkConfig walk_conf = WalkConfig{})
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    real_t p = conf.p;
    real_t q = conf.q;
    step_t walk_length = conf.walk_length;
    walker_id_t walker_num = conf.walker_num;

    VertexID local_vertex_begin = graph->get_local_vertex_begin();
    VertexID local_vertex_end = graph->get_local_vertex_end();

    graph->template process_vertices<VertexID>(
        [&](VertexID v_i)
        {
            std::sort(graph->csr->adj_lists[v_i].begin, graph->csr->adj_lists[v_i].end,
                      [](const AdjUnit<edge_data_t> a, const AdjUnit<edge_data_t> b)
                      { return a.neighbour < b.neighbour; });
            return 0;
        });
    real_t upperbound = std::max(1.0, 1.0 / q);
    real_t lowerbound = std::min(1.0 / p, std::min(1.0, 1.0 / q));
    WalkerConfig<edge_data_t, Node2vecState> walker_conf(walker_num, nullptr,
                                                         [&](Walker<Node2vecState>& walker, VertexID current_v, AdjUnit<edge_data_t>* edge)
                                                         { walker.data.previous_vertex = current_v; });

    Node2vecOutlierFuncWrapper::overflow_prob = 1.0 / p - upperbound;
    auto outlier_upperbound_func = Node2vecOutlierFuncWrapper::get_node2vec_outlier_upperbound_func(graph);
    auto outlier_search_func = Node2vecOutlierFuncWrapper::get_node2vec_outlier_search_func(graph);

    SecondOrderTransitionConfig<edge_data_t, Node2vecState, VertexID, bool> tr_conf(
        [&](Walker<Node2vecState>& walker, VertexID current_v) { return walker.step >= walk_length ? 0.0 : 1.0; },
        get_trivial_static_comp(graph),
        [&](Walker<Node2vecState>& walker, walker_id_t walker_idx, VertexID current_v, AdjUnit<edge_data_t>* edge)
        {
            if (walker.step != 0)
            {
                stateQuery<VertexID> query;
                query.src_v = current_v;
                query.walker_idx = walker_idx;
                query.data = edge->neighbour;
                graph->emit(walker.data.previous_vertex, query);
            }
        },
        [&](VertexID vtx, stateQuery<VertexID> query, AdjList<edge_data_t>* adj_list)
        {
            stateResponse<bool> response;
            response.walker_idx = query.walker_idx;
            AdjUnit<edge_data_t> target;
            target.neighbour = query.data;
            response.data = std::binary_search(adj_list->begin, adj_list->end, target,
                                               [](const AdjUnit<edge_data_t>& a, const AdjUnit<edge_data_t>& b)
                                               { return a.neighbour < b.neighbour; });
            graph->emit(query.src_v, response);
        },
        [&](Walker<Node2vecState>& walker, stateResponse<bool>& response, VertexID current_v, AdjUnit<edge_data_t>* edge)
        {
            if (walker.step == 0)
            {
                return upperbound;
            }
            else
            {
                if (walker.data.previous_vertex == edge->neighbour)
                {
                    return 1 / p;
                }
                else if (response.data)
                {
                    return (real_t)1;
                }
                else
                {
                    return 1 / q;
                }
            }
        },
        [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return upperbound; },
        [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return lowerbound; }, outlier_upperbound_func, outlier_search_func);
    graph->random_walk(&walker_conf, &tr_conf, walk_conf);

#if !defined(UNIT_TEST) && !defined(NDEBUG)
    printf("total time %lfs\n", timer.duration());
#endif
}