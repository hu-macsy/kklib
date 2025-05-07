#pragma once

#include <kklib/graph.hpp>
#include <kklib/metascheme.hpp>
#include <kklib/walk.hpp>

#include <mpi.h>

template <typename walker_state_t>
std::function<real_t(VertexID, AdjUnit<UnweightedMetaData>*)>
get_metapath_static_comp(WalkEngine<UnweightedMetaData, walker_state_t>* graph)
{
    return nullptr;
}

template <typename walker_state_t>
std::function<real_t(VertexID, AdjUnit<WeightedMetaData>*)>
get_metapath_static_comp(WalkEngine<WeightedMetaData, walker_state_t>* graph)
{
    auto static_comp = [&](VertexID v, AdjUnit<WeightedMetaData>* edge) { return edge->data.weight; };
    return static_comp;
}

template <typename edge_data_t>
void metapath(WalkEngine<edge_data_t, MetapathState>* graph,
              std::vector<std::vector<std::vector<bool>>> schemes,
              walker_id_t walker_num,
              step_t walk_length,
              WalkConfig walk_conf = WalkConfig{})
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    auto scheme_masks = get_scheme_mask(schemes);
    scheme_mask_t* vertex_masks = graph->template alloc_vertex_array<scheme_mask_t>();
    graph->template process_vertices<VertexID>(
        [&](VertexID v_id)
        {
            vertex_masks[v_id] = 0;
            for (auto p = graph->csr->adj_lists[v_id].begin; p < graph->csr->adj_lists[v_id].end; p++)
            {
                vertex_masks[v_id] |= (1 << p->data.get_meta());
            }
            return 0;
        });
    WalkerConfig<edge_data_t, MetapathState> walker_conf(
        walker_num,
        [&](Walker<MetapathState>& walker, VertexID start_vertex)
        {
            walker.data.scheme_id = graph->thread_local_random_engine_int()(0, schemes.size() - 1);
            walker.data.state = 0;
        },
        [&](Walker<MetapathState>& walker, VertexID current_v, AdjUnit<edge_data_t>* edge)
        { walker.data.state = (walker.data.state + 1) % schemes[walker.data.scheme_id].size(); });
    TransitionConfig<edge_data_t, MetapathState> tr_conf(
        [&](Walker<MetapathState>& walker, VertexID current_v)
        {
            return (walker.step >= walk_length ||
                    !(vertex_masks[current_v] & scheme_masks[walker.data.scheme_id][walker.data.state])) ?
                       0.0 :
                       1.0;
        },
        get_metapath_static_comp(graph),
        [&](Walker<MetapathState>& walker, VertexID current_v, AdjUnit<edge_data_t>* edge)
        {
            if (schemes[walker.data.scheme_id][walker.data.state][edge->data.get_meta()])
            {
                return 1.0;
            }
            else
            {
                return 0.0;
            }
        },
        [&](VertexID v_id, AdjList<edge_data_t>* adj_lists) { return 1.0; });
    graph->random_walk(&walker_conf, &tr_conf, walk_conf);
    graph->dealloc_vertex_array(vertex_masks);

#if !defined(UNIT_TEST) && !defined(NDEBUG)
    printf("total time %lfs\n", timer.duration());
#endif
}