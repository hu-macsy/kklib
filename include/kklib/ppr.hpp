#pragma once

#include <kklib/graph.hpp>
#include <kklib/static_comp.hpp>
#include <kklib/walk.hpp>

#include <mpi.h>

template <typename edge_data_t>
void ppr(WalkEngine<edge_data_t, EmptyData>* graph,
         walker_id_t walker_num,
         real_t terminate_prob,
         std::vector<VertexID>* start_vertices = nullptr,
         WalkConfig walk_conf = WalkConfig{})
{
    MPI_Barrier(MPI_COMM_WORLD);
    Timer timer;

    real_t extension_comp = 1 - terminate_prob;
    WalkerConfig<edge_data_t, EmptyData> walker_conf(walker_num, nullptr, nullptr,
                                                     [&](walker_id_t walker_id)
                                                     {
                                                         if (start_vertices == nullptr)
                                                         {
                                                             return (VertexID)walker_id % graph->get_vertex_num();
                                                         }
                                                         else
                                                         {
                                                             return (
                                                                 VertexID)(*start_vertices)[walker_id % start_vertices->size()];
                                                         }
                                                     });
    TransitionConfig<edge_data_t, EmptyData> tr_conf([&](Walker<EmptyData>& walker, VertexID current_v)
                                                     { return extension_comp; },
                                                     get_trivial_static_comp(graph));
    graph->random_walk(&walker_conf, &tr_conf, walk_conf);

#if !defined(UNIT_TEST) && !defined(NDEBUG)
    printf("total time %lfs\n", timer.duration());
#endif
}