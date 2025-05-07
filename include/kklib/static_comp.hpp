#pragma once

#include <kklib/graph_data_structure.hpp>
#include <kklib/walk.hpp>

#include <functional>

template <typename T> std::function<real_t(VertexID, AdjUnit<T>*)> get_trivial_static_comp()
{
    printf("[error] Undefined trivial static component\n");
    exit(1);
}

template <typename walker_state_t>
std::function<real_t(VertexID, AdjUnit<EmptyData>*)> get_trivial_static_comp(WalkEngine<EmptyData, walker_state_t>* graph)
{
    return nullptr;
}

template <typename walker_state_t>
std::function<real_t(VertexID, AdjUnit<real_t>*)> get_trivial_static_comp(WalkEngine<real_t, walker_state_t>* graph)
{
    auto static_comp = [&](VertexID v, AdjUnit<real_t>* edge) { return edge->data; };
    return static_comp;
}
