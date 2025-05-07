#pragma once

#include <kklib/type.hpp>

#include <algorithm>

struct EmptyData
{
};

template <typename edge_data_t> struct Edge
{
    VertexID src;
    VertexID dst;
    edge_data_t data;

    Edge() {}

    Edge(VertexID _src, VertexID _dst, edge_data_t _data)
    : src(_src)
    , dst(_dst)
    , data(_data)
    {
    }

    bool friend operator==(const Edge<edge_data_t>& a, const Edge<edge_data_t>& b)
    {
        return (a.src == b.src && a.dst == b.dst && a.data == b.data);
    }

    void transpose() { std::swap(src, dst); }
};

template <> struct Edge<EmptyData>
{
    VertexID src;

    //! Making this a union must have the effect, that either dst or data is an
    //! active member. Does that mean that data must include the destination or
    //! how is that handled? Since in this case EmptyData is an empty struct it
    //! must be of VertexID. VertexID is 32bit. Therefore, the minimum size of
    //! this union is 32bit.
    //!
    //! So this union used to hide the >1 byte allocation for EmptyData?!
    union
    {
        VertexID dst;
        EmptyData data;
    };

    Edge() {}

    Edge(VertexID _src, VertexID _dst)
    : src(_src)
    , dst(_dst)
    {
    }

    bool friend operator==(const Edge<EmptyData>& a, const Edge<EmptyData>& b)
    {
        return (a.src == b.src && a.dst == b.dst);
    }

    void transpose() { std::swap(src, dst); }
};

template <typename edge_data_t> struct AdjUnit
{
    VertexID neighbour;
    edge_data_t data;
};

template <> struct AdjUnit<EmptyData>
{
    union
    {
        VertexID neighbour;
        EmptyData data;
    };
};

template <typename edge_data_t> struct AdjList
{
    AdjUnit<edge_data_t>* begin;
    AdjUnit<edge_data_t>* end;
    void init()
    {
        begin = nullptr;
        end = nullptr;
    }
};

// comprised column row
template <typename edge_data_t> struct EdgeContainer
{
    AdjList<edge_data_t>* adj_lists;
    AdjUnit<edge_data_t>* adj_units;
    EdgeContainer()
    : adj_lists(nullptr)
    , adj_units(nullptr)
    {
    }
    ~EdgeContainer()
    {
        if (adj_lists != nullptr)
        {
            delete[] adj_lists;
        }
        if (adj_units != nullptr)
        {
            delete[] adj_units;
        }
    }
};

template <typename EdgeDataT>
EdgeContainer<EdgeDataT>* build_edge_container(Edge<EdgeDataT> const* const edges,
                                               edge_id_t local_edge_num,
                                               VertexID* vertex_out_degree,
                                               VertexID const v_num,
                                               partition_id_t const local_partition_id,
                                               VertexID* vertex_partition_begin,
                                               VertexID* vertex_partition_end)
{
    EdgeContainer<EdgeDataT>* ec = new EdgeContainer<EdgeDataT>();
    ec->adj_lists = new AdjList<EdgeDataT>[v_num];
    ec->adj_units = new AdjUnit<EdgeDataT>[local_edge_num];

    edge_id_t chunk_edge_idx = 0;
    for (VertexID v_i = vertex_partition_begin[local_partition_id]; v_i < vertex_partition_end[local_partition_id]; v_i++)
    {
        ec->adj_lists[v_i].begin = ec->adj_units + chunk_edge_idx;
        ec->adj_lists[v_i].end = ec->adj_lists[v_i].begin;
        chunk_edge_idx += vertex_out_degree[v_i];
    }

    for (edge_id_t e_i = 0; e_i < local_edge_num; e_i++)
    {
        auto e = edges[e_i];
        auto ep = ec->adj_lists[e.src].end++;
        ep->neighbour = e.dst;
        if (!std::is_same<EdgeDataT, EmptyData>::value)
        {
            ep->data = e.data;
        }
    }

    return ec;
}

namespace kklib
{

using EdgeNoData = Edge<EmptyData>;

} // namespace kklib