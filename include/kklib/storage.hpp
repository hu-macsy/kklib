#pragma once

#include <kklib/graph.hpp>
#include <kklib/mpi_helper.hpp>
#include <kklib/type.hpp>

#include <gdsb/graph_input.h>

#include <cassert>
#include <cstdio>
#include <fstream>
#include <type_traits>
#include <vector>

enum GraphFormat
{
    GF_Binary,
    GF_Edgelist
};

template <typename T> void read_graph(const char* fname, int partition_id, int partition_num, T*& edge, size_t& e_num)
{
    FILE *f = fopen(fname, "r");
    assert(f != NULL);
    fseek(f, 0, SEEK_END);
    size_t total_size = ftell(f);
    size_t total_e_num = total_size / sizeof(T);
    e_num = total_e_num / partition_num;
    if (partition_id == partition_num -1)
    {
        e_num += total_e_num % partition_num;
    }
    size_t f_offset = sizeof(T) * (total_e_num / partition_num) * partition_id; 

    edge = new T[e_num];
    fseek(f, f_offset, SEEK_SET);
    auto ret = fread(edge, sizeof(T), e_num, f);
    assert(ret == e_num);
    fclose(f);
}

namespace kklib
{
template <typename EdgeData> struct LoadedGraphData
{
    VertexID v_num_param;
    Edge<EdgeData>* read_edges;
    edge_id_t read_e_num;
};

template <typename Edge, typename EdgeData>
void read_graph(const char* fname, uint32_t partition_id, uint32_t partition_size, kklib::LoadedGraphData<EdgeData>& graph_data)
{
    size_t current_edge_id = 0;

    auto read_f_unweighted_32 = [&](std::ifstream& input)
    {
        gdsb::Edge32 tmp_edge;
        input.read((char*)&graph_data.read_edges[current_edge_id].src, sizeof(gdsb::Vertex32));
        input.read((char*)&graph_data.read_edges[current_edge_id].dst, sizeof(gdsb::Vertex32));

        current_edge_id++;
        return true;
    };

    auto read_f_weighted_32 = [&](std::ifstream& input)
    {
        gdsb::Edge32 tmp_edge;
        input.read((char*)&graph_data.read_edges[current_edge_id].src, sizeof(gdsb::Vertex32));
        input.read((char*)&graph_data.read_edges[current_edge_id].dst, sizeof(gdsb::Vertex32));
        input.read((char*)&graph_data.read_edges[current_edge_id].data, sizeof(gdsb::Weight));

        current_edge_id++;
        return true;
    };

    auto read_f_unweighted_undirected_32 = [&](std::ifstream& input)
    {
        gdsb::Vertex32 u;
        gdsb::Vertex32 v;

        input.read((char*)&u, sizeof(gdsb::Vertex32));
        input.read((char*)&v, sizeof(gdsb::Vertex32));

        graph_data.read_edges[current_edge_id].src = u;
        graph_data.read_edges[current_edge_id].dst = v;
        ++current_edge_id;

        graph_data.read_edges[current_edge_id].src = v;
        graph_data.read_edges[current_edge_id].dst = u;
        ++current_edge_id;

        return true;
    };

    auto read_f_weighted_undirected_32 = [&](std::ifstream& input)
    {
        gdsb::Vertex32 u;
        gdsb::Vertex32 v;
        EdgeData data;

        input.read((char*)&u, sizeof(gdsb::Vertex32));
        input.read((char*)&v, sizeof(gdsb::Vertex32));
        input.read((char*)&data, sizeof(gdsb::Weight));

        graph_data.read_edges[current_edge_id].src = u;
        graph_data.read_edges[current_edge_id].dst = v;
        graph_data.read_edges[current_edge_id].data = data;
        ++current_edge_id;

        graph_data.read_edges[current_edge_id].src = v;
        graph_data.read_edges[current_edge_id].dst = u;
        graph_data.read_edges[current_edge_id].data = data;
        ++current_edge_id;

        return true;
    };

    std::ifstream binary_graph(fname);
    gdsb::BinaryGraphHeader header = gdsb::read_binary_graph_header(binary_graph);
    uint64_t partition_edge_count = gdsb::partition_batch_count(header.edge_count, partition_id, partition_size);
    partition_edge_count = partition_edge_count + !header.directed * partition_edge_count;

    graph_data.read_edges = new Edge[partition_edge_count];

    auto const [vertex_count, edge_count] = [&]()
    {
        constexpr size_t size_of_edge_unweighted = sizeof(gdsb::Vertex32) + sizeof(gdsb::Vertex32);
        constexpr size_t size_of_edge_weighted = sizeof(gdsb::Vertex32) + sizeof(gdsb::Vertex32) + sizeof(gdsb::Weight);

        if constexpr (std::is_same_v<EdgeData, real_t>)
        {
            if (header.directed)
            {
                return gdsb::read_binary_graph_partition(binary_graph, header, std::move(read_f_weighted_32),
                                                         size_of_edge_weighted, partition_id, partition_size);
            }
            else
            {
                return gdsb::read_binary_graph_partition(binary_graph, header, std::move(read_f_weighted_undirected_32),
                                                         size_of_edge_weighted, partition_id, partition_size);
            }
        }

        if (header.directed)
        {
            return gdsb::read_binary_graph_partition(binary_graph, header, std::move(read_f_unweighted_32),
                                                     size_of_edge_unweighted, partition_id, partition_size);
        }
        else
        {

            return gdsb::read_binary_graph_partition(binary_graph, header, std::move(read_f_unweighted_undirected_32),
                                                     size_of_edge_unweighted, partition_id, partition_size);
        }
    }();

    assert(partition_edge_count == (edge_count + !header.directed * edge_count));

    graph_data.v_num_param = vertex_count;
    graph_data.read_e_num = edge_count + !header.directed * edge_count;
}

} // namespace kklib

template<typename T>
void write_graph(const char* fname, const T* es, const size_t e_num)
{
    FILE *out_f = fopen(fname, "w");
    assert(out_f != NULL);

    auto ret = fwrite(es, sizeof(T), e_num, out_f);
    assert(ret == e_num);
    fclose(out_f);
}

size_t next_endline_pos(FILE *f)
{
    size_t current_pos = ftell(f);
    while (true)
    {
        char ch;
        auto ret = fread(&ch, 1, 1, f);
        if (ret != 1 || ch == '\n')
        {
            break;
        }
        current_pos++;
    }
    return current_pos;
}

std::vector<size_t> partition_text_file(const char* fname, int partition_num)
{
    std::vector<size_t> partition_end;
    FILE *f = fopen(fname, "r");
    assert(f != NULL);
    fseek(f, 0, SEEK_END);
    size_t total_size = ftell(f);
    for (int p_i = 0; p_i < partition_num; p_i++)
    {
        size_t f_offset = total_size / partition_num * (p_i + 1);
        if (f_offset >= total_size)
        {
            partition_end.push_back(f_offset);
        } else
        {
            fseek(f, f_offset, SEEK_SET);
            partition_end.push_back(next_endline_pos(f));
        }
    }
    fclose(f);
    return partition_end;
}

bool read_edge_txt(FILE* f, kklib::EdgeNoData* edge) { return (2 == fscanf(f, "%u %u", &edge->src, &edge->dst)); }

bool read_edge_txt(FILE* f, Edge<real_t>* edge)
{
    return (3 == fscanf(f, "%u %u %f", &edge->src, &edge->dst, &edge->data));
}

template<typename T>
bool read_edge_txt(FILE* f, Edge<T>* edge)
{
    fprintf(stderr, "Edge type doesn't support reading from text\n");
    exit(1);
}

template<typename T>
void read_edgelist(const char* fname, int partition_id, int partition_num, Edge<T>* &edge, size_t &e_num)
{
    std::vector<size_t> partition_end = partition_text_file(fname, partition_num);
    size_t begin = (partition_id == 0 ? 0 : partition_end[partition_id - 1]);
    size_t end = partition_end[partition_id];
    FILE *f = fopen(fname, "r");
    assert(f != NULL);

    Edge<T> temp;
    e_num = 0;
    fseek(f, begin, SEEK_SET);
    while (ftell(f) < end)
    {
        if (read_edge_txt(f, &temp))
        {
            e_num++;
        }
    }
    edge = new Edge<T>[e_num];

    size_t e_i = 0;
    fseek(f, begin, SEEK_SET);
    while (ftell(f) < end)
    {
        if (read_edge_txt(f, &temp))
        {
            edge[e_i] = temp;
            e_i++;
        }
    }

    fclose(f);
}

void print_edge(FILE* f, const kklib::EdgeNoData* edge) { fprintf(f, "%u %u\n", edge->src, edge->dst); }

void print_edge(FILE* f, const Edge<real_t>* edge)
{
   fprintf(f, "%u %u %f\n", edge->src, edge->dst, edge->data);
}

template<typename T>
void print_edge(FILE* f, const Edge<T>* edge)
{
    fprintf(stderr, "Edge type doesn't support writing to text\n");
    exit(1);
}

template<typename T>
void write_edgelist(const char* fname, const Edge<T>* es, const size_t e_num)
{
    FILE *out_f = fopen(fname, "w");
    assert(out_f != NULL);
    for (size_t e_i = 0; e_i < e_num; e_i++)
    {
        print_edge(out_f, es + e_i);
    }
    fclose(out_f);
}

namespace kklib
{

template <typename EdgeData>
LoadedGraphData<EdgeData>
load_graph(VertexID v_num_param, char const* graph_path, bool load_as_undirected = false, GraphFormat graph_format = GF_Binary)
{
    //! In case the graph format is GF_Binary, the graph will be read using the
    //! kklib::read_graph() routine. If the graph is loaded as undirected, the
    //! edges in both directions will be added in place to
    //! graph_data.read_edges. This way, we remove unnecessary computation and
    //! allocation time loading the graph spend in the old if
    //! (load_as_undirected) branch. If graph_format is GF_Edgelist, we keep the
    //! old logic for compatibility and legacy reasons.
    load_as_undirected = graph_format == GF_Edgelist && load_as_undirected;

    LoadedGraphData<EdgeData> graph_data;
    partition_id_t const partition_num = get_mpi_size();
    partition_id_t const local_partition_id = get_mpi_rank();

    if (graph_format == GF_Binary)
    {
        kklib::read_graph<Edge<EdgeData>, EdgeData>(graph_path, local_partition_id, partition_num, graph_data);
    }
    else if (graph_format == GF_Edgelist)
    {
        read_edgelist(graph_path, local_partition_id, partition_num, graph_data.read_edges, graph_data.read_e_num);
    }
    else
    {
        if (kklib::is_master_process())
        {
            fprintf(stderr, "Unsupported graph formant");
        }
        exit(1);
    }

    if (load_as_undirected)
    {
        Edge<EdgeData>* undirected_edges = new Edge<EdgeData>[graph_data.read_e_num * 2];

#pragma omp parallel for
        for (edge_id_t e_i = 0; e_i < graph_data.read_e_num; e_i++)
        {
            undirected_edges[e_i * 2] = graph_data.read_edges[e_i];
            std::swap(graph_data.read_edges[e_i].src, graph_data.read_edges[e_i].dst);
            undirected_edges[e_i * 2 + 1] = graph_data.read_edges[e_i];
        }

        delete[] graph_data.read_edges;
        graph_data.read_edges = undirected_edges;
        graph_data.read_e_num *= 2;
    }

    return graph_data;
}

} // namespace kklib