#pragma once

#include <kklib/constants.hpp>

#include <algorithm>
#include <cstdint>
#include <mutex>

using VertexID = uint32_t;
using edge_id_t = uint64_t;
using partition_id_t = uint8_t;
using task_counter_t = uint32_t;
using walker_id_t = uint32_t;
using step_t = uint32_t;
using real_t = float;
using dist_counter_t = uint32_t;

enum MPIMessageTag
{
    Tag_ShuffleGraph,
    Tag_Msg
};

template <typename T> class Message
{
public:
    VertexID dst_vertex_id;
    T data;
};

struct DistributedExecutionCtx
{
    std::mutex phase_locks[kklib::distributed_execution_ctx_phasenum];
    int unlocked_phase;
    size_t** progress;

public:
    DistributedExecutionCtx() { progress = nullptr; }
};