#pragma once

#include <cstdint>

namespace kklib
{
constexpr uint32_t l1_cache_line_size = 64;

constexpr uint32_t thread_local_buf_capacity =
#ifdef UNIT_TEST
    16;
#else
    1024;
#endif

constexpr uint32_t omp_parallel_threshold =
#ifdef UNIT_TEST
    10;
#else
    4000;
#endif

constexpr uint32_t parallel_chunk_size =
#ifdef UNIT_TEST
    4;
#else
    128;
#endif

constexpr uint32_t phased_execution_threshold =
#ifdef UNIT_TEST
    100;
#else
    500000;
#endif

constexpr uint32_t distributed_execution_ctx_phasenum =
#ifdef UNIT_TEST
    5;
#else
    16;
#endif

constexpr uint32_t foot_print_chunk_size =
#ifdef UNIT_TEST
    16;
#else
    65536;
#endif

} // namespace kklib
