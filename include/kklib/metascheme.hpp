#pragma once

#include <kklib/type.hpp>

#include <assert.h>
#include <stdio.h>
#include <vector>

typedef int scheme_id_t;
typedef int meta_state_t;

struct MetapathState
{
    scheme_id_t scheme_id;
    meta_state_t state;
};

struct WeightedMetaData
{
    real_t weight;
    meta_state_t meta_info;
    meta_state_t get_meta() { return meta_info; }
    real_t get_weight() { return weight; }
};

struct UnweightedMetaData
{
    meta_state_t meta_info;
    meta_state_t get_meta() { return meta_info; }
    real_t get_weight() { return 1.0; }
};

std::vector<std::vector<std::vector<bool>>> read_metapath_schemes(const char* path);

void write_metapath_schemes(std::vector<std::vector<std::vector<bool>>> schemes, const char* path);

typedef uint8_t scheme_mask_t;

std::vector<std::vector<scheme_mask_t>> get_scheme_mask(std::vector<std::vector<std::vector<bool>>>& schemes);