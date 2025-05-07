#pragma once

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <limits>
#include <random>
#include <sys/mman.h>
#include <sys/resource.h>
#include <type_traits> // std::is_integral

#include <kklib/constants.hpp>
#include <kklib/type.hpp>

namespace kklib
{

//! Base copied on May 17th, 2024 from
//! https://en.cppreference.com/w/cpp/named_req/Allocator.
template <typename T> class MemMappedAllocator
{
public:
    typedef T value_type;

    MemMappedAllocator() = default;

    template <typename U> constexpr MemMappedAllocator(const MemMappedAllocator<U>&) noexcept {}

    // [[nodiscard]] T* allocate(std::size_t n)
    T* allocate(std::size_t const n)
    {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
        {
            throw std::bad_array_new_length();
        }

        void* mmap_allocated_mem = mmap(nullptr, sizeof(T) * n, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

        if (mmap_allocated_mem == MAP_FAILED)
        {
            throw std::bad_alloc{};
        }

        return reinterpret_cast<T*>(mmap_allocated_mem);
    }

    void deallocate(T* allocated_mem, std::size_t const n) noexcept { munmap(allocated_mem, sizeof(T) * n); }
};

template <typename T> void next_dealloc_array(T* array, size_t const num)
{
    if (array != nullptr)
    {
        munmap(array, sizeof(T) * num);
    }
}

template <typename T> T* next_alloc_array(size_t const num)
{
    void* alloc_mem = mmap(nullptr, sizeof(T) * num, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);

    if (alloc_mem == MAP_FAILED)
    {
        throw std::bad_alloc{};
    }

    T* array = reinterpret_cast<T*>(alloc_mem);

    return array;
}

//! Simple random engine class using standard c++ pseudo random number generator
//! classes.
template <typename T> class RandomEngine
{
public:
    RandomEngine() = default;

    RandomEngine(std::mt19937::result_type const seed)
    : m_random_device()
    , m_engine(m_random_device())
    {
        set_seed(seed);
    }

    ~RandomEngine() = default;

    //! Use this operator to generate a pseudo random number within interval
    //! [min, max].
    //!
    //! Implementation detail: the distribution objects are fairly cheap to
    //! construct providing the [min, max] interval. Due to the constexpr if
    //! expression, there is no runtime spent on runtime evaluation of T.
    T operator()(T const min, T const max)
    {
        if constexpr (std::is_integral<T>())
        {
            std::uniform_int_distribution<T> distribution{ min, max };
            return distribution(m_engine);
        }
        else
        {
            std::uniform_real_distribution<T> distribution{ min, max };
            return distribution(m_engine);
        }
    }

    void set_seed(std::mt19937::result_type seed) { m_engine.seed(seed); }

private:
    std::random_device m_random_device;
    std::mt19937 m_engine{ m_random_device() };
};

} // namespace kklib

// Timer is used for performance profiling
class Timer
{
public:
    void restart();

    double duration();

private:
    std::chrono::time_point<std::chrono::system_clock> m_start = std::chrono::system_clock::now();
};

class MessageBuffer
{
    char padding[kklib::l1_cache_line_size - sizeof(size_t) - sizeof(bool) - sizeof(VertexID) - sizeof(void*)];
    bool is_private_data;

public:
    size_t sz;
    size_t count;
    void* data;
    MessageBuffer()
    {
        sz = 0;
        count = 0;
        data = nullptr;
        is_private_data = false;
    }
    MessageBuffer(size_t _sz, void* _data = nullptr)
    : MessageBuffer()
    {
#ifdef UNIT_TEST
        assert(data == nullptr);
#endif
        alloc(_sz, _data);
    }
    void alloc(size_t _sz, void* _data = nullptr)
    {
        if (data != nullptr && is_private_data == true)
        {
            delete[](char*) data;
        }
        sz = _sz;
        count = 0;
        if (_data == nullptr)
        {
            data = new char[sz];
            is_private_data = true;
        }
        else
        {
            data = _data;
            is_private_data = false;
        }
    }

    ~MessageBuffer()
    {
        if (data != nullptr && is_private_data == true)
        {
            delete[](char*) data;
        }
    }

    template <typename data_t> void write(data_t* val)
    {
        ((data_t*)data)[count++] = *val;
#ifdef UNIT_TEST
        self_check<data_t>();
#endif
    }

    void clear() { count = 0; }

    template <typename data_t> void self_check() { assert(sz >= count * sizeof(data_t)); }
};