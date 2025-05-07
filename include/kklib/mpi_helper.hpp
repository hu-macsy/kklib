#pragma once

#include <stdexcept>
#include <type_traits>

#include <mpi.h>

namespace kklib
{

//! Simple wrapper to init and finalize MPI environment.
class MPI_Instance
{
public:
    MPI_Instance(int* argc, char*** argv);

    ~MPI_Instance() { MPI_Finalize(); }
};


//! Returns the corresponding MPI datatype given T. Inspired by gemini core
//! mpi.hpp function get_mpi_data_type(). The difference of this function is
//! that the result will be available at compile time due to the use of
//! constexpr. Only if the type is not supported, this function will throw a
//! std::invalid_argument exception (at runtime).
template <typename T> constexpr MPI_Datatype deduce_mpi_data_type()
{
    if constexpr (std::is_same<T, char>::value)
    {
        return MPI_CHAR;
    }
    else if constexpr (std::is_same<T, signed char>::value)
    {
        return MPI_SIGNED_CHAR;
    }
    else if constexpr (std::is_same<T, unsigned char>::value)
    {
        return MPI_UNSIGNED_CHAR;
    }
    // We are not handling byte or wchar!
    else if constexpr (std::is_same<T, short>::value)
    {
        return MPI_SHORT;
    }
    else if constexpr (std::is_same<T, unsigned short>::value)
    {
        return MPI_UNSIGNED_SHORT;
    }
    else if constexpr (std::is_same<T, int>::value)
    {
        return MPI_INT;
    }
    else if constexpr (std::is_same<T, unsigned int>::value)
    {
        return MPI_UNSIGNED;
    }
    else if constexpr (std::is_same<T, long>::value)
    {
        return MPI_LONG;
    }
    else if constexpr (std::is_same<T, unsigned long>::value)
    {
        return MPI_UNSIGNED_LONG;
    }
    else if constexpr (std::is_same<T, float>::value)
    {
        return MPI_FLOAT;
    }
    else if constexpr (std::is_same<T, double>::value)
    {
        return MPI_DOUBLE;
    }
    else if constexpr (std::is_same<T, long double>::value)
    {
        return MPI_LONG_DOUBLE;
    }
    else if constexpr (std::is_same<T, long long int>::value)
    {
        return MPI_LONG_LONG_INT;
    }
    else if constexpr (std::is_same<T, unsigned long long>::value)
    {
        return MPI_UNSIGNED_LONG_LONG;
    }
    else
    {
        throw std::invalid_argument("Type not supported by MPI.");
    }
}

bool is_master_process();

} // namespace kklib

int get_mpi_rank();

int get_mpi_size();
