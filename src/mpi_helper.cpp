#include <kklib/mpi_helper.hpp>

#include <cassert>

kklib::MPI_Instance::MPI_Instance(int* argc, char*** argv)
{
    constexpr int required_thread_support = MPI_THREAD_MULTIPLE;
    int provided_thread_support = required_thread_support;
    MPI_Init_thread(argc, argv, required_thread_support, &provided_thread_support);
    assert(required_thread_support == provided_thread_support);
}


bool kklib::is_master_process() { return get_mpi_rank() == 0; }

int get_mpi_rank()
{
    int a;
    MPI_Comm_rank(MPI_COMM_WORLD, &a);
    return a;
}

int get_mpi_size()
{
    int a;
    MPI_Comm_size(MPI_COMM_WORLD, &a);
    return a;
}
