#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <cstdlib>

#include "TestOklahoma.h"


int main(int argc, char* argv[])
{
#ifdef USE_NVCC
    cudaSetDevice(0);
#endif

    // MPI_Init //
#if 0
//    int  required = MPI_THREAD_FUNNELED;
//    int  required = MPI_THREAD_SERIALIZED;
    int  required = MPI_THREAD_MULTIPLE;
    int  provided;
    MPI_Init_thread(&argc, &argv, required, &provided);
    if (required != provided) { std::cout << "error : MPI_Init_thread" << std::endl; exit(0); }
#else
    MPI_Init(&argc, &argv);
#endif

    int  ncpu, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);

    for (int i=0; i<ncpu; i++) {
        if (i == rank) { std::cout << "rank / ncpu = " << rank << " / " << ncpu << std::endl; }
        MPI_Barrier(MPI_COMM_WORLD);
    }

#ifdef USE_NVCC
    Parser  parser(argc, argv);
    cudaSetDevice(rank%parser.optionParser()->gpu_per_node());
#endif


#ifdef USE_SHARED_MEMORY__
    if (rank == 0) { std::cout << "use_shared_memory" << std::endl; }
#endif


    // validation test //
    TestOklahoma  testOklahoma;
    testOklahoma.Flow(argc, argv);


    // finalize //
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) { std::cout << "finish" << std::endl; }
    MPI_Finalize();

    return  EXIT_SUCCESS;
}
