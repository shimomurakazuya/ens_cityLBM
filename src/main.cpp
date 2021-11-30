
#include "defineAMR.h"
#include "defineCal.h"
#include "defineCUDA.h"
#include "defineFilenames.h"
#include "defineFluidProperty.h"
#include "defineLBM.h"
#include "defineMemory.h"
#include "defineObjBC.h"
#include "definePrecision.h"
#include "defineSGS.h"

#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <ios>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>

#include "mpi_wrapper.hpp"
#include "MPICommEnsemble.h"
#include "runtime_error.hpp"

#include "VVTarget.h"

#include "check_citylbm_defines.h"


int main(int argc, char* argv[])
{
    auto&& opt = OptionParser("main", "[options...]");
    runtime_assert(1 == opt.parse_args(argc, argv), "InvalidConfiguration: OptionParser::parse_args failed");
    runtime_assert(opt.check_narguments(0), "InvalidConfiguration: OptionParser::check_narguments failed");

    // mpi init //
    #ifdef GPU_CALCULATION__ 
    cudaSetDevice(0);
    #endif
    MPI_Init(&argc, &argv);
    auto&& comm = MPICommEnsemble(opt.n_ensemble_members(), opt.ofs_ensemble_idx());
    comm.cout_info();
    if(comm.is_rank0()) { 
        std::cout << "target: " << typeid(VVTarget).name() << std::endl; 
        check_defines::inspect(std::cout); 
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // distribute multi-gpus by rank_citylbm //
    #ifdef GPU_CALCULATION__ 
    cudaSetDevice(comm.world().rank() % opt.gpu_per_node());
    #endif


    // validation test //
    try {
        ::VVTarget(opt, comm).Flow();
    } catch(const std::runtime_error& e) {
        if(comm.is_rank0()) {
            std::cerr << "!! program failed due to runtime_error: " << std::endl << e.what() << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    } catch(const std::exception& e) {
        if(comm.is_rank0()) {
            std::cerr << "!! program failed due to std::exception: " << std::endl << e.what() << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 2);
    } catch(...) {
        if(comm.is_rank0()) {
            std::cerr << "!! program failed due to unkown error" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 3);
    }

    // finalize //
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "finish" << std::endl;
    MPI_Finalize();

    return  EXIT_SUCCESS;
} 
