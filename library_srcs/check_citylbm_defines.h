#ifndef CHECK_CITYLBM_DEFINES_H_
#define CHECK_CITYLBM_DEFINES_H_

#ifdef CHECK_DEFINES

// invalid options
#if defined(USE_CONST_BOUNDARY_WOWRF) && (!defined(CONST_BOUNDARY_WO_WRF_U) || !defined(CONST_BOUNDARY_WO_WRF_V))
#error USE_CONST_BOUNDARY_WOWRF: error: CONST_BOUNDARY_WO_WRF_U and/or CONST_BOUNDARY_WO_WRF_V are not defined
#endif

#if defined(USE_PBVR) && defined(CPU_CALCULATION__)
#error USE_PBVR: error: CPU compilation is not supported
#endif

#if defined(USE_PBVR) && !defined(USE_VALUE_STAT)
#error USE_PBVR: error: pbvr needs USE_VALUE_STAT
#endif

#if defined(USE_VALUE_STAT) && defined(MESHVALUES_IO_MASK_HALO)
#error USE_VALUE_STAT: cannot be used with `MESHVALUES_IO_MASK_HALO`
#endif

// valid but not recommended
#ifdef VVTARGET_TestOklahoma
#ifndef EXPECT_DX_MESH
#warning `EXPECT_DX_MESH` is not defined. Using dxf=4m by default
#endif
#endif

#ifdef REGRESSION_TEST
#warning `REGRESSION_TEST` is tentative option for debug
#endif

#ifdef USE_WRF0
#warning `USE_WRF0` is tentative option for PBVR
#endif

#if defined(USE_TEMPORAL_BLOCKING) && !defined(TB_COUNT)
#warning `USE_TEMPORAL_BLOCKING` needs `TB_COUNT=N` (N = 1 or 3)
#endif

// print info
#ifdef CHECK_DEFINES_INFO

#include <iostream>
#include <type_traits>

namespace check_defines {

template<class ostream>
void inspect(ostream&& stream) {

    stream << "defined macro:" << std::endl;

    #ifdef USE_NVCC
    stream << "  USE_NVCC" << std::endl;
    #endif

    #ifdef GPU_CALCULATION__
    stream << "  GPU_CALCULATION__" << std::endl;
    #ifdef USE_SHARED_MEMORY__
    std::cout << "  USE_SHARED_MEMORY__" << std::endl;
    #endif
    #endif

    #ifdef CPU_CALCULATION__
    stream << "  CPU_CALCULATION__" << std::endl;
    #endif

    #ifdef SINGLE_PRECISION_CALCULATION
    stream << "  SINGLE_PRECISION_CALCULATION" << std::endl;
    #endif
    #ifdef DOUBLE_PRECISION_CALCULATION
    stream << "  DOUBLE_PRECISION_CALCULATION" << std::endl;
    #endif

    #ifdef FLOAT16_CAL
    stream << "  FLOAT16_CAL" << std::endl;
    #endif
    #ifdef BFLOAT16_CAL
    stream << "  BFLOAT16_CAL" << std::endl;
    #endif

    #ifdef USE_TEMPORAL_BLOCKING
    stream << "  USE_TEMPORAL_BLOCKING" << std::endl;
    #endif
    #ifdef TB_COUNT
    stream << "  TB_COUNT=" << TB_COUNT << std::endl;
    #endif

    #ifdef USE_PBVR
    stream << "  USE_PBVR" << std::endl;
    #endif

    #ifdef CHECK_MPI_TIME
    stream << "  CHECK_MPI_TIME" << std::endl;
    #endif
    
    #ifdef REGRESSION_TEST
    stream << "  REGRESSION_TEST: this is tentative option for regression test. (not recommended for release build)" << std::endl;
    #endif

    #ifdef USE_WRF0
    stream << "  USE_WRF0: this is tentative option for PBVR" << std::endl;
    #endif
    
    #ifdef NO_FIELD_WRITEDAT
    stream << "  NO_FIELD_WRITEDAT: disable meshValues_io_ etc. (deprecated)" << std::endl;
    #endif

    #ifdef EXPECT_DX_MESH
    stream << "  EXPECT_DX_MESH=" << EXPECT_DX_MESH << ": check dx_f. assersion fails if dx_f != EXPECT_DX_MESH" << std::endl;
    #endif

    #ifdef NO_IOTHREAD
    stream << "  NO_IOTHREAD: disable thread fork for I/O" << std::endl;
    #endif
    
    #ifdef NO_POSTPROCESS_MONITOR
    stream << "  NO_POSTPROCESS_MONITOR" << std::endl;
    #endif

    #ifdef FASTPOST_IGNORE_MISSING
    stream << "  FASTPOST_IGNORE_MISSING: FastPostprocess won't throw runtime_error when some monitors are missing" << std::endl;
    #endif

    #ifdef USE_VALUE_STAT
    stream << "  USE_VALUE_STAT" << std::endl;
    #endif
    
    #ifdef NO_IODATA_PARAVIEW
    stream << "  NO_IODATA_PARAVIEW" << std::endl;
    #else
    #ifdef PARAVIEW_FULL
    stream << "  PARAVIEW_FULL: output full-size vtk data" << std::endl;
    #endif
    #ifdef PARAVIEW_DOWNSIZE_2
    stream << "  PARAVIEW_DOWNSIZE_2: output half-pruned vtk data" << std::endl;
    #endif
    #ifdef PARAVIEW_DOWNSIZE_4
    stream << "  PARAVIEW_DOWNSIZE_4: output 3/4-pruned vtk data" << std::endl;
    #endif
    #ifdef PARAVIEW_SLICE
    stream << "  PARAVIEW_SLICE: output vtk data on horizontal and vertical slices" << std::endl;
    #endif
    #ifdef PARAVIEW_CELLDATA
    stream << "  PARAVIEW_CELLDATA: switch PointData to CellData. (for chaging visibility)" << std::endl;
    #endif
    #ifdef PARAVIEW_ENS0
    stream << "  PARAVIEW_ENS0: output vtk data on ens0, but skip them on ens1,2,..." << std::endl;
    #endif

    #endif // NO_IODATA_PARAVIEW

    #ifdef NO_IODATA_RESTART
    stream << "  NO_IODATA_RESTART" << std::endl;
    #endif
    
    #ifdef NO_CITYBOX
    stream << "  NO_CITYBOX: disable citybox refinement for Oklahoma. (recommended)" << std::endl;
    #endif
    
    #ifdef ENSEMBLE_SPIKE
    stream << "  ENSEMBLE_SPIKE: perturbate arrangement of roughness blocks. (recommended)" << std::endl;
    #endif

    #ifdef ENSEMBLE_NUDGING
    stream << "  ENSEMBLE_NUDGING: perturbate nudging parameter" << std::endl;
    #endif

    #ifdef ENSEMBLE_UVWT
    stream << "  ENSEMBLE_UVWT: perturbate input variable from WRF in ensemble simulations" << std::endl;
    #endif
    
    #ifdef SCALAR_NON_PRIODIC
    stream << "  SCALAR_NON_PRIODIC: zeroset scala_n on NWSE bound. (recommended)" << std::endl;
    #endif
    
    #ifdef SC_SCALAR_ALWAYS
    stream << "  SC_SCALAR_ALWAYS: releasing scalar always. (for PBVR, ML, or debug)" << std::endl;
    #endif
    
    #ifdef MESHVALUES_IO_MASK_HALO
    stream << "  MESHVALUES_IO_MASK_HALO: set NA on halo grids in meshValues_io_. (for debug)" << std::endl;
    #endif

    #ifdef COUT_MONITOR_PRECISION
    stream << "  COUT_MONITOR_PRECISION=" << COUT_MONITOR_PRECISION << std::endl;
    #endif

    stream << std::endl;

} // void inspect
} // namespace

#endif // ifdef CHECK_DEFINES_INFO

#endif // ifdef CHECK_DEFINES




#endif // include guard
