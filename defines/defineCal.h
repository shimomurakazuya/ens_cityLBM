#pragma once
#ifndef DEFINECAL_H_
#define DEFINECAL_H_


#ifdef USE_NVCC
  #define GPU_CALCULATION__
  #define USE_SHARED_MEMORY__
#else
  #define CPU_CALCULATION__
#endif


#define STREAM_COLLISION_OBJ_OPT
//#define OPT_MPI_COMM


#ifdef USE_NVCC
//    constexpr int  opt_comm_send_recv_size      = 131072; // KB //
//    constexpr int  opt_comm_send_recv_size      = 262144; // KB //
//    constexpr int  opt_comm_send_recv_size      = 393216; // KB //
//    constexpr int  opt_comm_send_recv_size      = 524288; // KB //
//    constexpr int  opt_comm_send_recv_size      = 1048576; // KB //
//    constexpr int  opt_comm_send_recv_size      = 2097152; // KB //
//    constexpr int  opt_comm_send_recv_size      = 4194304; // KB //
    constexpr int  opt_comm_send_recv_size      = 8388608; // KB // opt
//    constexpr int  opt_comm_send_recv_size      = 12582912; // KB // opt
//    constexpr int  opt_comm_send_recv_size      = 16777216; // KB //
//    constexpr int  opt_comm_send_recv_size      = 33554432; // KB // slow //

//    constexpr int  opt_comm_send_recv_size      = 1073741824; // KB //
#else
    // CPU //
    constexpr int  opt_comm_send_recv_size      = 1073741824; // KB //
#endif


// cuda //
#ifdef USE_NVCC
  #define __GLOBAL__  __global__
  #define __HOST__    __host__
  #define __DEVICE__  __device__
#else
  #define __GLOBAL__
  #define __HOST__
  #define __DEVICE__
#endif

#define __HD__  __HOST__ __DEVICE__ 

//const int CUDA_THREAD_MAX = 32;
//const int CUDA_THREAD_MAX = 64;
//const int CUDA_THREAD_MAX = 128;
const int CUDA_THREAD_MAX = 256;
//const int CUDA_THREAD_MAX = 512;


namespace OuterBoundaryConditions {
    // by defined macros DEFINECAL_PERIODIC_{X,Y,Z}
    // or, set Oklahoma config as defualt

    #ifdef DEFINECAL_PERIODIC_X
    constexpr bool  OuterPeriodicCommX = DEFINECAL_PERIODIC_X;
    #else
    constexpr bool  OuterPeriodicCommX = true;
    #endif

    #ifdef DEFINECAL_PERIODIC_Y
    constexpr bool  OuterPeriodicCommY = DEFINECAL_PERIODIC_Y;
    #else
    constexpr bool  OuterPeriodicCommY = true;
    #endif

    #ifdef DEFINECAL_PERIODIC_Z
    constexpr bool  OuterPeriodicCommZ = DEFINECAL_PERIODIC_Z;
    #else
    constexpr bool  OuterPeriodicCommZ = false;
    #endif
};


#endif
