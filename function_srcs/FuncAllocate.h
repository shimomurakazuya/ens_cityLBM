#pragma once
#ifndef FUNCALLOCATE_H_
#define FUNCALLOCATE_H_


#include <iostream>
#include <cstdlib>
#include <numeric>
#include "defineCal.h"
#include "defineCUDA.h"
#include "defineMemory.h"
#include "bytes_type.hpp"
#include "FuncMath.h"

#ifdef GPU_CALCULATION__
  #if defined(USE_NVCC)
    #include <cuda_runtime_api.h>
    #include "cuda_safe_call.hpp"
  #elif defined(ENABLE_HIP)
    #include "hip_safe_call.hpp"
  #endif
#endif

namespace  FuncAllocate {

template <typename  T>
void  fill_value(
          T*        val,
    const T         ini,
    const int       n,
    const MemType   memType
    )
{
    if      (memType == MemType::Host || memType == MemType::Managed) {
        std::fill(val, val+n, ini);
    }
    else if (memType == MemType::Device) {
        #ifdef GPU_CALCULATION__
        #endif
    }
    else {
    }
}


template <typename  T>
void  allocate_value(
          T**       val,
    const int       n,
    const MemType   memType
    )
{
    if      (memType == MemType::Host) {
        *val = (T *) malloc (sizeof(T) * n);
    }
    else if (memType == MemType::Device) {
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            cudaMalloc(val, sizeof(T) * n);
          #elif defined(ENABLE_HIP)
            hipMalloc(val, sizeof(T) * n);
          #endif
        #endif
    }
    else if (memType == MemType::Pinned) {
        #ifdef CPU_CALCULATION__
        *val = (T *) malloc (sizeof(T) * n);
        #endif
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            cudaMallocHost(val, sizeof(T) * n);
          #elif defined(ENABLE_HIP)
            hipHostMalloc(val, sizeof(T) * n);
          #endif
        #endif
    }
    else if (memType == MemType::Managed) {
        #ifdef CPU_CALCULATION__
        *val = (T *) malloc (sizeof(T) * n);
        #endif
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            cudaMallocManaged(val, sizeof(T) * n);
          #elif defined(ENABLE_HIP)
            // [RK] Experimental 
            // hip device memory serves as managed memory on AMD GPU
            // so we use device memory in place of managed memory
            // managed memory on AMD seems fairly slow
            // hipMallocManaged(val, sizeof(T) * n);
            hipMalloc(val, sizeof(T) * n);
          #endif
        #endif
    }
    else {
        *val = nullptr;
    }
    #ifdef REGRESSION_TEST
    /// set uninitialized value = 0xffff...
    util::bytes_union<T> tmp; tmp.b = -1;
    fill_value( *val, tmp.value, n, memType);
    #endif
}


template <typename  T>
void  release_value(
          T*&       val,
    const MemType   memType
    )
{
    if (val == nullptr)  { return; }

    if      (memType == MemType::Host) {
        free(val);
        val = nullptr;
    }
    else if (memType == MemType::Device) {
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            cudaFree(val);
          #elif defined(ENABLE_HIP)
            hipFree(val);
          #endif
        val = nullptr;
        #endif
    }
    else if (memType == MemType::Pinned) {
        #ifdef CPU_CALCULATION__
        free(val); val = nullptr;
        #endif
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            cudaFree(val); val = nullptr;
          #elif defined(ENABLE_HIP)
            hipFree(val); val = nullptr;
          #endif
        #endif
    }
    else if (memType == MemType::Managed) {
        #ifdef CPU_CALCULATION__
        free(val); val = nullptr;
        #endif
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            cudaFree(val); val = nullptr;
          #elif defined(ENABLE_HIP)
            hipFree(val); val = nullptr;
          #endif
        #endif
    }
    else {
    }
}

template<typename T>
void copy_values(
          T*        dst,
    const T*        src,
    const int       n,
    const MemType   memType
)
{
    if(memType == MemType::Host || memType == MemType::Pinned) {
        FuncMath::copy_array(dst, src, n);
    } else if(memType == MemType::Device || memType == MemType::Managed) {
        #ifdef GPU_CALCULATION__
          #if defined(USE_NVCC)
            CUDA_SAFE_CALL(cudaMemcpy(dst, src, sizeof(T)*n, cudaMemcpyDefault));
          #elif defined(ENABLE_HIP)
            HIP_SAFE_CALL(hipMemcpy(dst, src, sizeof(T)*n, hipMemcpyDefault));
          #endif
        #else
        FuncMath::copy_array(dst, src, n);
        #endif
    }
}

};


#endif
