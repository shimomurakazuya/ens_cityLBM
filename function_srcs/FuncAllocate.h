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
#include <cuda_runtime_api.h>
#include "cuda_safe_call.hpp"
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
        cudaMalloc(val, sizeof(T) * n);
        #endif
    }
    else if (memType == MemType::Pinned) {
        #ifdef CPU_CALCULATION__
        *val = (T *) malloc (sizeof(T) * n);
        #endif
        #ifdef GPU_CALCULATION__
        cudaMallocHost(val, sizeof(T) * n);
        #endif
    }
    else if (memType == MemType::Managed) {
        #ifdef CPU_CALCULATION__
        *val = (T *) malloc (sizeof(T) * n);
        #endif
        #ifdef GPU_CALCULATION__
        cudaMallocManaged(val, sizeof(T) * n);
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
        cudaFree(val);
        val = nullptr;
        #endif
    }
    else if (memType == MemType::Pinned) {
        #ifdef CPU_CALCULATION__
        free(val); val = nullptr;
        #endif
        #ifdef GPU_CALCULATION__
        cudaFree(val); val = nullptr;
        #endif
    }
    else if (memType == MemType::Managed) {
        #ifdef CPU_CALCULATION__
        free(val); val = nullptr;
        #endif
        #ifdef GPU_CALCULATION__
        cudaFree(val); val = nullptr;
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
        CUDA_SAFE_CALL(cudaMemcpy(dst, src, sizeof(T)*n, cudaMemcpyDefault));
        #else
        FuncMath::copy_array(dst, src, n);
        #endif
    }
}

};


#endif
