#pragma once
#ifndef PACKUNPACKGPU_H_
#define PACKUNPACKGPU_H_


#include <cmath>
#include "defineCal.h"
#include "definePrecision.h"
#include "FuncFP.h"


namespace PackUnpackDef {
//    constexpr int TheadMAX = 1024;
    constexpr int TheadMAX = 512;
};


namespace  PackUnpackCPU {


template <typename T0, typename T1>
void pack_cpu(
    const int* list,
    const T0*  val,
          T1*  buff,
    const int  num_list
    ) noexcept
{
#pragma omp parallel for
    for (int i=0; i<num_list; i++) { buff[i] = val[list[i]]; }
}



#ifdef USE_NVCC
template <typename T>
void pack_half_cpu(
    const int*  list,
    const T*    val,
          fp16* buff,
    const int   num_list
    ) noexcept
{
#pragma omp parallel for
    for (int i=0; i<num_list; i++) { buff[i] = FP::r<fp16>( val[list[i]] ); }
}
#endif


template <typename T0, typename T1>
void unpack_cpu(
    const int* list,
          T0*  val,
    const T1*  buff,
    const int  num_list
    ) noexcept
{
#pragma omp parallel for
    for (int i=0; i<num_list; i++) { val[list[i]] = buff[i]; }
}


#ifdef USE_NVCC
template <typename T>
void unpack_half_cpu(
    const int*  list,
          T*    val,
    const fp16* buff,
    const int   num_list
    ) noexcept
{
#pragma omp parallel for
    for (int i=0; i<num_list; i++) { val[list[i]] = FP::r<T>( buff[i] ); }
}
#endif


}; // PackUnpackCPU //


namespace  PackUnpackGPU {
#ifdef USE_NVCC

template <typename T0, typename T1>
__global__
void pack_gpu(
    const int* list,
    const T0*  val,
          T1*  buff,
    const int  num_list
    )
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i >= num_list) { return; }

    buff[i] = FP::r<T1>( val[list[i]] );
}


template <typename T>
__global__
void pack_half_gpu(
    const int* list,
    const T*    val,
          fp16* buff,
    const int   num_list
    )
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i >= num_list) { return; }

    buff[i] = FP::r<fp16>( val[list[i]] );
}


template <typename T0, typename T1>
__global__
void unpack_gpu(
    const int* list,
          T0*  val,
    const T1*  buff,
    const int  num_list
    )
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i >= num_list) { return; }

    val[list[i]] = buff[i];
}


template <typename T>
__global__
void unpack_half_gpu(
    const int* list,
          T*    val,
    const fp16* buff,
    const int   num_list
    )
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i >= num_list) { return; }

    val[list[i]] = FP::r<T>( buff[i] );
}

#endif // ifdef USE_NVCC
}; // PackUnpackGPU //


namespace  PackUnpack {

// pack //
template <typename T0, typename T1>
void pack_host(
    const int* list,
    const T0*  val,
          T1*  buff,
    const int  num_list
    ) noexcept
{
    PackUnpackCPU::pack_cpu(list, val, buff, num_list);
}


#ifdef USE_NVCC
template <typename T>
void pack_half_host(
    const int*  list,
    const T*    val,
          fp16* buff,
    const int   num_list
    ) noexcept
{
    PackUnpackCPU::pack_half_cpu(list, val, buff, num_list);
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T0, typename T1>
void pack_device(
    const int* list,
    const T0*  val,
          T1*  buff,
    const int  num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::pack_gpu<T0, T1>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::pack_gpu <<< grid, block >>> (list, val, buff, num_list);
    cudaDeviceSynchronize();
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T0, typename T1>
void pack_device(
    const cudaStream_t& streamX,
    const int* list,
    const T0*  val,
          T1*  buff,
    const int  num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::pack_gpu<T0, T1>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::pack_gpu <<< grid, block, 0, streamX >>> (list, val, buff, num_list);
    cudaStreamSynchronize(streamX);
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T>
void pack_half_device(
    const int* list,
    const T*    val,
          fp16* buff,
    const int   num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::pack_half_gpu<T>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::pack_half_gpu <<< grid, block >>> (list, val, buff, num_list);
    cudaDeviceSynchronize();
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T>
void pack_half_device(
    const cudaStream_t& streamX,
    const int* list,
    const T*    val,
          fp16* buff,
    const int   num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::pack_half_gpu<T>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::pack_half_gpu <<< grid, block, 0, streamX >>> (list, val, buff, num_list);
    cudaStreamSynchronize(streamX);
}
#endif


template <int NX_LEAF, typename T0, typename T1>
void pack(
    const int* list,
    const T0*  val,
          T1*  buff,
    const int  num_list
    ) noexcept
{
#if defined(USE_NVCC) && defined(GPU_CALCULATION__)
    pack_device<NX_LEAF>(list, val, buff, num_list);
#else
    pack_host(list, val, buff, num_list);
#endif
}


// unpack //
template <typename T0, typename T1>
void unpack_host(
    const int* list,
          T0*  val,
    const T1*  buff,
    const int  num_list
    ) noexcept
{
    PackUnpackCPU::unpack_cpu(list, val, buff, num_list);
}


#ifdef USE_NVCC
template <typename T>
void unpack_half_host(
    const int*  list,
          T*    val,
    const fp16* buff,
    const int   num_list
    ) noexcept
{
    PackUnpackCPU::unpack_half_cpu(list, val, buff, num_list);
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T0, typename T1>
void unpack_device(
    const int* list,
          T0*  val,
    const T1*  buff,
    const int  num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::unpack_gpu<T0, T1>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::unpack_gpu <<< grid, block >>> (list, val, buff, num_list);
    cudaDeviceSynchronize();
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T0, typename T1>
void unpack_device(
    const cudaStream_t& streamX,
    const int* list,
          T0*  val,
    const T1*  buff,
    const int  num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::unpack_gpu<T0, T1>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::unpack_gpu <<< grid, block, 0, streamX >>> (list, val, buff, num_list);
    cudaStreamSynchronize(streamX);
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T>
void unpack_half_device(
    const int* list,
          T*    val,
    const fp16* buff,
    const int   num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::unpack_half_gpu<T>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::unpack_half_gpu <<< grid, block >>> (list, val, buff, num_list);
    cudaDeviceSynchronize();
}
#endif


#ifdef USE_NVCC
template <int NX_LEAF, typename T>
void unpack_half_device(
    const cudaStream_t& streamX,
    const int* list,
          T*    val,
    const fp16* buff,
    const int   num_list
    ) noexcept
{
    const int NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    constexpr int num_threads = ( NN_LEAF <= PackUnpackDef::TheadMAX ) ? NN_LEAF : PackUnpackDef::TheadMAX;

    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
    const dim3 grid( (num_list+num_threads-1)/num_threads, 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::unpack_half_gpu<T>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::unpack_half_gpu <<< grid, block, 0, streamX >>> (list, val, buff, num_list);
    cudaStreamSynchronize(streamX);
}
#endif


template <int NX_LEAF, typename T0, typename T1>
void unpack(
    const int* list,
          T0*  val,
    const T1*  buff,
    const int  num_list
    ) noexcept
{
#if defined(USE_NVCC) && defined(GPU_CALCULATION__)
    unpack_device<NX_LEAF>(list, val, buff, num_list);
#else
    unpack_host(list, val, buff, num_list);
#endif
}


} // PackUnpack //


#endif
