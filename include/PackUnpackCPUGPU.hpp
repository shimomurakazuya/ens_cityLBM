#include "PackUnpackCPUGPU.h"
#include "PackUnpackKernel.h"
#include <cmath>
#include "FuncLoop.h"


namespace  PackUnpack {


// pack //
template <typename T>
void pack_host(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    )
{
    PackUnpackCPU::pack_cpu(list, val, buff, num_list);
}


template <typename T>
void pack_device(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    )
{
#ifdef GPU_CALCULATION__ // GPU //
    constexpr int num_threads = ( DefAMR::NN_LEAF <= 1024 ) ? DefAMR::NN_LEAF : 1024;

    const dim3 block(num_threads, 1, 1);
    const dim3 grid( ceil(num_list/num_threads), 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::pack_gpu<T>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::pack_gpu <<< grid, block >>> (list, val, buff, num_list);
    cudaDeviceSynchronize();
#endif
}


template <typename T>
void pack(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    )
{
#ifdef CPU_CALCULATION__ // CPU : OMP //
//    PackUnpackCPU::pack_cpu(list, val, buff, num_list);
    pack_host(list, val, buff, num_list);
#endif

#ifdef GPU_CALCULATION__ // GPU //
//    constexpr int num_threads = ( DefAMR::NN_LEAF <= 1024 ) ? DefAMR::NN_LEAF : 1024;
//
//    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
//
//    cudaFuncSetCacheConfig( PackUnpackGPU::pack_gpu<T>,    cudaFuncCachePreferL1 );
//
//    PackUnpackGPU::pack_gpu <<< grid, block >>> (list, val, buff, num_list);
//    cudaDeviceSynchronize();
    pack_device(list, val, buff, num_list);
#endif
}


// unpack //
template <typename T>
void unpack_host(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    )
{
    PackUnpackCPU::unpack_cpu(list, val, buff, num_list);
}


template <typename T>
void unpack_device(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    )
{
#ifdef GPU_CALCULATION__ // GPU //
    constexpr int num_threads = ( DefAMR::NN_LEAF <= 1024 ) ? DefAMR::NN_LEAF : 1024;

    const dim3 block(num_threads, 1, 1);
    const dim3 grid( ceil(num_list/num_threads), 1, 1 );

    cudaFuncSetCacheConfig( PackUnpackGPU::unpack_gpu<T>,    cudaFuncCachePreferL1 );

    PackUnpackGPU::unpack_gpu<T> <<< grid, block >>> (list, val, buff, num_list);
    cudaDeviceSynchronize();
#endif
}


template <typename T>
void unpack(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    )
{
#ifdef CPU_CALCULATION__ // CPU : OMP //
//    PackUnpackCPU::unpack_cpu(list, val, buff, num_list);
    unpack_host(list, val, buff, num_list);
#endif

#ifdef GPU_CALCULATION__ // GPU //
//    constexpr int num_threads = ( DefAMR::NN_LEAF <= 1024 ) ? DefAMR::NN_LEAF : 1024;
//
//    const dim3 block(num_threads, 1, 1);
//    const dim3 grid( ceil(num_list/num_threads), 1, 1 );
//
//    cudaFuncSetCacheConfig( PackUnpackGPU::unpack_gpu<T>,    cudaFuncCachePreferL1 );
//
//    PackUnpackGPU::unpack_gpu<T> <<< grid, block >>> (list, val, buff, num_list);
//    cudaDeviceSynchronize();
    unpack_device(list, val, buff, num_list);
#endif
}


}; // PackUnpack //


namespace  PackUnpackCPU {


template <typename T>
void pack_cpu(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    )
{
    FuncLoop::Loop1d loop1d(0, num_list);

    loop1d.for_each_omp([&list, &val, &buff, num_list](int i) {
        PackUnpackKernel::pack_kernel(i, list, val, buff);
        });
}


template <typename T>
void unpack_cpu(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    )
{
    FuncLoop::Loop1d loop1d(0, num_list);

    loop1d.for_each_omp([&list, &val, &buff, num_list](int i) {
        PackUnpackKernel::unpack_kernel(i, list, val, buff);
        });
}


}; // PackUnpackCPU //


#ifdef USE_NVCC
namespace  PackUnpackGPU {


template <typename T>
__global__
void pack_gpu(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    )
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i >= num_list) { return; }

    PackUnpackKernel::pack_kernel(i, list, val, buff);
}


template <typename T>
__global__
void unpack_gpu(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    )
{
    const int i = threadIdx.x + blockDim.x*blockIdx.x;
    if (i >= num_list) { return; }

    PackUnpackKernel::unpack_kernel(i, list, val, buff);
}


}; // PackUnpackGPU //
#endif
