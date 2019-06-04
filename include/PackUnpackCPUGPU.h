#pragma once
#ifndef PACKUNPACKGPU_H_
#define PACKUNPACKGPU_H_


#include "defineCal.h"
#include "defineCUDA.h"


namespace  PackUnpack {


template <typename T>
void pack(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    );


template <typename T>
void unpack(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    );


}; // PackUnpack //


namespace  PackUnpackCPU {


template <typename T>
void pack_cpu(
    const int* list,
    const T*   val,
          T*   buff,
    const int  num_list
    );


template <typename T>
void unpack_cpu(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    );


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
    );


template <typename T>
__global__
void unpack_gpu(
    const int* list,
          T*   val,
    const T*   buff,
    const int  num_list
    );


}; // PackUnpackGPU //
#endif


#include "PackUnpackCPUGPU.hpp"


#endif
