#pragma once
#ifndef DEFINE_PRECISION_H_
#define DEFINE_PRECISION_H_


typedef unsigned char BYTE;


// CUDA
#ifdef USE_NVCC
  #include <cuda.h>
  #include <cuda_fp16.h>
  #include <cuda_bf16.h>
  
  #define FLOAT16_CAL
//  #define BFLOAT16_CAL
  
  #if   defined(FLOAT16_CAL)
    typedef half        fp16;
  #elif defined(BFLOAT16_CAL)
    typedef __nv_bfloat16 fp16;
  #endif
#endif

#define SINGLE_PRECISION_CALCULATION
//#define   DOUBLE_PRECISION_CALCULATION


// single precision
#ifdef SINGLE_PRECISION_CALCULATION
    typedef float   real;
    #define MFLOAT  MPI_FLOAT
#endif


// double precision
#ifdef DOUBLE_PRECISION_CALCULATION
    typedef double  real;
    #define MFLOAT  MPI_DOUBLE
#endif


#ifdef USE_NVCC
  #if   defined(FLOAT16_CAL)
    #define FP64toFP16(VAL) __double2half( (VAL) )
    #define FP32toFP16(VAL) __float2half( (VAL) )
    #define FP16toFP32(VAL) __half2float( (VAL) )
  #elif defined(BFLOAT16_CAL)
    #define FP64toFP16(VAL) __double2bfloat16( (VAL) )
    #define FP32toFP16(VAL) __float2bfloat16( (VAL) )
    #define FP16toFP32(VAL) __bfloat162float( (VAL) )
  #endif
#endif


#endif
