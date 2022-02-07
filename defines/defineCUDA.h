#pragma once
#ifndef DEFINECUDA_H_
#define DEFINECUDA_H_


#include "defineCal.h"


#ifdef GPU_CALCULATION__
  #if defined(USE_NVCC)
    #include <cuda.h>
  #elif defined(ENABLE_HIP)
    #include <hip/hip_runtime.h>
  #endif
#endif


#endif
