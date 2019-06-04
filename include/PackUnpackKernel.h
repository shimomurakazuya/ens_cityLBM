#pragma once
#ifndef PACKUNPACKKERNEL_H_
#define PACKUNPACKKERNEL_H_


#include "defineCUDA.h"


namespace  PackUnpackKernel {


template <typename T>
__HOST__ __DEVICE__
void pack_kernel(
    const int  i, // buff[i] //
    const int* list,
    const T*   val,
          T*   buff
    );


template <typename T>
__HOST__ __DEVICE__
void unpack_kernel(
    const int  i, // buff[i] //
    const int* list,
          T*   val,
    const T*   buff
    );


};


#include "PackUnpackKernel.hpp"


#endif
