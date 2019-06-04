#include "PackUnpackKernel.h"


namespace  PackUnpackKernel {


template <typename T>
__HOST__ __DEVICE__
void pack_kernel(
    const int  i, // buff[i] //
    const int* list,
    const T*   val,
          T*   buff
    )
{
    buff[i] = val[list[i]];
}


template <typename T>
__HOST__ __DEVICE__
void unpack_kernel(
    const int  i, // buff[i] //
    const int* list,
          T*   val,
    const T*   buff
    )
{
    val[list[i]] = buff[i];
}


};
