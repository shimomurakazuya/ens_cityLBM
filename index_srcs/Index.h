#pragma once
#ifndef INDEX_H_
#define INDEX_H_


#include "defineCUDA.h"
#include <iostream>
#include "defineCal.h"
#include "Array3D.h"


namespace  Index {


inline
__HOST__ __DEVICE__
int  idv(const int iv, const int jv, const int kv);


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int  idv();


inline
__HOST__ __DEVICE__
int  idv0();


inline
__HOST__ __DEVICE__
int  id(
    const int   ix,
    const int   iy,
    const int   iz
    );


inline
__HOST__ __DEVICE__
int id(
    const int  ids[],
    const int  in, // neighbor access //
    const int  jn,
    const int  kn
    );


template<int in, int jn, int kn>
__HOST__ __DEVICE__
int id(const int  ids[]);


inline
__HOST__ __DEVICE__
int id0(const int  ids[]);


inline
__HOST__ __DEVICE__
int  id0(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offset0
    );


inline
__HOST__ __DEVICE__
int  id0(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offsets3d[]
    );


inline
__HOST__ __DEVICE__
int  id(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offsets3d[]
    );


};


#include "Index.hpp"


#endif
