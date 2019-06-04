#pragma once
#ifndef FUNCSGS_H_
#define FUNCSGS_H_


#include "defineCUDA.h"
#include <iostream>
#include <cmath>
#include "definePrecision.h"


namespace  FuncSGS {


// struct //
struct  TensorVel {
    real  s11, s22, s33;
    real  s12, s23, s13;
    real  w12, w23, w13;
};


struct  DerivativeVel {
    real  ux, uy, uz;
    real  vx, vy, vz;
    real  wx, wy, wz;
};


// functions //
inline
__HOST__ __DEVICE__
real  SS (const TensorVel&  tVel);

inline
__HOST__ __DEVICE__
real  WW (const TensorVel&  tVel);

inline
__HOST__ __DEVICE__
real  Fcs(const TensorVel&  tVel);


inline
__HOST__ __DEVICE__
real
get_Fcs(
    const real* u,
    const real* v,
    const real* w,
    const int   stencil[]
    );


inline
__HOST__ __DEVICE__
DerivativeVel
get_derivatives(
    const real* u,
    const real* v,
    const real* w,
    const int   stencil[]
    );


inline
__HOST__ __DEVICE__
TensorVel
get_tensor(
    const real* u,
    const real* v,
    const real* w,
    const int   stencil[]
    );


};


#include "FuncSGS.hpp"


#endif
