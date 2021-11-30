#pragma once
#ifndef FUNCSGS_H_
#define FUNCSGS_H_


#include "defineCUDA.h"
#include <iostream>
#include <cmath>
#include "definePrecision.h"
#include "Index.h"


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
real  SS(const TensorVel&  tVel)
{
    return               ( tVel.s11*tVel.s11 + tVel.s22*tVel.s22 + tVel.s33*tVel.s33 )
             + (real)2.0*( tVel.s12*tVel.s12 + tVel.s23*tVel.s23 + tVel.s13*tVel.s13 );
}


inline
__HOST__ __DEVICE__
real  WW(const TensorVel&  tVel)
{
    return     (real)2.0*( tVel.w12*tVel.w12 + tVel.w23*tVel.w23 + tVel.w13*tVel.w13 );
}


inline
__HOST__ __DEVICE__
real  Fcs(const TensorVel&  tVel)
{
    // SS, WW //
    const real  ss = SS(tVel);
    const real  ww = WW(tVel);
    const real  ep = (real)1.0e-16;

    return  fabs( ( ww-ss ) / ( ww+ss + ep ) );
}


inline
__HOST__ __DEVICE__
DerivativeVel
get_derivatives(
    const real* u,
    const real* v,
    const real* w,
    const int   stencil[]
    )
{
    return  DerivativeVel{
                ( u[stencil[Index::idv(1, 0, 0)]] - u[stencil[Index::idv(-1, 0, 0)]] ) * (real)0.5,
                ( u[stencil[Index::idv(0, 1, 0)]] - u[stencil[Index::idv(0, -1, 0)]] ) * (real)0.5,
                ( u[stencil[Index::idv(0, 0, 1)]] - u[stencil[Index::idv(0, 0, -1)]] ) * (real)0.5,

                ( v[stencil[Index::idv(1, 0, 0)]] - v[stencil[Index::idv(-1, 0, 0)]] ) * (real)0.5,
                ( v[stencil[Index::idv(0, 1, 0)]] - v[stencil[Index::idv(0, -1, 0)]] ) * (real)0.5,
                ( v[stencil[Index::idv(0, 0, 1)]] - v[stencil[Index::idv(0, 0, -1)]] ) * (real)0.5,

                ( w[stencil[Index::idv(1, 0, 0)]] - w[stencil[Index::idv(-1, 0, 0)]] ) * (real)0.5,
                ( w[stencil[Index::idv(0, 1, 0)]] - w[stencil[Index::idv(0, -1, 0)]] ) * (real)0.5,
                ( w[stencil[Index::idv(0, 0, 1)]] - w[stencil[Index::idv(0, 0, -1)]] ) * (real)0.5
                };
}


inline
__HOST__ __DEVICE__
TensorVel
get_tensor(
    const real* u,
    const real* v,
    const real* w,
    const int   stencil[]
    )
{
    const DerivativeVel  dVel = get_derivatives(u,v,w, stencil);

    return  TensorVel{
                dVel.ux,
                dVel.vy,
                dVel.wz,

                (dVel.vx + dVel.uy) * (real)0.5,
                (dVel.wy + dVel.vz) * (real)0.5,
                (dVel.uz + dVel.wx) * (real)0.5,

                (dVel.vx - dVel.uy) * (real)0.5,
                (dVel.wy - dVel.vz) * (real)0.5,
                (dVel.uz - dVel.wx) * (real)0.5
                };
}


inline
__HOST__ __DEVICE__
real
get_Fcs(
    const real* u,
    const real* v,
    const real* w,
    const int   stencil[]
    )
{
    return  Fcs( get_tensor(u,v,w, stencil) );
}


} // namespace FuncSGS


#endif
