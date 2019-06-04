#pragma once
#ifndef FUNCSCALARADVECTION_H_
#define FUNCSCALARADVECTION_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include "definePrecision.h"


namespace  FuncScalarAdvection {


inline
__HOST__ __DEVICE__
void  scalar_advection(
    const int    ids[],
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real   dx,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  scalar_advection(
    const int    ids[],
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real*  sc_scalar,
    const real   dx,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  scalar_advection_high_order(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[],
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  sc_scalar,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  scalar_advection_2nd(
    const int    ids[],
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real   dx,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  scalar_advection_2nd(
    const int    ids[],
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real*  sc_scalar,
    const real   dx,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  scalar_advection_3rd(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[],
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  sc_scalar,
    const real   dt
    );


inline
__HOST__ __DEVICE__
real  dfxv_1st(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    );


inline
__HOST__ __DEVICE__
real  dfxv_2nd(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    );


inline
__HOST__ __DEVICE__
real  dfxv_3rd(
    const real  f5[],
    const real  u5[],
    const real  lv5[]
    );


inline
__HOST__ __DEVICE__
real f_average(
    const int    ids[],
    const real* scalar,
    const real* lv_obs
    );


inline
__HOST__ __DEVICE__
real f_average(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[],
    const real*  scalar,
    const real*  lv_obs
    );


inline
__HOST__ __DEVICE__
real superbee(const real r);


inline
__HOST__ __DEVICE__
real  dfxx(
    const real  f3[],
    const real  lv3[],
    const real  coef_vis
    );


};


#include "FuncScalarAdvection.hpp"


#endif
