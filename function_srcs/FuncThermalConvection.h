#pragma once
#ifndef FUNCTHERMALCONVECTION_H_
#define FUNCTHERMALCONVECTION_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include "defineObjBC.h"
#include "definePrecision.h"


namespace  FuncThermalConvection {


inline
__HOST__ __DEVICE__
void  cal_thermal_convection(
    const int    ids[],
    const real*  T,
          real*  T_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real*  T_obs,
    const real*  Tn_obs,
    const int*   bcTypes_T,
    const real   time_weight0,
    const real   coef_heatf,
    const real   dx,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  cal_thermal_convection_2nd(
    const int    ids[],
    const real*  T,
          real*  T_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real*  T_obs,
    const real*  Tn_obs,
    const int*   bcTypes_T,
    const real   time_weight0,
    const real   coef_heatf,
    const real   dx,
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
real  dfxv_adv_1st(
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
real  dfxv_adv_2nd(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    );


inline
__HOST__ __DEVICE__
real  dfxx(
    const real  f3[],
    const real  lv3[],
    const real  Tobs3[],
    const int   bcTypesT3[],
    const real  coef_heatf
    );


inline
__HOST__ __DEVICE__
real  rho_buoyancy(const real  T);

inline
__HOST__ __DEVICE__
real  T_buoyancy(const real  T);


inline
__HOST__ __DEVICE__
real  T_buoyancy(
    const real T,
    const real Txm, const real Txp,
    const real Tym, const real Typ,
    const real Tzm, const real Tzp);


inline
__HOST__ __DEVICE__
real f_average(
    const int    ids[],
    const real* val
    );


inline
__HOST__ __DEVICE__
real f_average_w_obs(
    const int    ids[],
    const real* val,
    const real* lv_obs
    );


inline
__HOST__ __DEVICE__
real f_average(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[],
    const real*  val,
    const real*  lv_obs
    );


inline
__HOST__ __DEVICE__
real f_time_weight(
    const real f0,
    const real fn,
    const real weight0
    );


};


#include "FuncThermalConvection.hpp"


#endif
