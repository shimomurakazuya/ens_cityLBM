#pragma once
#ifndef FUNCLBMKERNEL_H_
#define FUNCLBMKERNEL_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include "definePrecision.h"
#include "defineLBM.h"
#include "Array3D.h"


namespace  FuncLBMKernel {


inline
__HOST__ __DEVICE__
void  stream_collision_sgs(
    const int    ids[],
    const int    ids_lbm_base[],
    const real*  f,
          real*  fn,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   kvis,
    const real   c_ref,
    const real   dx,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  stream_collision_cumulant_sgs_wo_wall(
    const int    ids[],
    const int    ids_lbm_base[],
    const real*  f,
          real*  fn,
          real*  sgs_vis,
    const real*  T,
    const real*  Tn,
    const real*  u,
    const real*  v,
    const real*  w,
    const real   kvis,
    const real   c_ref,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  stream_collision_cumulant_sgs(
    const int    ids[],
    const int    ids_lbm_base[],
    const real*  f,
          real*  fn,
          real*  sgs_vis,
    const real*  T,
    const real*  Tn,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   kvis,
    const real   c_ref,
    const real   dt
    );


inline
__HOST__ __DEVICE__
void  to_euler_variables(
    const int    id,
    const int    id_lbm_base,
    const real*  f,
          real*  rho,
          real*  u,
          real*  v,
          real*  w
    );


inline
__HOST__ __DEVICE__
void  flbm_in_wall(
    const int    id,
    const int    id_lbm_base,
    const real*  f,
          real*  fn,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs
    );


inline
__HOST__ __DEVICE__
void  flbm_inflow_outflow(
    const int    ids[],
    const int    ids_lbm_base[],
    const real*  f,
          real*  fn,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
    // solid wall //
    const real*  lv_obs,
    const real*  rho_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight,
    const real*  viscosity_weight
    );


inline
__HOST__ __DEVICE__
void  flbm_data_assimilation2(
    const real   time_weight0,
    const int    ids[],
    const int    ids_lbm_base[],
    const real*  f,
          real*  fn,
          real*  Tn,
          real*  scalarn,
    // solid wall //
    const real*  lv_obs,
    const real*  rho_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real*  T_obs,
    const real*  scalar_obs,
    // new //
    const real*  rhon_obs,
    const real*  un_obs,
    const real*  vn_obs,
    const real*  wn_obs,
    const real*  Tn_obs,
    const real*  scalarn_obs,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight,
    const real*  viscosity_weight
    );


inline
__HOST__ __DEVICE__
real  average(
    const int    ids[],
    const real*  val
    );


inline
__HOST__ __DEVICE__
real  average(
    const int    ids[],
    const real*  val,
    const real*  lv_obs
    );


inline
__HOST__ __DEVICE__
real  average3x3x3(
    const int    ids[],
    const real*  val
    );


inline
__HOST__ __DEVICE__
void  average_uvw(
          real&  ua,
          real&  va,
          real&  wa,
    const int    ids[],
    const real*  u,
    const real*  v,
    const real*  w
    );


};


#include "FuncLBMKernel.hpp"


#endif
