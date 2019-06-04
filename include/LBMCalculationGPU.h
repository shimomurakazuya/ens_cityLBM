#pragma once
#ifndef LBMCALCULATIONGPU_H_
#define LBMCALCULATIONGPU_H_


#include "defineCUDA.h"
#include "definePrecision.h"
#include "ValueObjLS.h"


#ifdef USE_NVCC
namespace  LBMCalculationGPU {


__global__
void  to_euler_variables_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets,
    const real*  f_lbm,
          real*  rho,
          real*  u,
          real*  v,
          real*  w
    );


__global__
void  stream_collision_sgs_gpu_wo_wall(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
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


__global__
void  stream_collision_sgs_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
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


__global__
void  scalar_advection_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
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


__global__
void  scalar_advection_w_source_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
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


__global__
void  temperature_advection_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
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


__global__
void  LBM_L2F_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  f_lbmF,
    const real*  f_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    );


__global__
void  LBM_L2FA_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  f_lbmF,
    const real*  f_lbm,
    const real*  fn_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    );


__global__
void  Val_L2F_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  scalarF,
    const real*  scalar
    );


__global__
void  Val_L2FA_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  scalarF,
    const real*  scalar,
    const real*  scalarn
    );


__global__
void  LBM_F2L_gpu(
    const int*   id_tasks,
    const int*   id_parent,
    const int*   mesh_offsets3x3x3,
    const int    lvF,
    const int*   node_pos_x,
    const int*   node_pos_y,
    const int*   node_pos_z,
          real*  f_lbmC,
    const real*  f_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    );


__global__
void  LBM_F2LA_gpu(
    const int*   id_tasks,
    const int*   id_parent,
    const int*   mesh_offsets3x3x3,
    const int    lvF,
    const int*   node_pos_x,
    const int*   node_pos_y,
    const int*   node_pos_z,
          real*  f_lbmC,
          real*  fn_lbmC,
    const real*  f_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    );


__global__
void  Val_F2L_gpu(
    const int*   id_tasks,
    const int*   id_parent,
    const int*   mesh_offsets3x3x3,
    const int    lvF,
    const int*   node_pos_x,
    const int*   node_pos_y,
    const int*   node_pos_z,
          real*  scalarC,
    const real*  scalar
    );


__global__
void  boundary_conditions_in_wall_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets,
    const real*  f_lbm,
          real*  fn_lbm,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   c_ref,
    const real   dx,
    const real   dt
    );


__global__
void  boundary_conditions_inflow_outflow_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
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
    const real*  viscosity_weight,
    const real   c_ref,
    const real   dx,
    const real   dt
    );


__global__
void  boundary_condition_data_assimilation_gpu(
    const real   time_weight0,
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
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
    const real*  rhon_obs,
    const real*  un_obs,
    const real*  vn_obs,
    const real*  wn_obs,
    const real*  Tn_obs,
    const real*  scalarn_obs,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight,
    const real*  viscosity_weight,
    const real   c_ref,
    const real   dx,
    const real   dt
    );


};
#endif


#endif
