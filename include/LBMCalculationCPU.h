#pragma once
#ifndef LBMCALCULATIONCPU_H_
#define LBMCALCULATIONCPU_H_


#include "defineCUDA.h"
#include "definePrecision.h"
#include "ValueObjLS.h"


namespace  LBMCalculationCPU {


void  to_euler_variables_cpu(
    const int    num_tasks,
    const int*   id_tasks,
    const int*   mesh_offsets,
    const real*  f_lbm,
          real*  rho,
          real*  u,
          real*  v,
          real*  w
    );


void  stream_collision_sgs_cpu_wo_wall(
    const int    num_tasks,
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


void  stream_collision_sgs_cpu(
    const int    num_tasks,
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


void  scalar_advection_cpu(
    const int    num_tasks,
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  sc_scalar,
    const real   dx,
    const real   dt
    );


void  temperature_advection_cpu(
    const int    num_tasks,
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


void  LBM_L2F_cpu(
    const int    num_tasks,
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


void  LBM_L2FA_cpu(
    const int    num_tasks,
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


void  LBM_L2FR_cpu(
    const int    num_tasks,
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


void  Val_L2F_cpu(
    const int    num_tasks,
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  scalarF,
    const real*  scalar
    );


void  Val_L2FA_cpu(
    const int    num_tasks,
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  scalarF,
    const real*  scalar,
    const real*  scalarn
    );


void  LBM_F2L_cpu(
    const int    num_tasks,
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


void  Val_F2L_cpu(
    const int    num_tasks,
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


void  boundary_conditions_in_wall_cpu(
    const int    num_tasks,
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


void  boundary_conditions_inflow_outflow_cpu(
    const int    num_tasks,
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


};


#endif
