#pragma once
#ifndef FUNCLBM_H_
#define FUNCLBM_H_


#include "defineCUDA.h"
#include <cmath>
#include "definePrecision.h"


namespace  FuncLBM {


inline
__HOST__ __DEVICE__
real  drho_lbm(const real dp, const real c_ref);

inline
__HOST__ __DEVICE__
real  dp_lbm(const real drho, const real c_ref);


inline
__HOST__ __DEVICE__
real  rho_lbm(const real  fs[]);

inline
__HOST__ __DEVICE__
void  momentum_lbm(const real fs[], real& rho, real& mu, real& mv, real& mw);

inline
__HOST__ __DEVICE__
void  velocity_lbm(const real fs[], real& rho, real&  u, real&  v, real&  w);

//void  velocity_lbm(const real* f, const int id_lbm_base, real& rho, real&  u, real&  v, real&  w);

inline
__HOST__ __DEVICE__
real  kvis_lbm(const real kvis, const real c_ref, const real dt);

inline
__HOST__ __DEVICE__
real  relaxation_time(const real kvis_lbm);


inline
__HOST__ __DEVICE__
real
indexed_vel(
    const real  u,
    const real  v,
    const real  w,
    const int   i,
    const int   j,
    const int   k
    );


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
indexed_vel(
    const real  u,
    const real  v,
    const real  w
    );


inline
__HOST__ __DEVICE__
real
f_streaming(
    const real* f,
    const int   ids_lbm_base[],
    const int   iv, // read lbm velocity //
    const int   jv,
    const int   kv
    );


template <int iv, int jv, int kv>
__HOST__ __DEVICE__
real
f_streaming(
    const real* f,
    const int   ids_lbm_base[]
    );


inline
__HOST__ __DEVICE__
real
f_streaming(
    const real* f,
    const int   i,
    const int   j,
    const int   k,
    const int   offsets_lbm[],
    const int   iv,
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
real
f_bounce_back(
    const real* f,
    const int   ids[],
    const int   ids_lbm_base[],
    const int   iv,
    const int   jv,
    const int   kv,
    const real* lv_obs,
    const real* u_obs,
    const real* v_obs,
    const real* w_obs
    );


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
f_bounce_back(
    const real* f,
    const int   ids[],
    const int   ids_lbm_base[],
    const real* lv_obs,
    const real* u_obs,
    const real* v_obs,
    const real* w_obs
    );


inline
__HOST__ __DEVICE__
real
_f_bounce_back(
    const real* f,
    const int   ids[],
    const int   ids_lbm_base[],
    const int   iv,
    const int   jv,
    const int   kv,
    const real* lv_obs,
    const real* u_obs,
    const real* v_obs,
    const real* w_obs
    );


inline
__HOST__ __DEVICE__
real
factor_feq_D3Q27(
    const int   iv,
    const int   jv,
    const int   kv
    );


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
factor_feq_D3Q27();


inline
__HOST__ __DEVICE__
void
feq_D3Q27(
          real  fs[],
    const real  rho,
    const real  u,
    const real  v,
    const real  w
    );


inline
__HOST__ __DEVICE__
real
feq_D3Q27(
    const real  rho,
    const real  u,
    const real  v,
    const real  w,
    const int   i,
    const int   j,
    const int   k
    );


inline
__HOST__ __DEVICE__
real
force_x_D3Q27(
    const real  force,  // mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
real
force_y_D3Q27(
    const real  force,  // mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
real
force_z_D3Q27(
    const real  force,  // mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
real
force_D3Q27(
    const real  force[],  // force_xyz[3] : mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
real
force_D3Q27(
    const real  force[],  // force_xyz[3] : mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv,
    const real  u,
    const real  v,
    const real  w
    );


inline
__HOST__ __DEVICE__
real
force_on_boundary (
    const real  rb,
    const real  ub,
    const real  vb,
    const real  wb,
    const int   iv,
    const int   jv,
    const int   kv
    );


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
force_on_boundary (
    const real  rb,
    const real  ub,
    const real  vb,
    const real  wb
    );


inline
__HOST__ __DEVICE__
real
f_int(
    const real f0,
    const real f1,
    const real l0,
    const real l1
    );


};


#include "FuncLBM.hpp"


#endif
