#pragma once
#ifndef FUNCLBMKERNEL_H_
#define FUNCLBMKERNEL_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include "definePrecision.h"
#include "defineLBM.h"
#include "Array3D.h"

#include "defineAMR.h"
#include "defineSGS.h"
#include "FuncLBM.h"
#include "FuncLoop.h"
#include "FuncObj.h"
#include "FuncLBMSGS.h"
#include "FuncThermalConvection.h"
#include "FuncSGS.h"
#include "defineFluidProperty.h"

#include "Index.h"
#include "IndexLBM.h"
#include "FuncCumulantLBM.h"


namespace  FuncLBMKernel {


// + sgs model //
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
    const real*  lv_obj,
    const real*  u_obj,
    const real*  v_obj,
    const real*  w_obj,
    const real   kvis,
    const real   c_ref,
    const real   dx,
    const real   dt
    )
{
    // fluid or wall //
    if ( FuncObj::is_obj( lv_obj[Index::id0(ids)] ) ) { return; }


    // streaming //
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                // lbm //
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                // ns //
                const int  id_up = Index::id(ids, -iv,-jv,-kv);

                fs[idv_leaf] = ( FuncObj::is_fluid( lv_obj[id_up] ) ) ? FuncLBM::f_streaming  (f, ids_lbm_base, iv,jv,kv) :
                                                                        FuncLBM::f_bounce_back(f, ids, ids_lbm_base, iv,jv,kv, lv_obj, u_obj,v_obj,w_obj);
            }
        }
    }


    // collision //
    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    real  fs_deq[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                fs_deq[idv_leaf] = fs[idv_leaf] - FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);
            }
        }
    }


    const real  Fcs_sgs = FuncSGS::get_Fcs(u,v,w, ids);
    const real  Csgs    = defineSGS::COEF_SGS;

    const real  kvis_lbm = FuncLBM::kvis_lbm(kvis, c_ref, dt);
#ifndef NO_SGS
    const real  sgs_vis_lbm = FuncLBMSGS::sgs_viscosity_D3Q27(Fcs_sgs, Csgs, rhos, kvis_lbm, fs_deq);
#elif NO_SGS
    const real  sgs_vis_lbm = 0.0;
#endif

    // viscosity //
    const real vis_total = kvis_lbm + sgs_vis_lbm;

    const real tau   = ((real)3.0*vis_total + (real)0.5);
    const real omega = (real)1.0 / tau;


    // update //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(ids_lbm_base, iv,jv,kv);

                fn[id_lbm] = fs[idv_leaf] - omega*fs_deq[idv_leaf];
            }
        }
    }

}


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
    const real*  pad_obj,
    const real   kvis,
    const real   c_ref,
    const real   dt
    )
{
    // streaming //
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                // lbm //
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                fs[idv_leaf] = FuncLBM::f_streaming  (f, ids_lbm_base, iv,jv,kv);
            }
        }
    }


    // collision //
    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // Boussinesq approximation //
//    const real Ts = fluid_property::Temperature_bouyancy;
    const real Ts = T[Index::id0(ids)];
//    const real Ts = (T[Index::id0(ids)] + average( ids, T ))*(real)0.5;
//    const real Ts = Tn[Index::id0(ids)];
//    const real Ts = (T[Index::id0(ids)] + Tn[Index::id0(ids)])*(real)0.5;
    const real drhog = -FuncThermalConvection::T_buoyancy( Ts ) * air_property::beta;
#ifndef NO_BUOYANCY 
    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;
#elif NO_BUOYANCY 
    const real gravity_force  = 0.0;
#endif
//    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref + fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;

    const real vel_mag = sqrtf(us*us + vs*vs + ws*ws)*c_ref;
    const real Cd_pad  = air_property::Cd_pad;
    const real pad = pad_obj[Index::id0(ids)];
    const real pad_force[] = { -rhos*Cd_pad*pad*us*vel_mag*dt, -rhos*Cd_pad*pad*vs*vel_mag*dt, -rhos*Cd_pad*pad*ws*vel_mag*dt };

//    const real force_half[] = { 0.0, 0.0, gravity_force*(real)0.5 }; // z : gravity
//    const real force_all[] = { 0.0, 0.0, gravity_force }; // z : gravity
    const real force_all[] = { pad_force[0], pad_force[1], pad_force[2]+gravity_force }; // z : gravity
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

//                fs[idv_leaf] += FuncLBM::force_D3Q27(force_half, iv, jv, kv, us, vs, ws);
                fs[idv_leaf] += FuncLBM::force_D3Q27(force_all, iv, jv, kv, us, vs, ws);
            }
        }
    }
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // sgs //
    real  fs_deq[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                fs_deq[idv_leaf] = fs[idv_leaf] - FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);
            }
        }
    }
    const real  Fcs_sgs = FuncSGS::get_Fcs(u,v,w, ids);
    const real  Csgs    = defineSGS::COEF_SGS;

    const real  kvis_lbm = FuncLBM::kvis_lbm(kvis, c_ref, dt);
#ifndef NO_SGS
    const real  sgs_vis_lbm = FuncLBMSGS::sgs_viscosity_D3Q27(Fcs_sgs, Csgs, rhos, kvis_lbm, fs_deq);
#elif NO_SGS
        const real  sgs_vis_lbm = 0.0;
#endif

    sgs_vis[Index::id0(ids)] = sgs_vis_lbm * (c_ref*c_ref*dt);
    // sgs //

    // viscosity //
    const real vis_total = kvis_lbm + sgs_vis_lbm;
    const real tau   = FuncLBM::relaxation_time(vis_total);
    const real omega = (real)1.0 / tau;

    real  fIs[27];
    FuncCumulantLBM::fs_cumulant_lbm(fIs, fs, omega, rhos, us,vs,ws);

    // update //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(ids_lbm_base, iv,jv,kv);

//                fn[id_lbm] = fIs[idv_leaf] + FuncLBM::force_D3Q27(force_half, iv, jv, kv, us, vs, ws);
                fn[id_lbm] = fIs[idv_leaf];
            }
        }
    }
}


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
    const real*  lv_obj,
    const real*  u_obj,
    const real*  v_obj,
    const real*  w_obj,
    const real*  pad_obj,
    const real   kvis,
    const real   c_ref,
    const real   dt
    )
{
    // fluid or wall //
    if ( FuncObj::is_obj( lv_obj[Index::id0(ids)] ) ) { return; }


    // streaming //
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                // lbm //
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                // ns //
                const int  id_up = Index::id(ids, -iv,-jv,-kv);

                fs[idv_leaf] = ( FuncObj::is_fluid( lv_obj[id_up] ) ) ? FuncLBM::f_streaming  (f, ids_lbm_base, iv,jv,kv) :
                                                                        FuncLBM::f_bounce_back(f, ids, ids_lbm_base, iv,jv,kv, lv_obj, u_obj,v_obj,w_obj);
            }
        }
    }


    // collision //
    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // Boussinesq approximation //
//    const real Ts = fluid_property::Temperature_bouyancy;
    const real Ts = T[Index::id0(ids)];
//    const real Ts = (T[Index::id0(ids)] + average( ids, T ))*(real)0.5;
//    const real Ts = Tn[Index::id0(ids)];
//    const real Ts = (T[Index::id0(ids)] + Tn[Index::id0(ids)])*(real)0.5;
    const real drhog = -FuncThermalConvection::T_buoyancy( Ts ) * air_property::beta;
#ifndef NO_BUOYANCY 
    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;
#elif NO_BUOYANCY 
    const real gravity_force  = 0.0;
#endif
//    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref + fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;

    const real vel_mag = sqrtf(us*us + vs*vs + ws*ws)*c_ref;
    const real Cd_pad  = air_property::Cd_pad;
    const real pad = pad_obj[Index::id0(ids)];
    const real pad_force[] = { -rhos*Cd_pad*pad*us*vel_mag*dt, -rhos*Cd_pad*pad*vs*vel_mag*dt, -rhos*Cd_pad*pad*ws*vel_mag*dt };

//    const real force_half[] = { 0.0, 0.0, gravity_force*(real)0.5 }; // z : gravity
//    const real force_all[] = { 0.0, 0.0, gravity_force }; // z : gravity
    const real force_all[] = { pad_force[0], pad_force[1], pad_force[2]+gravity_force }; // z : gravity
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

//                fs[idv_leaf] += FuncLBM::force_D3Q27(force_half, iv, jv, kv, us, vs, ws);
                fs[idv_leaf] += FuncLBM::force_D3Q27(force_all, iv, jv, kv, us, vs, ws);
            }
        }
    }
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // sgs //
    real  fs_deq[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                fs_deq[idv_leaf] = fs[idv_leaf] - FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);
            }
        }
    }
    const real  Fcs_sgs = FuncSGS::get_Fcs(u,v,w, ids);
    const real  Csgs    = defineSGS::COEF_SGS;

    const real  kvis_lbm = FuncLBM::kvis_lbm(kvis, c_ref, dt);
#ifndef NO_SGS
    const real  sgs_vis_lbm = FuncLBMSGS::sgs_viscosity_D3Q27(Fcs_sgs, Csgs, rhos, kvis_lbm, fs_deq);
#elif NO_SGS
    const real  sgs_vis_lbm = 0.0;
#endif
    
    sgs_vis[Index::id0(ids)] = sgs_vis_lbm * (c_ref*c_ref*dt);
    // sgs //

    // viscosity //
    const real vis_total = kvis_lbm + sgs_vis_lbm;
    const real tau   = FuncLBM::relaxation_time(vis_total);
    const real omega = (real)1.0 / tau;

    real  fIs[27];
    FuncCumulantLBM::fs_cumulant_lbm(fIs, fs, omega, rhos, us,vs,ws);

    // update //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(ids_lbm_base, iv,jv,kv);

//                fn[id_lbm] = fIs[idv_leaf] + FuncLBM::force_D3Q27(force_half, iv, jv, kv, us, vs, ws);
                fn[id_lbm] = fIs[idv_leaf];
            }
        }
    }

}


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
    )
{
#if 1
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);

                fs[idv_leaf] = f[id_lbm];
            }
        }
    }
#else
    const real  fs[27] = {
        f[IndexLBM::id0<-1,-1,-1>(id_lbm_base)],
        f[IndexLBM::id0< 0,-1,-1>(id_lbm_base)],
        f[IndexLBM::id0< 1,-1,-1>(id_lbm_base)],
        f[IndexLBM::id0<-1, 0,-1>(id_lbm_base)],
        f[IndexLBM::id0< 0, 0,-1>(id_lbm_base)],
        f[IndexLBM::id0< 1, 0,-1>(id_lbm_base)],
        f[IndexLBM::id0<-1, 1,-1>(id_lbm_base)],
        f[IndexLBM::id0< 0, 1,-1>(id_lbm_base)],
        f[IndexLBM::id0< 1, 1,-1>(id_lbm_base)],
        f[IndexLBM::id0<-1,-1, 0>(id_lbm_base)],
        f[IndexLBM::id0< 0,-1, 0>(id_lbm_base)],
        f[IndexLBM::id0< 1,-1, 0>(id_lbm_base)],
        f[IndexLBM::id0<-1, 0, 0>(id_lbm_base)],
        f[IndexLBM::id0< 0, 0, 0>(id_lbm_base)],
        f[IndexLBM::id0< 1, 0, 0>(id_lbm_base)],
        f[IndexLBM::id0<-1, 1, 0>(id_lbm_base)],
        f[IndexLBM::id0< 0, 1, 0>(id_lbm_base)],
        f[IndexLBM::id0< 1, 1, 0>(id_lbm_base)],
        f[IndexLBM::id0<-1,-1, 1>(id_lbm_base)],
        f[IndexLBM::id0< 0,-1, 1>(id_lbm_base)],
        f[IndexLBM::id0< 1,-1, 1>(id_lbm_base)],
        f[IndexLBM::id0<-1, 0, 1>(id_lbm_base)],
        f[IndexLBM::id0< 0, 0, 1>(id_lbm_base)],
        f[IndexLBM::id0< 1, 0, 1>(id_lbm_base)],
        f[IndexLBM::id0<-1, 1, 1>(id_lbm_base)],
        f[IndexLBM::id0< 0, 1, 1>(id_lbm_base)],
        f[IndexLBM::id0< 1, 1, 1>(id_lbm_base)],
    };
#endif


    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // update //
    rho[id] = rhos;
    u  [id] = us;
    v  [id] = vs;
    w  [id] = ws;
}


inline
__HOST__ __DEVICE__
void  to_euler_variables_with_time_average(
    const int    id,
    const int    id_lbm_base,
    const real*  f,
          real*  rho,
          real*  u,
          real*  v,
          real*  w,
    const real*  T,
          real*  u_mean_1min,
          real*  v_mean_1min,
          real*  w_mean_1min,
          real*  T_mean_1min,
          real*  vel2_fluc_1min,
//          real*  uu_fluc_1min,
//          real*  vv_fluc_1min,
//          real*  ww_fluc_1min,
    const real   time_1min,
    const real   dt
    )
{
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);

                fs[idv_leaf] = f[id_lbm];
            }
        }
    }

    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // update //
    rho[id] = rhos;
    u  [id] = us;
    v  [id] = vs;
    w  [id] = ws;

    // time average
    const real Ts = T[id];

    const real Tsec_1min = 60.0 * time_1min;
    const real coef_1min = dt / Tsec_1min;

    const real un_mean_1min = coef_1min*us + (1.0-coef_1min)*u_mean_1min[id];
    const real vn_mean_1min = coef_1min*vs + (1.0-coef_1min)*v_mean_1min[id];
    const real wn_mean_1min = coef_1min*ws + (1.0-coef_1min)*w_mean_1min[id];
    const real Tn_mean_1min = coef_1min*Ts + (1.0-coef_1min)*T_mean_1min[id];

    const real _vel2_fluc_1min =   pow(us-u_mean_1min[id], 2)
                                 + pow(vs-v_mean_1min[id], 2)
                                 + pow(ws-w_mean_1min[id], 2);

    const real vel2n_fluc_1min = coef_1min*_vel2_fluc_1min + (1.0-coef_1min)*vel2_fluc_1min[id];
//    const real uun_fluc_1min = coef_1min*pow(us-u_mean_1min[id], 2) + (1.0-coef_1min)*uu_fluc_1min[id];
//    const real vvn_fluc_1min = coef_1min*pow(vs-v_mean_1min[id], 2) + (1.0-coef_1min)*vv_fluc_1min[id];
//    const real wwn_fluc_1min = coef_1min*pow(ws-w_mean_1min[id], 2) + (1.0-coef_1min)*ww_fluc_1min[id];

    u_mean_1min[id] = un_mean_1min;
    v_mean_1min[id] = vn_mean_1min;
    w_mean_1min[id] = wn_mean_1min;
    T_mean_1min[id] = Tn_mean_1min;

    vel2_fluc_1min[id] = vel2n_fluc_1min;
//    uu_fluc_1min[id] = uun_fluc_1min;
//    vv_fluc_1min[id] = vvn_fluc_1min;
//    ww_fluc_1min[id] = wwn_fluc_1min;
}


inline
__HOST__ __DEVICE__
void  flbm_in_wall(
    const int    id,
    const int    id_lbm_base,
          real*  fn,
    const real*  lv_obj,
    const real*  u_obj,
    const real*  v_obj,
    const real*  w_obj
    )
{
    // fluid or wall //
    if ( FuncObj::is_fluid( lv_obj[id] ) ) { return; }

    constexpr real  rhob = (real)1.0;
    const real  ub = u_obj[id];
    const real  vb = v_obj[id];
    const real  wb = w_obj[id];
//    const real  ub = 0.0;
//    const real  vb = 0.0;
//    const real  wb = 0.0;

    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);

                fn[id_lbm] = FuncLBM::feq_D3Q27(rhob, ub, vb, wb, iv, jv, kv);
            }
        }
    }

}


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
          real*  Tn,
    // solid wall //
    const real*  lv_obj,
    const real*  rho_obj,
    const real*  u_obj,
    const real*  v_obj,
    const real*  w_obj,
    const real*  T_obj,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight
    )
{
    const int  id = Index::id0(ids);
    if ( FuncObj::is_obj( lv_obj[id] ) ) { return; }
//    if ( bcTypes_f[id] == BCTypes::CalRegion ) { return; }
    if ( bcTypes_f[id] != BCTypes::BCRegion  ) { return; }

    const real rhob = rho_obj[id];
    const real ub   = u_obj  [id];
    const real vb   = v_obj  [id];
    const real wb   = w_obj  [id];
    const real Tb   = T_obj  [id];

    const real dir_weight = dirichlet_weight[id];

    const real Torg = Tn[id];

    // lbm //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(ids_lbm_base, iv,jv,kv);

                const real forg    = fn[id_lbm];
                const real feq_obj = FuncLBM::feq_D3Q27(rhob, ub, vb, wb, iv, jv, kv);

                fn[id_lbm] = ( (real)1.0 - dir_weight )*forg + dir_weight*feq_obj;
            }
        }
    }

    Tn[id] = ( (real)1.0 - dir_weight )*Torg      + dir_weight*Tb;
}


inline
__HOST__ __DEVICE__
void  flbm_data_assimilation2(
    const real   time_weight0,
    const int    n_scalars,
    const int    nn_max,
    const int    id,
    const int    id_lbm_base,
    const real*  f,
          real*  fn,
          real*  Tn,
          real*  scalarn,
    // solid wall //
    const real*  lv_obj,
    const real*  rho_obj,
    const real*  u_obj,
    const real*  v_obj,
    const real*  w_obj,
    const real*  T_obj,
    const real*  scalar_obj,
    // new //
    const real*  rhon_obj,
    const real*  un_obj,
    const real*  vn_obj,
    const real*  wn_obj,
    const real*  Tn_obj,
    const real*  scalarn_obj,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight
    )
{
//    if ( FuncObj::is_obj( lv_obj[id] )       ) { return; }
//    if ( bcTypes_f[id] == BCTypes::CalRegion ) { return; }
    if ( bcTypes_f[id] != BCTypes::BCRegion  ) { return; }

    const real rhob = time_weight0*rho_obj[id] + ((real)1.0-time_weight0)*rhon_obj[id];
//    const real rhob = rho[id];
    const real ub      = time_weight0*u_obj     [id] + ((real)1.0-time_weight0)*un_obj     [id];
    const real vb      = time_weight0*v_obj     [id] + ((real)1.0-time_weight0)*vn_obj     [id];
    const real wb      = time_weight0*w_obj     [id] + ((real)1.0-time_weight0)*wn_obj     [id];
    const real Tb      = time_weight0*T_obj     [id] + ((real)1.0-time_weight0)*Tn_obj     [id];

    // values //
    const real Ts      = Tn[id];
    const real dir_weight = dirichlet_weight[id];

    // lbm //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);

                const real forg    = fn[id_lbm];
                const real feq_obj = FuncLBM::feq_D3Q27(rhob, ub, vb, wb, iv, jv, kv);
//                const real feq     = FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);

                fn[id_lbm] = ( (real)1.0 - dir_weight )*forg + dir_weight*feq_obj;

                // unstable //
//                const real f_obj = (forg + FuncLBM::force_D3Q27(force, iv, jv, kv, us, vs, ws)) * rhob/rhos;
//                fn[id_lbm] = ( (real)1.0 - dir_weight )*forg + dir_weight*f_obj;
            }
        }
    }

    // ns //
//    Tn[id] = Tb;
    Tn     [id] = ( (real)1.0 - dir_weight )*Ts      + dir_weight*Tb;
    for(int n=0; n<n_scalars; n++) {
        const int id_sc = id + nn_max*n;
        const real scalarb = time_weight0*scalar_obj[id_sc] + ((real)1.0-time_weight0)*scalarn_obj[id_sc];
        const real scalars = scalarn[id_sc];
        scalarn[id_sc] = ( (real)1.0 - dir_weight )*scalars + dir_weight*scalarb;
    }
}


inline
__HOST__ __DEVICE__
void  flbm_data_assimilation2_force(
    const real   time_weight0,
    const int    n_scalars,
    const int    nn_max,
    const int    id,
    const int    id_lbm_base,
    const real*  f,
          real*  fn,
    const real*  u,
    const real*  v,
    const real*  w,
          real*  Tn,
          real*  scalarn,
    // solid wall //
    const real*  lv_obj,
    const real*  rho_obj,
    const real*  u_obj,
    const real*  v_obj,
    const real*  w_obj,
    const real*  T_obj,
    const real*  scalar_obj,
    // new //
    const real*  rhon_obj,
    const real*  un_obj,
    const real*  vn_obj,
    const real*  wn_obj,
    const real*  Tn_obj,
    const real*  scalarn_obj,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight
    )
{
//    if ( FuncObj::is_obj( lv_obj[id] )       ) { return; }
//    if ( bcTypes_f[id] == BCTypes::CalRegion ) { return; }
    if ( bcTypes_f[id] != BCTypes::BCRegion  ) { return; }

    const real rhob = time_weight0*rho_obj[id] + ((real)1.0-time_weight0)*rhon_obj[id];
//    const real rhob = rho[id];
    const real ub = time_weight0*u_obj     [id] + ((real)1.0-time_weight0)*un_obj     [id];
    const real vb = time_weight0*v_obj     [id] + ((real)1.0-time_weight0)*vn_obj     [id];
    const real wb = time_weight0*w_obj     [id] + ((real)1.0-time_weight0)*wn_obj     [id];
    const real Tb = time_weight0*T_obj     [id] + ((real)1.0-time_weight0)*Tn_obj     [id];

    // values //
    const real us = u[id];
    const real vs = v[id];
    const real ws = w[id];
    const real Ts = Tn[id];
    const real dir_weight = dirichlet_weight[id];

    // lbm //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);

                const real forg    = fn[id_lbm];
//                const real feq_obj = FuncLBM::feq_D3Q27(rhob, ub, vb, wb, iv, jv, kv);

                const real force[] = { (ub - us)*dir_weight, (vb - vs)*dir_weight, (wb - ws)*dir_weight };
                fn[id_lbm] = forg + FuncLBM::force_D3Q27(force, iv, jv, kv, us, vs, ws);
            }
        }
    }

    // ns //
//    Tn[id] = Tb;
    Tn     [id] = ( (real)1.0 - dir_weight )*Ts      + dir_weight*Tb;
    for(int n=0; n<n_scalars; n++) {
        const int id_sc = id + nn_max*n;
        const real scalarb = time_weight0*scalar_obj[id_sc] + ((real)1.0-time_weight0)*scalarn_obj[id_sc];
        const real scalars = scalarn[id_sc];
        scalarn[id_sc] = ( (real)1.0 - dir_weight )*scalars + dir_weight*scalarb;
    }
}


inline
__HOST__ __DEVICE__
void  boundary_condition_outer_region_kernel(
    const real   time_weight0,
    const int    ids[],
          real*  scalarn,
    // solid wall //
    const real*  lv_obj,
    const real*  scalar_obj,
    const real*  scalarn_obj,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight
    )
{
    const int  id = Index::id0(ids);
//    if ( FuncObj::is_obj( lv_obj[id] )       ) { return; }
//    if ( bcTypes_f[id] == BCTypes::CalRegion ) { return; }
    if ( bcTypes_f[id] != BCTypes::BCRegion  ) { return; }

    const real scalarb = time_weight0*scalar_obj[id] + ((real)1.0-time_weight0)*scalarn_obj[id];

    // values //
    const real scalars = scalarn[id];
    const real dir_weight = dirichlet_weight[id];

    scalarn[id] = ( (real)1.0 - dir_weight )*scalars + dir_weight*scalarb;
}


inline
__HOST__ __DEVICE__
real  average(
    const int    ids[],
    const real*  val
    )
{
    return  ( (real)6.0*val[ Index::id(ids,  0,  0,  0) ]
             +          val[ Index::id(ids, -1,  0,  0) ] + val[ Index::id(ids,  1,  0,  0) ]
             +          val[ Index::id(ids,  0, -1,  0) ] + val[ Index::id(ids,  0,  1,  0) ]
             +          val[ Index::id(ids,  0,  0, -1) ] + val[ Index::id(ids,  0,  0,  1) ] ) / (real)12.0;
}


inline
__HOST__ __DEVICE__
real  average(
    const int    ids[],
    const real*  val,
    const real*  lv_obj
    )
{
    int   id_tmp;
    real  f = (real)0.0;
    int   count = 0;

    id_tmp = Index::id(ids, 0, 0, 0);
    f += val[ id_tmp ];   count++;

    id_tmp = Index::id(ids, -1, 0, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  1, 0, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, -1, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0,  1, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, 0, -1);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, 0,  1);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    return  ( f/(real)count );
}


inline
__HOST__ __DEVICE__
real  average3x3x3(
    const int    ids[],
    const real*  val
    )
{
    real  tmp = (real)0.0;
    for (int kk=0; kk<3; kk++) {
    for (int jj=0; jj<3; jj++) {
    for (int ii=0; ii<3; ii++) {
        tmp += val[ Index::id(ids, ii,jj,kk) ];
    }
    }
    }
    return  tmp/(real)27.0;
}


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
    )
{
    ua = average(ids, u);
    va = average(ids, v);
    wa = average(ids, w);
}


} // namespace FuncLBMKernel //


#endif
