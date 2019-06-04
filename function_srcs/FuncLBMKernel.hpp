#include "FuncLBMKernel.h"
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
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   kvis,
    const real   c_ref,
    const real   dx,
    const real   dt
    )
{
    // fluid or wall //
    if ( FuncObj::is_obj( lv_obs[Index::id0(ids)] ) ) { return; }


    // streaming //
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                // lbm //
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                // ns //
                const int  id_up = Index::id(ids, -iv,-jv,-kv);

                fs[idv_leaf] = ( FuncObj::is_fluid( lv_obs[id_up] ) ) ? FuncLBM::f_streaming  (f, ids_lbm_base, iv,jv,kv) :
                                                                        FuncLBM::f_bounce_back(f, ids, ids_lbm_base, iv,jv,kv, lv_obs, u_obs,v_obs,w_obs);
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
    const real  sgs_vis_lbm = FuncLBMSGS::sgs_viscosity_D3Q27(Fcs_sgs, Csgs, rhos, kvis_lbm, fs_deq);
//    const real  sgs_vis_lbm = 0.0;


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
//    const real Ts = T[Index::id0(ids)];
//    const real Ts = (T[Index::id0(ids)] + average( ids, T ))*(real)0.5;
//    const real Ts = Tn[Index::id0(ids)];
    const real Ts = (T[Index::id0(ids)] + Tn[Index::id0(ids)])*(real)0.5;
    const real drhog = -FuncThermalConvection::T_buoyancy( Ts ) * air_property::beta;
    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;
//    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref + fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;

//    const real force_half[] = { 0.0, 0.0, gravity_force*(real)0.5 }; // z : gravity
    const real force_all[] = { 0.0, 0.0, gravity_force }; // z : gravity
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
    const real  sgs_vis_lbm = FuncLBMSGS::sgs_viscosity_D3Q27(Fcs_sgs, Csgs, rhos, kvis_lbm, fs_deq);

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
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   kvis,
    const real   c_ref,
    const real   dt
    )
{
    // fluid or wall //
    if ( FuncObj::is_obj( lv_obs[Index::id0(ids)] ) ) { return; }


    // streaming //
    real  fs[27];
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                // lbm //
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                // ns //
                const int  id_up = Index::id(ids, -iv,-jv,-kv);

                fs[idv_leaf] = ( FuncObj::is_fluid( lv_obs[id_up] ) ) ? FuncLBM::f_streaming  (f, ids_lbm_base, iv,jv,kv) :
                                                                        FuncLBM::f_bounce_back(f, ids, ids_lbm_base, iv,jv,kv, lv_obs, u_obs,v_obs,w_obs);
            }
        }
    }


    // collision //
    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    // Boussinesq approximation //
//    const real Ts = fluid_property::Temperature_bouyancy;
//    const real Ts = T[Index::id0(ids)];
//    const real Ts = (T[Index::id0(ids)] + average( ids, T ))*(real)0.5;
//    const real Ts = Tn[Index::id0(ids)];
    const real Ts = (T[Index::id0(ids)] + Tn[Index::id0(ids)])*(real)0.5;
    const real drhog = -FuncThermalConvection::T_buoyancy( Ts ) * air_property::beta;
    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;
//    const real gravity_force  = rhos * (fluid_property::GravitationalAcceleration/c_ref + fluid_property::GravitationalAcceleration/c_ref * drhog) * dt;

//    const real force_half[] = { 0.0, 0.0, gravity_force*(real)0.5 }; // z : gravity
    const real force_all[] = { 0.0, 0.0, gravity_force }; // z : gravity
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
    const real  sgs_vis_lbm = FuncLBMSGS::sgs_viscosity_D3Q27(Fcs_sgs, Csgs, rhos, kvis_lbm, fs_deq);

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
void  flbm_in_wall(
    const int    id,
    const int    id_lbm_base,
    const real*  f,
          real*  fn,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs
    )
{
    // fluid or wall //
    if ( FuncObj::is_fluid( lv_obs[id] ) ) { return; }

    constexpr real  rhob = (real)1.0;
    const real  ub = u_obs[id];
    const real  vb = v_obs[id];
    const real  wb = w_obs[id];
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
    )
{
    const int  id = Index::id0(ids);
    if ( FuncObj::is_obj( lv_obs[id] ) ) { return; }
    if ( bcTypes_f[id] == BCTypes::CalRegion ) { return; }

    const real rhob = rho_obs[id];
    const real ub   = u_obs  [id];
    const real vb   = v_obs  [id];
    const real wb   = w_obs  [id];

    const real rhos = rho[id];
    const real us   = u  [id];
    const real vs   = v  [id];
    const real ws   = w  [id];

    const real dir_weight = dirichlet_weight[id];
    const real vis_weight = viscosity_weight[id];

    // lbm //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(ids_lbm_base, iv,jv,kv);

                const real forg    = fn[id_lbm];
                const real feq_obs = FuncLBM::feq_D3Q27(rhob, ub, vb, wb, iv, jv, kv);
                const real feq     = FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);

                fn[id_lbm] = ( (real)1.0 - dir_weight - vis_weight )*forg + dir_weight*feq_obs + vis_weight*feq;
            }
        }
    }

}


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
    )
{
    const int  id = Index::id0(ids);
//    if ( FuncObj::is_obj( lv_obs[id] )       ) { return; }
//    if ( bcTypes_f[id] == BCTypes::CalRegion ) { return; }
    if ( bcTypes_f[id] != BCTypes::BCRegion  ) { return; }

    const real rhob = time_weight0*rho_obs[id] + ((real)1.0-time_weight0)*rhon_obs[id];
//    const real rhob = rho[id];
    const real ub      = time_weight0*u_obs     [id] + ((real)1.0-time_weight0)*un_obs     [id];
    const real vb      = time_weight0*v_obs     [id] + ((real)1.0-time_weight0)*vn_obs     [id];
    const real wb      = time_weight0*w_obs     [id] + ((real)1.0-time_weight0)*wn_obs     [id];
    const real Tb      = time_weight0*T_obs     [id] + ((real)1.0-time_weight0)*Tn_obs     [id];
    const real scalarb = time_weight0*scalar_obs[id] + ((real)1.0-time_weight0)*scalarn_obs[id];

    // values //
    const real Ts      = Tn[id];
    const real scalars = scalarn[id];
    const real dir_weight = dirichlet_weight[id];

//    real  fs[27];
//    for (int kv=-1; kv<=1; kv++) {
//        for (int jv=-1; jv<=1; jv++) {
//            for (int iv=-1; iv<=1; iv++) {
//                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
//                const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);
//
//                fs[idv_leaf] = f[id_lbm];
//            }
//        }
//    }
//    real  rhos, us, vs, ws;
//    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);
//    const real force[] = { ub-us, vb-vs, wb-ws };

    // lbm //
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                const int  id_lbm   = IndexLBM::id0(ids_lbm_base, iv,jv,kv);

                const real forg    = fn[id_lbm];
                const real feq_obs = FuncLBM::feq_D3Q27(rhob, ub, vb, wb, iv, jv, kv);
//                const real feq     = FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);

                fn[id_lbm] = ( (real)1.0 - dir_weight )*forg + dir_weight*feq_obs;

                // unstable //
//                const real f_obs = (forg + FuncLBM::force_D3Q27(force, iv, jv, kv, us, vs, ws)) * rhob/rhos;
//                fn[id_lbm] = ( (real)1.0 - dir_weight )*forg + dir_weight*f_obs;
            }
        }
    }

    // ns //
//    Tn[id] = Tb;
    Tn     [id] = ( (real)1.0 - dir_weight )*Ts      + dir_weight*Tb;
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
    const real*  lv_obs
    )
{
    int   id_tmp;
    real  f = (real)0.0;
    int   count = 0;

    id_tmp = Index::id(ids, 0, 0, 0);
    f += val[ id_tmp ];   count++;

    id_tmp = Index::id(ids, -1, 0, 0);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  1, 0, 0);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, -1, 0);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0,  1, 0);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, 0, -1);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, 0,  1);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

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


};
