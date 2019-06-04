#include "FuncLBM.h"
#include "defineAMR.h"
#include "defineSGS.h"
#include "FuncMath.h"
#include "FuncObj.h"

#include "Index.h"
#include "IndexLBM.h"


namespace  FuncLBM {


inline
__HOST__ __DEVICE__
real  drho_lbm(const real dp, const real c_ref)
{
    return  (real)3.0*dp/(c_ref*c_ref);
}


inline
__HOST__ __DEVICE__
real  dp_lbm(const real drho, const real c_ref)
{
    return drho*(c_ref*c_ref)/(real)3.0;
}


inline
__HOST__ __DEVICE__
real  rho_lbm(const real  fs[])
{
    return  fs[0] + fs[1] + fs[2] + fs[3] + fs[4] + fs[5] + fs[6] + fs[7] + fs[8] + fs[9] + fs[10] + fs[11] + fs[12] + fs[13] + fs[14] + fs[15] + fs[16] + fs[17] + fs[18] + fs[19] + fs[20] + fs[21] + fs[22] + fs[23] + fs[24] + fs[25] + fs[26];
}


inline
__HOST__ __DEVICE__
void  momentum_lbm(const real fs[], real& rho, real& mu, real& mv, real& mw)
{
    rho =  fs[0] + fs[1] + fs[2] + fs[3] + fs[4] + fs[5] + fs[6] + fs[7] + fs[8] + fs[9] + fs[10] + fs[11] + fs[12] + fs[13] + fs[14] + fs[15] + fs[16] + fs[17] + fs[18] + fs[19] + fs[20] + fs[21] + fs[22] + fs[23] + fs[24] + fs[25] + fs[26];

    mu  = -fs[0] + fs[2] - fs[3] + fs[5] - fs[6] + fs[8] - fs[9] + fs[11] - fs[12] + fs[14] - fs[15] + fs[17] - fs[18] + fs[20] - fs[21] + fs[23] - fs[24] + fs[26];
    mv  = -fs[0] - fs[1] - fs[2] + fs[6] + fs[7] + fs[8] - fs[9] - fs[10] - fs[11] + fs[15] + fs[16] + fs[17] - fs[18] - fs[19] - fs[20] + fs[24] + fs[25] + fs[26];
    mw  = -fs[0] - fs[1] - fs[2] - fs[3] - fs[4] - fs[5] - fs[6] - fs[7]  - fs[8]  + fs[18] + fs[19] + fs[20] + fs[21] + fs[22] + fs[23] + fs[24] + fs[25] + fs[26];
}


inline
__HOST__ __DEVICE__
void  velocity_lbm(const real fs[], real& rho, real&  u, real&  v, real&  w)
{
    real  mu, mv, mw;
    momentum_lbm(fs, rho, mu,mv,mw);
    const real  drho = (real)1.0/rho;

    u = mu*drho;
    v = mv*drho;
    w = mw*drho;
}


// slow //
//void  velocity_lbm(const real* f, const int id_lbm_base, real& rho, real&  u, real&  v, real&  w)
//{
//    rho = 0.0;
//    u = 0.0; v = 0.0; w = 0.0;
//    for (int kv=-1; kv<=1; kv++) {
//        for (int jv=-1; jv<=1; jv++) {
//            for (int iv=-1; iv<=1; iv++) {
//                const real  ftmp = f[id_lbm_base + IndexLBM::idv_lbm_leaf(iv, jv, kv)];
//
//                rho += ftmp;
//                u   += ftmp*iv;
//                v   += ftmp*jv;
//                w   += ftmp*kv;
//            }
//        }
//    }
//    u /= rho;
//    v /= rho;
//    w /= rho;
//}


inline
__HOST__ __DEVICE__
real  kvis_lbm(const real kvis, const real c_ref, const real dt)
{
    return kvis/(c_ref*c_ref*dt);
}


inline
__HOST__ __DEVICE__
real  relaxation_time(const real kvis_lbm)
{
    // tau //
    return  (real)3.0*kvis_lbm + (real)0.5;
}


inline
__HOST__ __DEVICE__
real
indexed_vel(
    const real  u,
    const real  v,
    const real  w,
    const int   iv,
    const int   jv,
    const int   kv
    )
{
    return  iv*u + jv*v + kv*w;
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
indexed_vel(
    const real  u,
    const real  v,
    const real  w
    )
{
    return  iv*u + jv*v + kv*w;
}


inline
__HOST__ __DEVICE__
real
f_streaming(
    const real* f,
    const int   ids_lbm_base[],
    const int   iv, // read lbm velocity //
    const int   jv,
    const int   kv
    )
{
    const int  id_lbm = IndexLBM::id(ids_lbm_base, -iv,-jv,-kv, iv,jv,kv);

    return  f[ id_lbm ];
}


template <int iv, int jv, int kv>
__HOST__ __DEVICE__
real
f_streaming(
    const real* f,
    const int   ids_lbm_base[]
    )
{
    const int  id_lbm = IndexLBM::id<-iv,-jv,-kv, iv,jv,kv>(ids_lbm_base);

    return  f[ id_lbm ];
}


inline
__HOST__ __DEVICE__
real
f_streaming(
    const real* f,
    const int   i,
    const int   j,
    const int   k,
    const int   offsets[],
    const int   iv, // read lbm velocity //
    const int   jv,
    const int   kv
    )
{
    const int  id_lbm = IndexLBM::id(i-iv,j-jv,k-kv, offsets, iv,jv,kv);

    return  f[ id_lbm ];
}


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
    )
{
    constexpr real  ep = (real)1.0e-7;

    const int   id0     = Index::id0( ids );
    const int   id_up   = Index::id ( ids, -iv,-jv,-kv);
    const int   id_down = Index::id ( ids,  iv, jv, kv);

    const real  fp_f  = f[ IndexLBM::id(ids_lbm_base, 0, 0, 0,  -iv,-jv,-kv ) ];
    const real  fm_f  = f[ IndexLBM::id(ids_lbm_base, 0, 0, 0,   iv, jv, kv ) ];
    const real  fp_ff = f[ IndexLBM::id(ids_lbm_base, iv,jv,kv, -iv,-jv,-kv ) ];

    // delta : 0 ~ 0.5 ~ 1.0 //
    const real  lv0     = lv_obs[id0];
    const real  lv_up   = lv_obs[id_up];
    const real  lv_down = lv_obs[id_down];

//    const real  delta   = fabs( lv0 ) / ( fabs(lv0) + fabs(lv_up) + ep );
    const real  delta =  ( FuncObj::is_fluid( lv_down ) )
                            ? fabs( lv0 ) / ( fabs(lv0) + fabs(lv_up) + ep )
                            : (real)0.5;

    // force //
//    const real  rb = (real)1.0; // density on boundary //
    constexpr real  rb = (real)1.0; // density on boundary //

    const real  ub = f_int( u_obs[id0], u_obs[id_up], lv0, lv_up );
    const real  vb = f_int( v_obs[id0], v_obs[id_up], lv0, lv_up );
    const real  wb = f_int( w_obs[id0], w_obs[id_up], lv0, lv_up );

    const real  force_bc = force_on_boundary(rb, ub,vb,wb, iv,jv,kv);

//    return  fp_ff;
    return  (delta < (real)0.5) ?
                ((real)2.0*delta*fp_f) + ((real)1.0 - (real)2.0*delta)*fp_ff + force_bc :
                ( fp_f + ((real)2.0*delta - (real)1.0)*fm_f + force_bc ) / ((real)2.0*delta);
}


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
    )
{
//    constexpr real  ep = (real)1.0e-6;

    const int   id0     = Index::id0( ids );
    const int   id_up   = Index::id<-iv,-jv,-kv>( ids );
    const int   id_down = Index::id< iv, jv, kv>( ids );

    const real  fp_f  = f[ IndexLBM::id<0, 0, 0,  -iv,-jv,-kv>(ids_lbm_base) ];
    const real  fm_f  = f[ IndexLBM::id<0, 0, 0,   iv, jv, kv>(ids_lbm_base) ];
    const real  fp_ff = f[ IndexLBM::id<iv,jv,kv, -iv,-jv,-kv>(ids_lbm_base) ];

    // delta : 0 ~ 0.5 ~ 1.0 //
    const real  lv0     = lv_obs[id0];
    const real  lv_up   = lv_obs[id_up];
    const real  lv_down = lv_obs[id_down];

//    const real  delta   = fabs( lv0 ) / ( fabs(lv0) + fabs(lv_up) + (real)1.0e-7 );
    const real  delta =  ( FuncObj::is_fluid( lv_down ) )
                            ? fabs( lv0 ) / ( fabs(lv0) + fabs(lv_up) + (real)1.0e-7 )
                            : (real)0.5;

    // force //
//    const real  rb = (real)1.0; // density on boundary //
    constexpr real  rb = (real)1.0; // density on boundary //

    const real  ub = f_int( u_obs[id0], u_obs[id_up], lv0, lv_up );
    const real  vb = f_int( v_obs[id0], v_obs[id_up], lv0, lv_up );
    const real  wb = f_int( w_obs[id0], w_obs[id_up], lv0, lv_up );

    const real  force_bc = force_on_boundary<iv,jv,kv>(rb, ub,vb,wb);

//    return  fp_ff;
    return  (delta < (real)0.5) ?
                ((real)2.0*delta*fp_f) + ((real)1.0 - (real)2.0*delta)*fp_ff + force_bc :
                ( fp_f + ((real)2.0*delta - (real)1.0)*fm_f + force_bc ) / ((real)2.0*delta);
}


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
    )
{
//    constexpr real  ep = (real)1.0e-6;

    const int   id0   = Index::id0( ids );
    const int   id_up = Index::id ( ids, -iv,-jv,-kv);

    const real  fp_f  = f[ IndexLBM::id(ids_lbm_base, 0, 0, 0,  -iv,-jv,-kv ) ];
    const real  fm_f  = f[ IndexLBM::id(ids_lbm_base, 0, 0, 0,   iv, jv, kv ) ];
    const real  fp_ff = f[ IndexLBM::id(ids_lbm_base, iv,jv,kv, -iv,-jv,-kv ) ];

    // delta : 0 ~ 0.5 ~ 1.0 //
    const real  lv0   = lv_obs[id0];
    const real  lv_up = lv_obs[id_up];

    const real  delta   = fabs( lv0 ) / ( fabs(lv0) + fabs(lv_up) + (real)1.0e-6 );
//    const real  delta   = fabs( lv0 ) / ( fabs(lv0) + fabs(lv_up) + ep );

    // force //
//    const real  rb = (real)1.0; // density on boundary //
    constexpr real  rb = (real)1.0; // density on boundary //

    const real  ub = f_int( u_obs[id0], u_obs[id_up], lv0, lv_up );
    const real  vb = f_int( v_obs[id0], v_obs[id_up], lv0, lv_up );
    const real  wb = f_int( w_obs[id0], w_obs[id_up], lv0, lv_up );

    const real  force_bc = force_on_boundary(rb, ub,vb,wb, iv,jv,kv);

//    return  fp_ff;
    return  (delta < (real)0.5) ?
                ((real)2.0*delta*fp_f) + ((real)1.0 - (real)2.0*delta)*fp_ff + force_bc :
                ( fp_f + ((real)2.0*delta - (real)1.0)*fm_f + force_bc ) / ((real)2.0*delta);
}


inline
__HOST__ __DEVICE__
real
factor_feq_D3Q27(const int iv, const int jv, const int kv)
{
    return   (FuncMath::abs(iv,jv,kv) == 0) ?  (real)64.0/(real)216.0
           : (FuncMath::abs(iv,jv,kv) == 1) ?  (real)16.0/(real)216.0
           : (FuncMath::abs(iv,jv,kv) == 2) ?  (real) 4.0/(real)216.0
           : (FuncMath::abs(iv,jv,kv) == 3) ?  (real) 1.0/(real)216.0
           :                                   (real) 0.0;
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
factor_feq_D3Q27()
{
    return   (FuncMath::abs<iv,jv,kv>() == 0) ?  (real)64.0/(real)216.0
           : (FuncMath::abs<iv,jv,kv>() == 1) ?  (real)16.0/(real)216.0
           : (FuncMath::abs<iv,jv,kv>() == 2) ?  (real) 4.0/(real)216.0
           : (FuncMath::abs<iv,jv,kv>() == 3) ?  (real) 1.0/(real)216.0
           :                                     (real) 0.0;
}


inline
__HOST__ __DEVICE__
void
feq_D3Q27(
          real  fs[],
    const real  rho,
    const real  u,
    const real  v,
    const real  w
    )
{
    for (int kv=-1; kv<=1; kv++) {
    for (int jv=-1; jv<=1; jv++) {
    for (int iv=-1; iv<=1; iv++) {
        const int  idv = IndexLBM::idv_lbm(iv,jv,kv);

        fs[idv] = feq_D3Q27(rho,u,v,w, iv,jv,kv);
    }
    }
    }
}


inline
__HOST__ __DEVICE__
real
feq_D3Q27(
    const real  rho,
    const real  u,
    const real  v,
    const real  w,
    const int   iv,
    const int   jv,
    const int   kv
    )
{
    // ref. in: http://www.me.osakafu-u.ac.jp/htlab/html-e/pdf/ascht09.pdf //
    const real  ivel = indexed_vel(u,v,w, iv,jv,kv);
    const real  vel2 = u*u + v*v + w*w;

    // 2nd //
//    return  ( 1.0 + 3.0*ivel + (4.5*ivel*ivel - 1.5*vel2) )  * rho * factor_feq_D3Q27(iv,jv,kv);

    // 3rd //
    return  ( (real)1.0  + (real)3.0*ivel + ((real)4.5*ivel*ivel - (real)1.5*vel2 ) + ((real)4.5*ivel*ivel*ivel - (real)4.5*vel2*ivel) ) * rho * factor_feq_D3Q27(iv,jv,kv);
}


inline
__HOST__ __DEVICE__
real
force_x_D3Q27(
    const real  force,  // mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    )
{
    // ????? sign : + or - //
    return  (real)3.0*(force*iv) * factor_feq_D3Q27 (iv, jv, kv);
}


inline
__HOST__ __DEVICE__
real
force_y_D3Q27(
    const real  force,  // mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    )
{
    // ????? sign : + or - //
    return  (real)3.0*(force*jv) * factor_feq_D3Q27 (iv, jv, kv);
}


inline
__HOST__ __DEVICE__
real
force_z_D3Q27(
    const real  force,  // mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    )
{
    // ????? sign : + or - //
    return  (real)3.0*(force*kv) * factor_feq_D3Q27 (iv, jv, kv);
}


inline
__HOST__ __DEVICE__
real
force_D3Q27(
    const real  force[],  // force_xyz[3] : mass * gravity //
    const int   iv,
    const int   jv,
    const int   kv
    )
{
    // ????? sign : + or - //
    return  (real)3.0*(force[0]*iv + force[1]*jv + force[2]*kv) * factor_feq_D3Q27 (iv, jv, kv);
}


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
    )
{
//    return  ( (real)3.0*( force[0]*(iv - u) + force[1]*(jv - v) + force[2]*(kv - w) ) ) * factor_feq_D3Q27 (iv, jv, kv);

    return  ( (real)3.0*( force[0]*(iv - u) + force[1]*(jv - v) + force[2]*(kv - w) )
            + (real)9.0*( force[0]*iv + force[1]*jv + force[2]*kv )*( iv*u + jv*v + kv*w ) ) * factor_feq_D3Q27 (iv, jv, kv);
}


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
    )
{
    const real  coef = (real)6.0 * rb;
    return  indexed_vel(ub,vb,wb, iv,jv,kv) * factor_feq_D3Q27(iv,jv,kv) * coef;
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
real
force_on_boundary (
    const real  rb,
    const real  ub,
    const real  vb,
    const real  wb
    )
{
    const real  coef = (real)6.0 * rb;
    return  indexed_vel<iv,jv,kv>(ub,vb,wb) * factor_feq_D3Q27<iv,jv,kv>() * coef;
}


inline
__HOST__ __DEVICE__
real
f_int(
    const real f0,
    const real f1,
    const real l0,
    const real l1
    )
{
    constexpr real  ep = (real)1.0e-8;
    return  ( f0*fabs(l1) + f1*fabs(l0) ) / ( fabs(l0) + fabs(l1) + ep );
}


};
