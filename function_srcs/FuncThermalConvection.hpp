#include "FuncThermalConvection.h"
#include "FuncObj.h"
#include "Index.h"
#include "FuncMath.h"

#include "defineFluidProperty.h"


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
    )
{
    cal_thermal_convection_2nd(
        ids,
        T, T_n,
        u,v,w,
        sgs_vis,
        lv_obs,
        T_obs, Tn_obs, bcTypes_T, time_weight0,
        coef_heatf,
        dx, dt
        );
}


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
    )
{
    const int  id0 = Index::id0(ids);
    if ( FuncObj::is_obj( lv_obs[id0] ) ) {
        if      ( bcTypes_T[id0] == BCTypes::BCDirichlet ) { T_n[id0] = T_obs[id0];          }
//        else if ( bcTypes_T[id0] == BCTypes::BCRegion    ) { T_n[id0] = T[id0];          }
        else if ( bcTypes_T[id0] == BCTypes::BCRegion    ) { T_n[id0] = f_average( ids, T );          }
        else                                               { T_n[id0] = f_average_w_obs( ids, T, lv_obs ); }
        return;
    }

    const real  f3x[] = { T[Index::id(ids, -1,  0,  0)], T[Index::id(ids,  0,  0,  0)], T[Index::id(ids,  1,  0,  0)] };
    const real  f3y[] = { T[Index::id(ids,  0, -1,  0)], T[Index::id(ids,  0,  0,  0)], T[Index::id(ids,  0,  1,  0)] };
    const real  f3z[] = { T[Index::id(ids,  0,  0, -1)], T[Index::id(ids,  0,  0,  0)], T[Index::id(ids,  0,  0,  1)] };

    const real  u3[] = { u[Index::id(ids, -1,  0,  0)], u[Index::id(ids,  0,  0,  0)], u[Index::id(ids,  1,  0,  0)] };
    const real  v3[] = { v[Index::id(ids,  0, -1,  0)], v[Index::id(ids,  0,  0,  0)], v[Index::id(ids,  0,  1,  0)] };
    const real  w3[] = { w[Index::id(ids,  0,  0, -1)], w[Index::id(ids,  0,  0,  0)], w[Index::id(ids,  0,  0,  1)] };

    const real  lv3x[] = { lv_obs[Index::id(ids, -1,  0,  0)], lv_obs[Index::id(ids,  0,  0,  0)], lv_obs[Index::id(ids,  1,  0,  0)] };
    const real  lv3y[] = { lv_obs[Index::id(ids,  0, -1,  0)], lv_obs[Index::id(ids,  0,  0,  0)], lv_obs[Index::id(ids,  0,  1,  0)] };
    const real  lv3z[] = { lv_obs[Index::id(ids,  0,  0, -1)], lv_obs[Index::id(ids,  0,  0,  0)], lv_obs[Index::id(ids,  0,  0,  1)] };

    const real  Tobs3x[] = { f_time_weight(T_obs[Index::id(ids, -1,  0,  0)], Tn_obs[Index::id(ids, -1,  0,  0)], time_weight0), f_time_weight(T_obs[Index::id(ids,  0,  0,  0)], Tn_obs[Index::id(ids,  0,  0,  0)], time_weight0), f_time_weight(T_obs[Index::id(ids,  1,  0,  0)], Tn_obs[Index::id(ids,  1,  0,  0)], time_weight0) };
    const real  Tobs3y[] = { f_time_weight(T_obs[Index::id(ids,  0, -1,  0)], Tn_obs[Index::id(ids,  0, -1,  0)], time_weight0), f_time_weight(T_obs[Index::id(ids,  0,  0,  0)], Tn_obs[Index::id(ids,  0,  0,  0)], time_weight0), f_time_weight(T_obs[Index::id(ids,  0,  1,  0)], Tn_obs[Index::id(ids,  0,  1,  0)], time_weight0) };
    const real  Tobs3z[] = { f_time_weight(T_obs[Index::id(ids,  0,  0, -1)], Tn_obs[Index::id(ids,  0,  0, -1)], time_weight0), f_time_weight(T_obs[Index::id(ids,  0,  0,  0)], Tn_obs[Index::id(ids,  0,  0,  0)], time_weight0), f_time_weight(T_obs[Index::id(ids,  0,  0,  1)], Tn_obs[Index::id(ids,  0,  0,  1)], time_weight0) };

    const int   bcTypesT3x[] = { bcTypes_T[Index::id(ids, -1,  0,  0)], bcTypes_T[Index::id(ids,  0,  0,  0)], bcTypes_T[Index::id(ids,  1,  0,  0)] };
    const int   bcTypesT3y[] = { bcTypes_T[Index::id(ids,  0, -1,  0)], bcTypes_T[Index::id(ids,  0,  0,  0)], bcTypes_T[Index::id(ids,  0,  1,  0)] };
    const int   bcTypesT3z[] = { bcTypes_T[Index::id(ids,  0,  0, -1)], bcTypes_T[Index::id(ids,  0,  0,  0)], bcTypes_T[Index::id(ids,  0,  0,  1)] };

//    // 2d //
//    const real  Tx1st = (   dfxv_1st(f3x, u3, lv3x)
//                          + dfxv_1st(f3y, v3, lv3y) );
//
//    const real  Tx2nd = (   dfxv_2nd(f3x, u3, lv3x)
//                          + dfxv_2nd(f3y, v3, lv3y) );
//
//    const real  dTxx  = (   dfxx(f3x, lv3x, Tobs3x, bcTypesT3x, coef_heatf)
//                          + dfxx(f3y, lv3y, Tobs3y, bcTypesT3y, coef_heatf)
//                          + dfxx(f3z, lv3z, Tobs3z, bcTypesT3z, coef_heatf) );

    // 3d //
#if 0
    const real  Tx1st = (   dfxv_1st(f3x, u3, lv3x)
                          + dfxv_1st(f3y, v3, lv3y)
                          + dfxv_1st(f3z, w3, lv3z) );

    const real  Tx2nd = (   dfxv_2nd(f3x, u3, lv3x)
                          + dfxv_2nd(f3y, v3, lv3y)
                          + dfxv_2nd(f3z, w3, lv3z) );
#else
    const real  Tx1st = (   dfxv_adv_1st(f3x, u3, lv3x)
                          + dfxv_adv_1st(f3y, v3, lv3y)
                          + dfxv_adv_1st(f3z, w3, lv3z) );

    const real  Tx2nd = (   dfxv_adv_2nd(f3x, u3, lv3x)
                          + dfxv_adv_2nd(f3y, v3, lv3y)
                          + dfxv_adv_2nd(f3z, w3, lv3z) );
#endif

//    const real coef_heatf_total =  coef_heatf / (dx*dx); // because dx of the fxx set to 1 //
    const real coef_heatf_sgs   = (sgs_vis[id0]/air_property::Pr);
    const real _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
    const real coef_heatf_total = ( _coef_heatf_total < (real)1.0/(real)6.0/dt ) ? _coef_heatf_total : (real)1.0/(real)6.0/dt;

    const real  dTxx  = (   dfxx(f3x, lv3x, Tobs3x, bcTypesT3x, coef_heatf_total)
                          + dfxx(f3y, lv3y, Tobs3y, bcTypesT3y, coef_heatf_total)
                          + dfxx(f3z, lv3z, Tobs3z, bcTypesT3z, coef_heatf_total) );


    const real  ww   = (real)0.0;
//    const real  ww   = (real)0.05;
//    const real  ww   = (real)0.20;
//    const real  ww   = (real)0.50;
//    const real  ww   = (real)0.90;
    const real  dAdv = ((real)1.0-ww)*Tx1st + ww*Tx2nd;
    const real  dDif = dTxx * dt;

    T_n  [id0] =  T[id0] + dAdv + dDif;
//    T_n  [id0] =  T[id0]        + dDif;
//    T_n  [id0] =  T[id0];
}


inline
__HOST__ __DEVICE__
real  dfxv_1st(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    )
{
    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
                         ( FuncObj::is_obj( lv3[2] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };

    const real  fs[] = { ( us[0] >= (real)0.0 ) ? f3[0] : f3[1],
                         ( us[1] >= (real)0.0 ) ? f3[1] : f3[2] };

    const real  dft  = - ( us[1]*fs[1] - us[0]*fs[0] );

    return  dft;
}


inline
__HOST__ __DEVICE__
real  dfxv_adv_1st(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    )
{
    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
                         ( FuncObj::is_obj( lv3[2] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };

    const real  fxs[] = { (f3[1] - f3[0]),
                          (f3[2] - f3[1]) };

#if 0
    const real  ufx_central = ( us[0]*fxs[0] + us[1]*fxs[1] ) * (real)0.5;
    const real  fxx = ( f3[2] - (real)2.0*f3[1] + f3[0] );
//    const real  ufx_upwind = ufx_central - fabs( us[0] + us[1] ) * (real)0.5 * fxx * (real)0.5;
    const real  ufx_upwind = ufx_central - ( fabs(us[0]) + fabs(us[1]) ) * (real)0.5 * fxx * (real)0.5;
//    const real  ufx_upwind = ufx_central - fabs( u3[1] ) * fxx * (real)0.5;
#else // unstable //
//    const real  ufx_upwind = (u3[1] > (real)0.0) ? us[0]*fxs[0] : us[1]*fxs[1];
    const real  ufx_upwind = (u3[1] > (real)0.0) ? u3[1]*fxs[0] : u3[1]*fxs[1];
#endif

    const real  dft  = - ufx_upwind;

    return  dft;
}


inline
__HOST__ __DEVICE__
real  dfxv_2nd(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    )
{
    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
                         ( FuncObj::is_obj( lv3[2] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };

    const real  fs[] = { (f3[1] + f3[0]) * (real)0.5,
                         (f3[2] + f3[1]) * (real)0.5 };

    const real  fx[] = { f3[1] - f3[0],
                         f3[2] - f3[1] };

    const real  dft  = - ( us[1]*fs[1] - us[0]*fs[0] );
    const real  dftt =   ( us[1]*us[1]*fx[1] - us[0]*us[0]*fx[0] );


    return  dft + (real)0.5*dftt;
}


inline
__HOST__ __DEVICE__
real  dfxv_adv_2nd(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    )
{
    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
                         ( FuncObj::is_obj( lv3[2] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };

    const real  fxs[] = { f3[1] - f3[0],
                          f3[2] - f3[1] };

    const real  ufx_central = ( us[0]*fxs[0] + us[1]*fxs[1] ) * (real)0.5;

    const real  dft  = - ( ufx_central );
//    const real  dftt =   ( us[1]*us[1]*fxs[1] - us[0]*us[0]*fxs[0] );
//    const real  dftt =  u3[1] * ( us[1]*fxs[1] - us[0]*fxs[0] );
//    const real  dftt =  (us[0] + us[1])*(real)0.5 * ( us[1]*fxs[1] - us[0]*fxs[0] );
//    const real  dftt =  (u3[1]*u3[1])*(real)0.5 * ( fxs[1] - fxs[0] );
    const real  dftt =  (u3[1]*u3[1])*(real)0.5 * ( fxs[1] - fxs[0] );


    return  dft + (real)0.5*dftt;
}


inline
__HOST__ __DEVICE__
real  dfxx(
    const real  f3[],
    const real  lv3[],
    const real  Tobs3[],
    const int   bcTypesT3[],
    const real  coef_heatf
    )
{
    const real  _f3[] = { (FuncObj::is_obj( lv3[0] ) && bcTypesT3[0] == BCTypes::BCDirichlet) ? Tobs3[0] : f3[0],
                          (FuncObj::is_obj( lv3[1] ) && bcTypesT3[1] == BCTypes::BCDirichlet) ? Tobs3[1] : f3[1],
                          (FuncObj::is_obj( lv3[2] ) && bcTypesT3[2] == BCTypes::BCDirichlet) ? Tobs3[2] : f3[2]  };

    const real  __f3[] = { (FuncObj::is_obj( lv3[0] ) && bcTypesT3[0] == BCTypes::BCNeumann) ? _f3[1] : _f3[0],
                           (FuncObj::is_obj( lv3[1] ) && bcTypesT3[1] == BCTypes::BCNeumann) ? _f3[1] : _f3[1],
                           (FuncObj::is_obj( lv3[2] ) && bcTypesT3[2] == BCTypes::BCNeumann) ? _f3[1] : _f3[2]  };

//    const real  fx[] = { __f3[1] - __f3[0],
//                         __f3[2] - __f3[1] };

    // dirichlet : wall is placed at the cell surface (dx -> 0.5*dx)
    const real  _fx[] = { __f3[1] - __f3[0],
                          __f3[2] - __f3[1] };
    const real  fx[] = { FuncObj::is_obj( lv3[0] ) ? _fx[0]*(real)2.0 : _fx[0],
                         FuncObj::is_obj( lv3[2] ) ? _fx[1]*(real)2.0 : _fx[1]  };

//    // cell center ?? //
//    const real  fx[] = { _fx[0], _fx[1]  };


    const real  dft = coef_heatf*(fx[1] - fx[0]);

    return  dft;
}


inline
__HOST__ __DEVICE__
real  rho_buoyancy(const real  T)
{
//    return  (real)1.0;
    return  (real)1.0 - air_property::beta*(T - fluid_property::Temperature_bouyancy);
}


inline
__HOST__ __DEVICE__
real  T_buoyancy(const real  T)
{
    return  (T - fluid_property::Temperature_bouyancy);
}


inline
__HOST__ __DEVICE__
real  T_buoyancy(
    const real T,
    const real Txm, const real Txp,
    const real Tym, const real Typ,
    const real Tzm, const real Tzp
    )
{
    const real Taverage = ( (real)6.0*T + Txm + Txp + Tym + Typ + Tzm + Tzp ) / (real)12.0; // diverge //
    return  (Taverage - fluid_property::Temperature_bouyancy);
}


inline
__HOST__ __DEVICE__
real f_average(
    const int    ids[],
    const real* val
    )
{
    int   id_tmp;
    real  f = (real)0.0;
    int   count = 0;

    id_tmp = Index::id(ids, 0, 0, 0);
    f += val[ id_tmp ];   count++;

    id_tmp = Index::id(ids, -1, 0, 0);
    f += val[ id_tmp ]; count++;

    id_tmp = Index::id(ids,  1, 0, 0);
    f += val[ id_tmp ]; count++;

    id_tmp = Index::id(ids,  0, -1, 0);
    f += val[ id_tmp ]; count++;

    id_tmp = Index::id(ids,  0,  1, 0);
    f += val[ id_tmp ]; count++;

    id_tmp = Index::id(ids,  0, 0, -1);
    f += val[ id_tmp ]; count++;

    id_tmp = Index::id(ids,  0, 0,  1);
    f += val[ id_tmp ]; count++;

    return  ( f/(real)count );
}


inline
__HOST__ __DEVICE__
real f_average_w_obs(
    const int    ids[],
    const real* val,
    const real* lv_obs
    )
{
    int   id_tmp;
    real  f      = (real)0.0;
    real  weight = (real)0.0;
    real  tmp_weight;

    constexpr real cF = (real)1000.0;
    constexpr real cO = (real)1.0;

    for (int kk=-1; kk<=1; kk++) {
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        id_tmp = Index::id(ids, ii, jj, kk);
        tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
        f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;
    }
    }
    }

/*
    id_tmp = Index::id(ids, 0, 0, 0);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;


    id_tmp = Index::id(ids, -1, 0, 0);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;


    id_tmp = Index::id(ids,  1, 0, 0);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;


    id_tmp = Index::id(ids,  0, -1, 0);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;


    id_tmp = Index::id(ids,  0,  1, 0);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;


    id_tmp = Index::id(ids,  0, 0, -1);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;


    id_tmp = Index::id(ids,  0, 0,  1);
    tmp_weight = ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) ? cF : cO;
    f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;
*/

    return  ( f/weight );
}


inline
__HOST__ __DEVICE__
real f_average(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[],
    const real*  val,
    const real*  lv_obs
    )
{
    int   id_tmp;
    real  f = (real)0.0;
    int   count = 0;

    id_tmp = Index::id(i, j, k, offset3d);
    f += val[ id_tmp ];   count++;

    id_tmp = Index::id(i-1, j, k, offset3d);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(i+1, j, k, offset3d);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j-1, k, offset3d);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j+1, k, offset3d);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j, k-1, offset3d);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j, k+1, offset3d);
    if ( FuncObj::is_fluid( lv_obs[ id_tmp ] ) ) { f += val[ id_tmp ]; count++; }

    return  ( f/(real)count );
}


inline
__HOST__ __DEVICE__
real f_time_weight(
    const real f0,
    const real fn,
    const real weight0
    )
{
    return weight0*f0 + ((real)1.0-weight0)*fn;
}

};
