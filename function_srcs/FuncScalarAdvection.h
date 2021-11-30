#pragma once
#ifndef FUNCSCALARADVECTION_H_
#define FUNCSCALARADVECTION_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include "definePrecision.h"
#include "FuncObj.h"
#include "Index.h"
#include "FuncMath.h"


namespace  FuncScalarAdvection {


inline
__HOST__ __DEVICE__
real  dfxv_1st(
    const real  f3[],
    const real  u3[],
    const real  lv3[]
    )
{
//    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
//                         ( FuncObj::is_obj( lv3[1] ) || FuncObj::is_obj( lv3[2] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };
    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
                         ( FuncObj::is_obj( lv3[2] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };

    const real  fs[] = { ( us[0] >= (real)0.0 ) ? f3[0] : f3[1],
                         ( us[1] >= (real)0.0 ) ? f3[1] : f3[2] };

    const real  dft  = - ( us[1]*fs[1] - us[0]*fs[0] );

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
//    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) || FuncObj::is_obj( lv3[1] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
//                         ( FuncObj::is_obj( lv3[1] ) || FuncObj::is_obj( lv3[2] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };
    const real  us[] = { ( FuncObj::is_obj( lv3[0] ) ) ? (real)0.0 : (u3[0] + u3[1])*(real)0.5,
                         ( FuncObj::is_obj( lv3[2] ) ) ? (real)0.0 : (u3[1] + u3[2])*(real)0.5 };

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
real  dfxv_3rd(
    const real  f5[],
    const real  u5[],
    const real  lv5[]
    )
{
    const real  us[] = { ( FuncObj::is_obj( lv5[1] ) ) ? (real)0.0 : (u5[1] + u5[2])*(real)0.5,
                         ( FuncObj::is_obj( lv5[3] ) ) ? (real)0.0 : (u5[2] + u5[3])*(real)0.5 };

    const real  fs1st[] = { ( us[0] >= (real)0.0 ) ? f5[1] : f5[2],
                            ( us[1] >= (real)0.0 ) ? f5[2] : f5[3] };

//    const real  fs2nd[] = { (f5[2] + f5[1]) * (real)0.5,
//                            (f5[3] + f5[2]) * (real)0.5 };

    const real  fs3rd[] = { ( us[0] >= (real)0.0 ) ? ( (real)3.0*f5[2] + (real)6.0*f5[1] - f5[0] ) * (real)0.125 :  ( (real)3.0*f5[1] + (real)6.0*f5[2] - f5[3] ) * (real)0.125,
                            ( us[1] >= (real)0.0 ) ? ( (real)3.0*f5[3] + (real)6.0*f5[2] - f5[1] ) * (real)0.125 :  ( (real)3.0*f5[2] + (real)6.0*f5[3] - f5[4] ) * (real)0.125 };


    const real  fs_tmp[] = { ( FuncObj::is_obj( lv5[0] ) || FuncObj::is_obj( lv5[1] ) || FuncObj::is_obj( lv5[2] ) || FuncObj::is_obj( lv5[3] ) ) ? fs1st[0] : fs3rd[0] ,
                             ( FuncObj::is_obj( lv5[1] ) || FuncObj::is_obj( lv5[2] ) || FuncObj::is_obj( lv5[3] ) || FuncObj::is_obj( lv5[4] ) ) ? fs1st[1] : fs3rd[1] };

    const real  ww[] = { 0.95, 0.95 };
//    const real  ep   = (real)1.0e-8;
//    const real  ww[] = { ( (real)1.0 - fabs( (f5[0] - f5[1]) - (f5[2] - f5[3]) ) / ( fabs(f5[0] - f5[1]) + fabs(f5[2] - f5[3]) + ep ) ) * (real)1.0,
//                         ( (real)1.0 - fabs( (f5[1] - f5[2]) - (f5[3] - f5[4]) ) / ( fabs(f5[1] - f5[2]) + fabs(f5[3] - f5[4]) + ep ) ) * (real)1.0 };

    const real  fs[] = { ww[0]*fs_tmp[0] + ((real)1.0-ww[0])*fs1st[0],
                         ww[1]*fs_tmp[1] + ((real)1.0-ww[1])*fs1st[1] };


    const real  fx[] = { f5[2] - f5[1],
                         f5[3] - f5[2] };

    const real  dft  = - ( us[1]*fs[1] - us[0]*fs[0] );
    const real  dftt =   ( us[1]*us[1]*fx[1] - us[0]*us[0]*fx[0] );


    return  dft + (real)0.5*dftt;
}


inline
__HOST__ __DEVICE__
real f_average(
    const int    ids[],
    const real* scalar,
    const real* lv_obj
    )
{
    int   id_tmp;
    real  f = (real)0.0;
    int   count = 0;

    id_tmp = Index::id(ids, 0, 0, 0);
    f += scalar[ id_tmp ];   count++;

    id_tmp = Index::id(ids, -1, 0, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  1, 0, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, -1, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0,  1, 0);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, 0, -1);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(ids,  0, 0,  1);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    return  ( f/(real)count );
}


inline
__HOST__ __DEVICE__
real f_average(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[],
    const real*  scalar,
    const real*  lv_obj
    )
{
    int   id_tmp;
    real  f = (real)0.0;
    int   count = 0;

    id_tmp = Index::id(i, j, k, offset3d);
    f += scalar[ id_tmp ];   count++;

    id_tmp = Index::id(i-1, j, k, offset3d);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(i+1, j, k, offset3d);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j-1, k, offset3d);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j+1, k, offset3d);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j, k-1, offset3d);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    id_tmp = Index::id(i, j, k+1, offset3d);
    if ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) { f += scalar[ id_tmp ]; count++; }

    return  ( f/(real)count );
}


inline
__HOST__ __DEVICE__
real superbee(const real r)
{
    return FuncMath::max3( (real)0.0, FuncMath::min((real)2.0*r, (real)1.0), FuncMath::min(r, (real)2.0) );
}


inline
__HOST__ __DEVICE__
real  dfxx(
    const real  f3[],
    const real  lv3[],
    const real  coef_vis
    )
{
    const real  _f3[] = { (FuncObj::is_obj( lv3[0] ) ) ? f3[1] : f3[0],
                          (FuncObj::is_obj( lv3[1] ) ) ? f3[1] : f3[1],
                          (FuncObj::is_obj( lv3[2] ) ) ? f3[1] : f3[2]  };

    const real  _fx[] = { _f3[1] - _f3[0],
                          _f3[2] - _f3[1] };
//    const real  fx[] = { FuncObj::is_obj( lv3[0] ) ? _fx[0]*(real)2.0 : _fx[0],
//                         FuncObj::is_obj( lv3[2] ) ? _fx[1]*(real)2.0 : _fx[1]  };
    const real  fx[] = { _fx[0], _fx[1]  };


    const real  dft = coef_vis*(fx[1] - fx[0]);

    return  dft;
}


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
    const real*  lv_obj,
    const real   dx,
    const real   dt
    )
{
    const int  id0 = Index::id0(ids);
    if ( FuncObj::is_obj( lv_obj[id0] ) ) { scalar_n[id0] = f_average( ids, scalar, lv_obj ); return; }

    const real  f3x[] = { scalar[Index::id(ids, -1,  0,  0)], scalar[Index::id(ids,  0,  0,  0)], scalar[Index::id(ids,  1,  0,  0)] };
    const real  f3y[] = { scalar[Index::id(ids,  0, -1,  0)], scalar[Index::id(ids,  0,  0,  0)], scalar[Index::id(ids,  0,  1,  0)] };
    const real  f3z[] = { scalar[Index::id(ids,  0,  0, -1)], scalar[Index::id(ids,  0,  0,  0)], scalar[Index::id(ids,  0,  0,  1)] };

    const real  u3[] = { u[Index::id(ids, -1,  0,  0)], u[Index::id(ids,  0,  0,  0)], u[Index::id(ids,  1,  0,  0)] };
    const real  v3[] = { v[Index::id(ids,  0, -1,  0)], v[Index::id(ids,  0,  0,  0)], v[Index::id(ids,  0,  1,  0)] };
    const real  w3[] = { w[Index::id(ids,  0,  0, -1)], w[Index::id(ids,  0,  0,  0)], w[Index::id(ids,  0,  0,  1)] };

    const real  lv3x[] = { lv_obj[Index::id(ids, -1,  0,  0)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  1,  0,  0)] };
    const real  lv3y[] = { lv_obj[Index::id(ids,  0, -1,  0)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  0,  1,  0)] };
    const real  lv3z[] = { lv_obj[Index::id(ids,  0,  0, -1)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  0,  0,  1)] };

    const real  dscalar1 = (   dfxv_1st(f3x, u3, lv3x)
                             + dfxv_1st(f3y, v3, lv3y)
                             + dfxv_1st(f3z, w3, lv3z) );

    const real  dscalar2 = (   dfxv_2nd(f3x, u3, lv3x)
                             + dfxv_2nd(f3y, v3, lv3y)
                             + dfxv_2nd(f3z, w3, lv3z) );

    const real  ww      = (real)0.0;
//    const real  ww      = (real)0.5;
    const real  dAdv = ((real)1.0-ww)*dscalar1 + ww*dscalar2;


    const real coef_vis_sgs   = (sgs_vis[id0]/10.0);
    const real _coef_vis_total = (coef_vis_sgs) / (dx*dx); // because dx of the fxx set to 1 //
    const real coef_vis_total = ( _coef_vis_total < (real)1.0/(real)6.0/dt ) ? _coef_vis_total : (real)1.0/(real)6.0/dt;

    const real  dscalarxx  = (   dfxx(f3x, lv3x, coef_vis_total)
                               + dfxx(f3y, lv3y, coef_vis_total)
                               + dfxx(f3z, lv3z, coef_vis_total) );

    const real  dDif = 0.0;
//    const real  dDif = dscalarxx * dt;

    scalar_n[id0] =  scalar[id0] + dAdv + dDif;
}


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
    const real*  lv_obj,
    const real*  sc_scalar, const bool flag_source,
    const real   dx,
    const real   dt
    )
{
    const int  id0 = Index::id0(ids);
    if ( FuncObj::is_obj( lv_obj[id0] ) ) { scalar_n[id0] = f_average( ids, scalar, lv_obj ); return; }

    const real  f3x[] = { scalar[Index::id(ids, -1,  0,  0)], scalar[Index::id(ids,  0,  0,  0)], scalar[Index::id(ids,  1,  0,  0)] };
    const real  f3y[] = { scalar[Index::id(ids,  0, -1,  0)], scalar[Index::id(ids,  0,  0,  0)], scalar[Index::id(ids,  0,  1,  0)] };
    const real  f3z[] = { scalar[Index::id(ids,  0,  0, -1)], scalar[Index::id(ids,  0,  0,  0)], scalar[Index::id(ids,  0,  0,  1)] };

    const real  u3[] = { u[Index::id(ids, -1,  0,  0)], u[Index::id(ids,  0,  0,  0)], u[Index::id(ids,  1,  0,  0)] };
    const real  v3[] = { v[Index::id(ids,  0, -1,  0)], v[Index::id(ids,  0,  0,  0)], v[Index::id(ids,  0,  1,  0)] };
    const real  w3[] = { w[Index::id(ids,  0,  0, -1)], w[Index::id(ids,  0,  0,  0)], w[Index::id(ids,  0,  0,  1)] };

    const real  lv3x[] = { lv_obj[Index::id(ids, -1,  0,  0)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  1,  0,  0)] };
    const real  lv3y[] = { lv_obj[Index::id(ids,  0, -1,  0)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  0,  1,  0)] };
    const real  lv3z[] = { lv_obj[Index::id(ids,  0,  0, -1)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  0,  0,  1)] };

    const real  dscalar1 = (   dfxv_1st(f3x, u3, lv3x)
                             + dfxv_1st(f3y, v3, lv3y)
                             + dfxv_1st(f3z, w3, lv3z) );

    const real  dscalar2 = (   dfxv_2nd(f3x, u3, lv3x)
                             + dfxv_2nd(f3y, v3, lv3y)
                             + dfxv_2nd(f3z, w3, lv3z) );

    const real  ww   = (real)0.0;
//    const real  ww   = (real)0.5;
    const real  dAdv = ((real)1.0-ww)*dscalar1 + ww*dscalar2;

    const real  source  = flag_source ? (sc_scalar[id0] * dt) : 0;


///    const real coef_vis_sgs   = (sgs_vis[id0]/10.0);
///    const real _coef_vis_total = (coef_vis_sgs) / (dx*dx); // because dx of the fxx set to 1 //
///    const real coef_vis_total = ( _coef_vis_total < (real)1.0/(real)6.0/dt ) ? _coef_vis_total : (real)1.0/(real)6.0/dt;
///
///    const real  dscalarxx  = (   dfxx(f3x, lv3x, coef_vis_total)
///                               + dfxx(f3y, lv3y, coef_vis_total)
///                               + dfxx(f3z, lv3z, coef_vis_total) );
//    const real  dDif = dscalarxx * dt;

    const real  dDif = 0.0;

    scalar_n[id0] =  scalar[id0] + dAdv + dDif + source;
}


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
    const real*  lv_obj,
    const real*  sc_scalar, const bool flag_source,
    const real   dt
    )
{
    const int  id0 = Index::id0(i, j, k, offset3d);
//    if ( FuncObj::is_obj( lv_obj[id0] ) ) { scalar_n[id0] = 0.0; return; }
    if ( FuncObj::is_obj( lv_obj[id0] ) ) { scalar_n[id0] = f_average(i,j,k,offset3d, scalar, lv_obj); return; }

    const real  f5x[] = { scalar[Index::id(i-2, j  , k  , offset3d)], scalar[Index::id(i-1, j  , k  , offset3d)], scalar[Index::id( i, j, k, offset3d)], scalar[Index::id(i+1, j  , k  , offset3d)], scalar[Index::id(i+2, j  , k  , offset3d)] };
    const real  f5y[] = { scalar[Index::id(i  , j-2, k  , offset3d)], scalar[Index::id(i  , j-1, k  , offset3d)], scalar[Index::id( i, j, k, offset3d)], scalar[Index::id(i  , j+1, k  , offset3d)], scalar[Index::id(i  , j+2, k  , offset3d)] };
    const real  f5z[] = { scalar[Index::id(i  , j  , k-2, offset3d)], scalar[Index::id(i  , j  , k-1, offset3d)], scalar[Index::id( i, j, k, offset3d)], scalar[Index::id(i  , j  , k+1, offset3d)], scalar[Index::id(i  , j  , k+2, offset3d)] };

    const real  u5[] = { u[Index::id(i-2, j  , k  , offset3d)], u[Index::id(i-1, j  , k  , offset3d)], u[Index::id( i, j, k, offset3d)], u[Index::id(i+1, j  , k  , offset3d)], u[Index::id(i+2, j  , k  , offset3d)] };
    const real  v5[] = { v[Index::id(i  , j-2, k  , offset3d)], v[Index::id(i  , j-1, k  , offset3d)], v[Index::id( i, j, k, offset3d)], v[Index::id(i  , j+1, k  , offset3d)], v[Index::id(i  , j+2, k  , offset3d)] };
    const real  w5[] = { w[Index::id(i  , j  , k-2, offset3d)], w[Index::id(i  , j  , k-1, offset3d)], w[Index::id( i, j, k, offset3d)], w[Index::id(i  , j  , k+1, offset3d)], w[Index::id(i  , j  , k+2, offset3d)] };

    const real  lv5x[] = { lv_obj[Index::id(i-2, j  , k  , offset3d)], lv_obj[Index::id(i-1, j  , k  , offset3d)], lv_obj[Index::id( i, j, k, offset3d)], lv_obj[Index::id(i+1, j  , k  , offset3d)] , lv_obj[Index::id(i+2, j  , k  , offset3d)] };
    const real  lv5y[] = { lv_obj[Index::id(i  , j-2, k  , offset3d)], lv_obj[Index::id(i  , j-1, k  , offset3d)], lv_obj[Index::id( i, j, k, offset3d)], lv_obj[Index::id(i  , j+1, k  , offset3d)] , lv_obj[Index::id(i  , j+2, k  , offset3d)] };
    const real  lv5z[] = { lv_obj[Index::id(i  , j  , k-2, offset3d)], lv_obj[Index::id(i  , j  , k-1, offset3d)], lv_obj[Index::id( i, j, k, offset3d)], lv_obj[Index::id(i  , j  , k+1, offset3d)] , lv_obj[Index::id(i  , j  , k+2, offset3d)] };

    const real  dscalar3 = (   dfxv_3rd(f5x, u5, lv5x)
                             + dfxv_3rd(f5y, v5, lv5y)
                             + dfxv_3rd(f5z, w5, lv5z) );

    const real  dscalar = dscalar3;
//    const real  dscalar = (real)0.0;

    const real  source  = flag_source ? (sc_scalar[id0] * dt) : 0;

    scalar_n[id0] =  scalar[id0] + dscalar + source;
}


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
    const real*  lv_obj,
    const real   dx,
    const real   dt
    )
{
    scalar_advection_2nd(
        ids,
        scalar, scalar_n,
        u,v,w,
        sgs_vis,
        lv_obj,
        dx, dt
        );
}


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
    const real*  lv_obj,
    const real*  sc_scalar, const bool flag_source,
    const real   dx,
    const real   dt
    )
{
    scalar_advection_2nd(
        ids,
        scalar, scalar_n,
        u,v,w,
        sgs_vis,
        lv_obj,
        sc_scalar, flag_source,
        dx, dt
        );
}


inline
__HOST__ __DEVICE__
void  scalar_advection_high_order(
    const int    i,
    const int    j,
    const int    k,
    const int    offset3d[], // not ids //
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obj,
    const real*  sc_scalar, const bool flag_source,
    const real   dt
    )
{
    scalar_advection_3rd(
        i, j, k,
        offset3d,
        scalar, scalar_n,
        u,v,w,
        lv_obj,
        sc_scalar, flag_source,
        dt
        );
}


} // namespace FuncScalarAdvection


#endif
