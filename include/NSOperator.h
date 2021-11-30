#pragma once
#ifndef NSOPERATOR_H_
#define NSOPERATOR_H_


#include "IndexMG.h"
#include "FuncObj.h"


namespace  NSOperator {


template <int NX_LEAF, typename T>
inline __HD__
T average_neighbor(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T val_av  = 0.0;
    val_av   += val[ID( 0, 0, 0)];
    val_av   += val[ID(-1, 0, 0)];
    val_av   += val[ID( 1, 0, 0)];
    val_av   += val[ID( 0,-1, 0)];
    val_av   += val[ID( 0, 1, 0)];
    val_av   += val[ID( 0, 0,-1)];
    val_av   += val[ID( 0, 0, 1)];

    constexpr T inv7 = 1.0/7.0;

    return  (val_av * inv7);
}


template <int NX_LEAF, typename T>
inline __HD__
T average_neighbor_w_obj(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val,
    const T*  lv_obj
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    int count = 0;
    
    T val_av  = 0.0; count++;
    if ( FuncObj::is_fluid( lv_obj[ID( 0, 0, 0)] ) ) { val_av   += val[ID( 0, 0, 0)];  count++; }
    if ( FuncObj::is_fluid( lv_obj[ID(-1, 0, 0)] ) ) { val_av   += val[ID(-1, 0, 0)];  count++; }
    if ( FuncObj::is_fluid( lv_obj[ID( 1, 0, 0)] ) ) { val_av   += val[ID( 1, 0, 0)];  count++; }
    if ( FuncObj::is_fluid( lv_obj[ID( 0,-1, 0)] ) ) { val_av   += val[ID( 0,-1, 0)];  count++; }
    if ( FuncObj::is_fluid( lv_obj[ID( 0, 1, 0)] ) ) { val_av   += val[ID( 0, 1, 0)];  count++; }
    if ( FuncObj::is_fluid( lv_obj[ID( 0, 0,-1)] ) ) { val_av   += val[ID( 0, 0,-1)];  count++; }
    if ( FuncObj::is_fluid( lv_obj[ID( 0, 0, 1)] ) ) { val_av   += val[ID( 0, 0, 1)];  count++; }

    const T inv = 1.0/(T)count;

    return  (val_av * inv);
}


template <int NX_LEAF, typename T>
inline __HD__
T average(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T val_av    = 0.0;
    T weight_av = 0.0;
    for (int kk=-1; kk<=1; kk++) {
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        val_av    += val[ID(ii,jj,kk)];
    }
    }
    }
    constexpr T inv27 = 1.0/27.0;

    return  (val_av * inv27);
}


template <int NX_LEAF, typename T>
inline __HD__
T average_w_obj(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val,
    const T*  lv_obj
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T val_av    = 0.0;
    T weight_av = 0.0;
    for (int kk=-1; kk<=1; kk++) {
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        const T _weight = FuncObj::is_obj( lv_obj[ID(ii,jj,kk)] )  ? 1.0 : 1.0e5;
    
        val_av    += _weight * val[ID(ii,jj,kk)];
        weight_av += _weight;
    }
    }
    }

    return  (val_av / weight_av);
}


template <int NX_LEAF, typename T>
inline __HD__
void offsets_w_obj(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  lv_obj,
          T&  ix_off,
          T&  iy_off,
          T&  iz_off
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T weight_av = 0.0;

    ix_off = 0.0;
    iy_off = 0.0;
    iz_off = 0.0;
    for (int kk=-1; kk<=1; kk++) {
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        const T _weight = FuncObj::is_obj( lv_obj[ID(ii,jj,kk)] )  ? 1.0 : 1.0e5;
    
        weight_av += _weight;

        ix_off += _weight * ii;
        iy_off += _weight * jj;
        iz_off += _weight * kk;
    }
    }
    }
    ix_off /= weight_av;
    iy_off /= weight_av;
    iz_off /= weight_av;
}


template <int NX_LEAF, typename T>
inline __HD__
T xy_average_w_obj(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val,
    const T*  lv_obj
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T val_av    = 0.0;
    T weight_av = 0.0;

    const int kk = 0;
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        const T _weight = FuncObj::is_obj( lv_obj[ID(ii,jj,kk)] )  ? 1.0 : 1.0e5;
    
        val_av    += _weight * val[ID(ii,jj,kk)];
        weight_av += _weight;
    }
    }

    return  (val_av / weight_av);
}


template <int NX_LEAF, typename T>
inline __HD__
T average_w_noslip(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val,
    const T*  lv_obj
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T val_av    = 0.0;
    for (int kk=-1; kk<=1; kk++) {
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        const T _val = FuncObj::is_obj( lv_obj[ID(ii,jj,kk)] )  ? 0.0 : val[ID(ii,jj,kk)];
    
        val_av    += _val;
    }
    }
    }

    return  (val_av / 27.0);
}


template <int NX_LEAF, typename T>
inline __HD__
T xy_average_w_noslip(
    const int offset3d[], // int  offset3d[27];
    const int i,
    const int j,
    const int k,
    const T*  val,
    const T*  lv_obj
    )
{
    auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };
    
    T val_av    = 0.0;

    const int kk = 0;
    for (int jj=-1; jj<=1; jj++) {
    for (int ii=-1; ii<=1; ii++) {
        const T _val = FuncObj::is_obj( lv_obj[ID(ii,jj,kk)] )  ? 0.0 : val[ID(ii,jj,kk)];
    
        val_av    += _val;
    }
    }

    return  (val_av / 9.0);
}


} // NSOperator

#endif
