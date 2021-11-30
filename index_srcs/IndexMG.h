#pragma once
#ifndef INDEXMG_H_
#define INDEXMG_H_


#include <iostream>
#include "defineCal.h"


namespace  IndexMG {


inline
__HOST__ __DEVICE__
int  idv(const int iv, const int jv, const int kv) noexcept
{
    return  ( (iv+1) + 3*(jv+1) + 9*(kv+1) );
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
constexpr int  idv() noexcept
{
    return  ( (iv+1) + 3*(jv+1) + 9*(kv+1) );
}


template<int NX_LEAF>
__HOST__ __DEVICE__
int  id(
    const int   ix,
    const int   iy,
    const int   iz
    ) noexcept
{
    return  ix + NX_LEAF*iy + NX_LEAF*NX_LEAF*iz;
}


template<int NX_LEAF, int ix, int iy, int iz>
__HOST__ __DEVICE__
constexpr int  id() noexcept
{
    return  ix + NX_LEAF*iy + NX_LEAF*NX_LEAF*iz;
}


inline
__HOST__ __DEVICE__
int  id(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx,
    const int   ny,
    const int   nz
    ) noexcept
{
    return  ix + nx*iy + nx*ny*iz;
}


inline
__HOST__ __DEVICE__
int id(
    const int  ids[],
    const int  in,  const int  jn, const int  kn
    ) noexcept
{
    return  ids[ idv(in,jn,kn) ];
}


template<int in, int jn, int kn>
__HOST__ __DEVICE__
int id(const int  ids[]) noexcept
{
    return  ids[ idv<in,jn,kn>() ];
}


template<int NX_LEAF>
__HOST__ __DEVICE__
int  id(
    const int   ix, const int   iy, const int   iz,
    const int   offsets3d[]
    ) noexcept
{
    // inside //
    const int  ix_tmp = (ix+NX_LEAF)%NX_LEAF;
    const int  iy_tmp = (iy+NX_LEAF)%NX_LEAF;
    const int  iz_tmp = (iz+NX_LEAF)%NX_LEAF;
    const int  id_tmp = id<NX_LEAF>(ix_tmp, iy_tmp, iz_tmp);

    // outside //
    const int  nix_tmp = (ix < 0) ?  -1 : (ix < NX_LEAF) ? 0 : 1 ;
    const int  niy_tmp = (iy < 0) ?  -1 : (iy < NX_LEAF) ? 0 : 1 ;
    const int  niz_tmp = (iz < 0) ?  -1 : (iz < NX_LEAF) ? 0 : 1 ;

    const int  offset_tmp  = offsets3d[ idv(nix_tmp, niy_tmp, niz_tmp) ];

    // global offset //
    return  id_tmp + offset_tmp;
}


template<int NX_LEAF>
__HOST__ __DEVICE__
int  id0(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offset0
    ) noexcept
{
    return  id<NX_LEAF>(ix,iy,iz) + offset0;
}


template<int nx, int ny, int nz>
__HOST__ __DEVICE__
int  index(
    const int   i,
    const int   j,
    const int   k
    ) noexcept
{
    return i + nx*j + nx*ny*k;
}

// 7 stencils
inline
__HOST__ __DEVICE__
int  idv_offset7(const int iv, const int jv, const int kv) noexcept
{
    return (iv ==  0 && jv ==  0 && kv ==  0) ? 0 :
           (iv == -1 && jv ==  0 && kv ==  0) ? 1 :
           (iv ==  1 && jv ==  0 && kv ==  0) ? 2 :
           (iv ==  0 && jv == -1 && kv ==  0) ? 3 :
           (iv ==  0 && jv ==  1 && kv ==  0) ? 4 :
           (iv ==  0 && jv ==  0 && kv == -1) ? 5 :
           (iv ==  0 && jv ==  0 && kv ==  1) ? 6 :
                                                0 ;
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
constexpr int  idv_offset7() noexcept
{
    return (iv ==  0 && jv ==  0 && kv ==  0) ? 0 :
           (iv == -1 && jv ==  0 && kv ==  0) ? 1 :
           (iv ==  1 && jv ==  0 && kv ==  0) ? 2 :
           (iv ==  0 && jv == -1 && kv ==  0) ? 3 :
           (iv ==  0 && jv ==  1 && kv ==  0) ? 4 :
           (iv ==  0 && jv ==  0 && kv == -1) ? 5 :
           (iv ==  0 && jv ==  0 && kv ==  1) ? 6 :
                                                0 ;
}


template<int NX_LEAF>
__HOST__ __DEVICE__
int  id_offset7(
    const int   ix, const int   iy, const int   iz,
    const int   offsets7[]
    ) noexcept
{
    // inside //
    const int  ix_tmp = (ix+NX_LEAF)%NX_LEAF;
    const int  iy_tmp = (iy+NX_LEAF)%NX_LEAF;
    const int  iz_tmp = (iz+NX_LEAF)%NX_LEAF;
    const int  id_tmp = id<NX_LEAF>(ix_tmp, iy_tmp, iz_tmp);

    // outside //
    const int  nix_tmp = (ix < 0) ?  -1 : (ix < NX_LEAF) ? 0 : 1 ;
    const int  niy_tmp = (iy < 0) ?  -1 : (iy < NX_LEAF) ? 0 : 1 ;
    const int  niz_tmp = (iz < 0) ?  -1 : (iz < NX_LEAF) ? 0 : 1 ;

    const int  offset_tmp  = offsets7[ idv_offset7(nix_tmp, niy_tmp, niz_tmp) ];

    // global offset //
    return  id_tmp + offset_tmp;
}


template<int NX_LEAF>
__HOST__ __DEVICE__
int  id_offset7_in_leaf(
    const int   ix, const int   iy, const int   iz,
    const int   offsets7[]
    ) noexcept
{
    // global offset //
    return  id<NX_LEAF>(ix, iy, iz) + offsets7[ idv_offset7<0,0,0>() ];
}


};


#endif
