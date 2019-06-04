#pragma once
#ifndef INDEXT_H_
#define INDEXT_H_


#include "Array3D.h"


namespace IndexT {


inline
int  id(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx_leaf
    )
{
    return  ix + nx_leaf*iy + nx_leaf*nx_leaf*iz;
}


inline
int  id(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx_leaf,
    const Array3D<int>&  offsets3d
    )
{
    // inside //
    const int  ix_tmp = (ix+nx_leaf)%nx_leaf;
    const int  iy_tmp = (iy+nx_leaf)%nx_leaf;
    const int  iz_tmp = (iz+nx_leaf)%nx_leaf;
    const int  id_tmp = id(ix_tmp, iy_tmp, iz_tmp, nx_leaf);

    // outside //
    const int  nix_tmp = (ix < 0) ?      -1 :
                         (ix < nx_leaf) ? 0 :
                                          1 ;

    const int  niy_tmp = (iy < 0) ?      -1 :
                         (iy < nx_leaf) ? 0 :
                                          1 ;

    const int  niz_tmp = (iz < 0) ?      -1 :
                         (iz < nx_leaf) ? 0 :
                                          1 ;

    const int  offset_tmp  = offsets3d.offset(nix_tmp, niy_tmp, niz_tmp);

    // global offset //
    return  id_tmp + offset_tmp;
}


};


#endif
