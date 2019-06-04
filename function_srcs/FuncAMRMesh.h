#pragma once
#ifndef FUNCAMRMESH_H_
#define FUNCAMRMESH_H_


#include <iostream>
#include "defineAMR.h"
#include "Array3D.h"
#include "IndexT.h"


namespace  FuncAMRMesh {


inline
Array3D <int>
stencil3d(
    const int  i,
    const int  j,
    const int  k,
    const Array3D<int>& offsets
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;

    Array3D <int>  stencil3d;
    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                stencil3d.val[kv+1][jv+1][iv+1] = IndexT::id(i+iv,j+jv,k+kv, nx_leaf, offsets);
            }
        }
    }
    return  stencil3d;
}


template <typename T>
T  RawData(
    const T*    val,
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx_leaf,
    const Array3D<int>&  mesh_offsets
    )
{
    return  val[ IndexT::id(ix, iy, iz, nx_leaf, mesh_offsets) ];
}


template <typename T>
T  NodeToCell(
    const T*    val,
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx_leaf,
    const Array3D<int>&  mesh_offsets
    )
{
    T  tmp = 0;
    for (int k=0; k<2; k++) {
        for (int j=0; j<2; j++) {
            for (int i=0; i<2; i++) {
                tmp += val[ IndexT::id(ix+i, iy+j, iz+k, nx_leaf, mesh_offsets) ];
            }
        }
    }
    tmp /= (T)8;

    return  tmp;
}


template <typename T>
T  CellToNode(
    const T*    val,
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx_leaf,
    const Array3D<int>&  mesh_offsets
    )
{
    T  tmp = 0;
    for (int k=0; k<2; k++) {
        for (int j=0; j<2; j++) {
            for (int i=0; i<2; i++) {
                tmp += val[ IndexT::id(ix+i-1, iy+j-1, iz+k-1, nx_leaf, mesh_offsets) ];
            }
        }
    }
    tmp /= (T)8;

    return  tmp;
}


};


#endif
