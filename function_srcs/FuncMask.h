#pragma once
#ifndef FUNCMASK_H_
#define FUNCMASK_H_


#include <iostream>
#include "MPIDatatypeInfo.h"


namespace FuncMask {

template<int NX_LEAF>
bool mask_grid(
    const int i,  const int j,  const int k, // grid in leaf
    const int nlayer,
    const int iv, const int jv, const int kv // direction : iv,jv,kv = -1,0,1 (-1:leaf, 1:right)
    )
{
    const bool flag_x = (iv == -1) ? (i < nlayer)
                      : (iv ==  1) ? (i > NX_LEAF-1 - nlayer)
                      :              true;

    const bool flag_y = (jv == -1) ? (j < nlayer)
                      : (jv ==  1) ? (j > NX_LEAF-1 - nlayer)
                      :              true;

    const bool flag_z = (kv == -1) ? (k < nlayer)
                      : (kv ==  1) ? (k > NX_LEAF-1 - nlayer)
                      :              true;

    return ( flag_x && flag_y ) && flag_z;
}


template<int NX_LEAF>
void mask_leaf(
    std::vector<bool>& flag,
    const int nlayer,
    const v3x3x3&  v3d
    )
{
    for (int k=0; k<NX_LEAF; k++) {
    for (int j=0; j<NX_LEAF; j++) {
    for (int i=0; i<NX_LEAF; i++) {
        const int id = i + NX_LEAF*j + NX_LEAF*NX_LEAF*k;

        for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
        for (int iv=-1; iv<=1; iv++) {
            if (v3d.v[iv+1][jv+1][kv+1]) {
                const bool flag_tmp = mask_grid<NX_LEAF>(i,j,k, nlayer, iv,jv,kv);
                flag[id] = ( flag[id] || flag_tmp );
//                flag[id] = true;
            }
        }
        }
        }
    }
    }
    }
}


}


#endif
