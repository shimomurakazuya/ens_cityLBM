#pragma once
#ifndef ARRAY3D_H_
#define ARRAY3D_H_


#include "defineLBM.h"


template <typename T>
struct Array3D {
    T val[3][3][3];

public:
    T offset0() const { return val[1][1][1]; }
    T offset(const int i, const int j, const int k) const { return val[k+1][j+1][i+1]; }

    template<int i, int j, int k>
    T offset() const { return val[k+1][j+1][i+1]; }


    // lbm //
    T offset_lbm0() const { return val[1][1][1] * (LBM_velocity_model::nQ); }
    T offset_lbm(const int i, const int j, const int k) const { return val[k+1][j+1][i+1] * (LBM_velocity_model::nQ); }

    template<int i, int j, int k>
    T offset_lbm() const { return val[k+1][j+1][i+1] * (LBM_velocity_model::nQ); }


    Array3D<T>  offset3d_lbm()
    const
    {
        constexpr int nQ = LBM_velocity_model::nQ;

        return  Array3D<T>
        {
            val[0][0][0]*nQ, val[0][0][1]*nQ, val[0][0][2]*nQ,
            val[0][1][0]*nQ, val[0][1][1]*nQ, val[0][1][2]*nQ,
            val[0][2][0]*nQ, val[0][2][1]*nQ, val[0][2][2]*nQ,

            val[1][0][0]*nQ, val[1][0][1]*nQ, val[1][0][2]*nQ,
            val[1][1][0]*nQ, val[1][1][1]*nQ, val[1][1][2]*nQ,
            val[1][2][0]*nQ, val[1][2][1]*nQ, val[1][2][2]*nQ,

            val[2][0][0]*nQ, val[2][0][1]*nQ, val[2][0][2]*nQ,
            val[2][1][0]*nQ, val[2][1][1]*nQ, val[2][1][2]*nQ,
            val[2][2][0]*nQ, val[2][2][1]*nQ, val[2][2][2]*nQ,
        };
    }
};


#endif
