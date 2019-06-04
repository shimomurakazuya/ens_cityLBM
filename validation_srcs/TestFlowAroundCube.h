#pragma once
#ifndef TESTFLOWAROUNDCUBE_H_
#define TESTFLOWAROUNDCUBE_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


#if 0
class  TestFlowAroundCube {
public:
    int  rank_;

    TestFlowAroundCube(){ MPI_Comm_rank(MPI_COMM_WORLD, &rank_); }
    ~TestFlowAroundCube(){}

    const real  box_min_[3] = { -0.05, -0.05, -0.10 };
    const real  box_max_[3] = {  0.05,  0.05,  0.10 };

    const real  xyz_min_[3] = { -19.6, -1.05, 0.00 };
    const real  xyz_max_[3] = {  19.6,  1.05, 2.10 };

    const real  bc_layer_min_[3] = { 0.50, 0.10, 0.00 };
    const real  bc_layer_max_[3] = { 0.50, 0.10, 0.25 };

    const real  bc_dirichlet_min_[3] = { 0.02, 0.00, 0.0 };
    const real  bc_dirichlet_max_[3] = { 0.02, 0.00, 0.0 };

    // fin //
    const int   num_fins_  = 9;
    const real  fin_height_ = 0.4;
    const real  fin_width_ = 0.08;
    const real  fin_pos_   = xyz_min_[0] + bc_layer_min_[0] + 0.1;


    // source //
    const real  sc_xyz_[3] = { -0.150, 0.00, 0.050 };
    const real  sc_dscalar_ = 0.1 * 1016.0 / 60.0 * 0.001; //   kg / s  (1cc = 1g)  // 10% //

public:
    void  FlowAroundCube(int argc, char* argv[]);


private:
    void  Init(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  InitValueFlowAroundCube(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitObjectFlowAroundCube(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues,
        const real          box_min[],
        const real          box_max[]
        );

    void  InitBCAroundCube(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitSource(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    real  levelset_from_Inflow_fin(const real lv_obs, const real xyz[]) const;
};


namespace FuncFlowAroundCube {


inline
__HOST__ __DEVICE__
real  inflow_velocity(const real z)
{
//    constexpr real  us = 2.14; // default
//    constexpr real  us = 2.24;
    constexpr real  us = 2.30;
    constexpr real  zs = 0.5;
    constexpr real  alpha = 1.0/7.0;

    const real  ztmp = (z < 2.0) ? z : 2.0;

    return  us * pow( (ztmp/zs), alpha );
}


};
#endif


#endif
