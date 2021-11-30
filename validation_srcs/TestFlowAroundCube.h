#pragma once
#ifndef TESTFLOWAROUNDCUBE_H_
#define TESTFLOWAROUNDCUBE_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"
#include "defineFluidProperty.h"
#include "FuncMapData.h"
#include "MPICommEnsemble.h"

class  TestFlowAroundCube {
    const MPICommEnsemble comm_;
    const OptionParser& opt_;

public:

    TestFlowAroundCube() = delete;

    TestFlowAroundCube(const OptionParser& opt, const MPICommEnsemble comm): opt_(opt), comm_(comm) {}
    ~TestFlowAroundCube(){}

    const real TemperatureGround_{ fluid_property::Temperature0 };
    const real Zmin_{ 0.0 };

    const real  box_min_[3] = { -0.05, -0.05,  0.00 };
    const real  box_max_[3] = {  0.05,  0.05,  0.10 };
 
    const real  xyz_min_[3] = { -20.0, -0.96, -0.08 };
    const real  xyz_max_[3] = {  20.0,  0.96,  1.80 };
 
    const real  bc_layer_min_[3] = { 0.50, 0.10, 0.00 };
    const real  bc_layer_max_[3] = { 0.50, 0.10, 0.25 };
 
    const real  bc_dirichlet_min_[3] = { 0.02, 0.00, 0.0 };
    const real  bc_dirichlet_max_[3] = { 0.02, 0.00, 0.0 };
 
    // fin //
    const int   num_fins_  = 9;
    const real  fin_height_ = 0.4;
    const real  fin_width_ = 0.08;
    //const real  fin_width_ = 0.16;
    const real  fin_pos_   = xyz_min_[0] + bc_layer_min_[0] + 5.0;

    // source //
    //const real  sc_xyz_[3] = { -0.150, 0.00, 0.050 }; 
    const real  sc_dscalar_ = 0.1 * 1016.0 / 60.0 * 0.001; //   kg / s  (1cc = 1g)  // 10% //

    int start_minutes0_ { -99999 };
    int  start_minutes_restart_{ -99999 };

    MapData mapData_;

public:
    void  Flow();


private:
    void  Init(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitValue(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitObjectFloarRoof(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues,
        const real          box_min[],
        const real          box_max[]
        );

    void  InitSource(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    real levelset_from_Inflow_fin(const real lv_obs, const real xyz[]) const;

    void  initBCAroundCube(
         const Grid*         grids,
         const Tree&         tree,
         const Parameters&   parameters,
               MeshValue*    meshValues
         );
};


namespace FuncFlowAroundCube {


inline
__HOST__ __DEVICE__
real  inflow_velocity(const real z)
{
    constexpr real  us = 2.30;
    constexpr real  zs = 0.5;
    constexpr real  alpha = 1.0/7.0;

    const real  ztmp = (z < 2.0) ? z : 2.0;
    //const real  ztmp = z;

    return  us * pow( (ztmp/zs), alpha );
}


};
#endif
