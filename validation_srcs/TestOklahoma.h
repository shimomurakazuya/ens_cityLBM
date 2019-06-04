#pragma once
#ifndef TESTOKLAHOMA_H_
#define TESTOKLAHOMA_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"
#include "defineFluidProperty.h"


class  TestOklahoma {
public:
    int rank_;
    TestOklahoma() {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }
    ~TestOklahoma(){}

    const real TemperatureGround_{ fluid_property::Temperature0 };
    const real Zmin_{ 0.0 };

    const real  bc_layer_min_[3] = { 64.0, 64.0,   0.0 };
    const real  bc_layer_max_[3] = { 64.0, 64.0, 256.0 };
//    const real  bc_layer_min_[3] = { 256.0, 256.0,   0.0 };
//    const real  bc_layer_max_[3] = { 256.0, 256.0, 256.0 };
//    const real  bc_layer_min_[3] = { 512.0, 512.0,   0.0 };
//    const real  bc_layer_max_[3] = { 512.0, 512.0, 512.0 };
//    const real  bc_layer_min_[3] = { 768.0, 768.0,   0.0 };
//    const real  bc_layer_max_[3] = { 768.0, 768.0, 512.0 };

    const real  bc_dirichlet_min_[3] = { 256.0, 256.0,    0.0 };
    const real  bc_dirichlet_max_[3] = { 256.0, 256.0, 1024.0 };

    int  start_minutes_{ 360 };

public:
    void  Flow(int argc, char* argv[]);


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

    void  InitValueOklahoma(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitObjectFloarRoof(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  DigitalizeObject(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitObjectMapTokyo(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitObjectMapOklahoma(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitDecoBoco(
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


    void  InitObservationDataOklahoma(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  InitObservationData(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  InitObservationData2(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  _InitObservationData(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  ReadOklahoma_uvwTE(float* u, float* v, float* w, float* T, float* E,
            int nx, int ny, int nz,
            int time_num);

    void  Refinement_uvwTE(float* uf, float* vf, float* wf, float* Tf, float* Ef,
            const float* u, const float* v, const float* w, const float* T, const float* E,
            int nx, int ny, int nz, int cx, int cy, int cz);

};


namespace FuncOklahoma {


inline
__HOST__ __DEVICE__
real  inflow_velocity(const real z)
{
    const real umax = 10.0;

    const real zH = 100.0;
    const real z0 =  2.0;
    const real A  =  umax / log10( zH/z0 );

    const real  ztmp = (z < zH) ? z : zH;

    return A * log10( ztmp/z0 );
}


};


#endif
