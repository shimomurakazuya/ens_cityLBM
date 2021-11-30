#pragma once
#ifndef TESTMAPFLOW_H_
#define TESTMAPFLOW_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


class  TestMapFlow {
public:
    int  rank_;

    TestMapFlow(){ 
        //MPI_Comm_rank_(MPI_COMM_WORLD, &rank_);
        rank_ = -1;
    }
    ~TestMapFlow(){}

    const real  bc_layer_min_[3] = { 64.0, 64.0,   0.0 };
    const real  bc_layer_max_[3] = { 64.0, 64.0, 128.0 };

    const real  bc_dirichlet_min_[3] = { 16.0, 16.0, 0.0 };
    const real  bc_dirichlet_max_[3] = { 16.0, 16.0, 16.0 };


    // source //
    const real  sc_dscalar_ = 1016.0 / 60.0; //   kg / s  (1cc = 1g)  // 10% //

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

    void  InitObjectFloarRoof(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitObjectMap(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitBC(
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


    void  ReconstructLevelset(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

};


namespace FuncMapFlow {


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
