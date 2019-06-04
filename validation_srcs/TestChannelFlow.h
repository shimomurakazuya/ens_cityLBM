#pragma once
#ifndef TESTCHANNELFLOW_H_
#define TESTCHANNELFLOW_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


class  TestChannelFlow {
public:
    int  rank_;

    TestChannelFlow(){ MPI_Comm_rank(MPI_COMM_WORLD, &rank_); }
    ~TestChannelFlow(){}

    // source //
    const real  sc_dscalar_ = 0.0 / 60.0 * 0.001; //   kg / s  (1cc = 1g)  // 10% //

public:
    void  Flow(int argc, char* argv[]);


private:
    void  Init(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  InitValueChannelFlow(
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

//    void  InitBC(
//        const Grid*         grids,
//        const Tree&         tree,
//        const Parameters&   parameters,
//              MeshValue*    meshValues
//        );

    void  InitSource(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );
};


#endif
