#pragma once
#ifndef TESTCAVITYFLOW_H_
#define TESTCAVITYFLOW_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


class  TestCavityFlow {
public:
    int rank_;
    TestCavityFlow(){ MPI_Comm_rank(MPI_COMM_WORLD, &rank_); }
    ~TestCavityFlow(){}

public:
    void  CavityFlow2D(int argc, char* argv[]);


private:
    void  Init(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  Init_Value_Cavity2D(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  Init_Object_Cavity2D(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues,
        const real          box_min[],
        const real          box_max[]
        );
};


#endif
