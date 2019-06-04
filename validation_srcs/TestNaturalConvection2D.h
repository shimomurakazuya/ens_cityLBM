#pragma once
#ifndef TESTNATURALCONVECTION2D_H_
#define TESTNATURALCONVECTION2D_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


class  TestNaturalConvection2D {
public:
    int rank_;
    TestNaturalConvection2D(){ MPI_Comm_rank(MPI_COMM_WORLD, &rank_); }
    ~TestNaturalConvection2D(){}

public:
    void  Flow(int argc, char* argv[]);


private:
    void  Init(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );


    void  Init_Value(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  Init_Object(
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues,
        const real          box_min[],
        const real          box_max[]
        );
};


#endif
