#pragma once
#ifndef TESTNATURALCONVECTION3D_H_
#define TESTNATURALCONVECTION3D_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


class  TestNaturalConvection3D {
public:
    int rank_;
    TestNaturalConvection3D() { 
        //MPI_Comm_rank_(MPI_COMM_WORLD, &rank_); 
        rank_ = -1;

        // heated: 420.0 K //
        TemperatureInit_     = 323.0;
        TemperatureHeated_   = 420.0;
        TemperatureSideWall_ = 312.0;
        TemperatureTopWall_  = 309.5;
        TemperatureSideWall_ax_ = -17.131;
        TemperatureSideWall_b_  = 320.27;

//        // heated: 476.0 K //
//        TemperatureInit_     = 343.0;
//        TemperatureHeated_   = 476.0;
//        TemperatureSideWall_ = 327.0;
//        TemperatureTopWall_  = 323.0;
//        TemperatureSideWall_ax_ = -24.924;
//        TemperatureSideWall_b_  = 340.62;
//
//        // heated: 573.0 K //
//        TemperatureInit_     = 380.0;
//        TemperatureHeated_   = 573.0;
//        TemperatureSideWall_ = 357.0;
//        TemperatureTopWall_  = 343.0;
//        TemperatureSideWall_ax_ = -58.436;
//        TemperatureSideWall_b_  = 387.7;

    }
    ~TestNaturalConvection3D(){}

    real TemperatureInit_;
    real TemperatureHeated_;
    real TemperatureSideWall_;
    real TemperatureTopWall_;
    real TemperatureSideWall_ax_;
    real TemperatureSideWall_b_;

    const real z_min_{0.000};
    const real z_max_{0.800};

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
