#pragma once
#ifndef TESTCHANNELFLOS_WITH_HF_H_
#define TESTCHANNELFLOS_WITH_HF_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"


class  TestChannelFlow_with_HF {
public:
    int rank_;
    TestChannelFlow_with_HF() { 
        //MPI_Comm_rank_(MPI_COMM_WORLD, &rank_); 
        rank  = -1;

        // heated: 420.0 K //
        TemperatureInit_     = 0.5;
        TemperatureBotoomWall = 0.0;
        TemperatureSideWall_  = 0.5;
        TemperatureTopWall_   = 1.0;

    }
    ~TestChannelFlow_with_HF(){}

    real TemperatureInit_;
    real TemperatureBotoomWall;
    real TemperatureSideWall_;
    real TemperatureTopWall_;

    const real z_min_{-1.0};
    const real z_max_{ 1.0};

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

    void  Init_ValueChannelFlow(
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
