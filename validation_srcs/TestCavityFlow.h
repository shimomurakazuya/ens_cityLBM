#pragma once
#ifndef TESTCAVITYFLOW_H_
#define TESTCAVITYFLOW_H_


#include <iostream>
#include <string>
#include <fstream>
#include "defineCUDA.h"
#include "Field.h"
#include "defineLBM.h"
#include "MPICommEnsemble.h"

class  TestCavityFlow {
    const MPICommEnsemble comm_;
    const OptionParser& opt_;

public:
    int rank_;
    TestCavityFlow(const OptionParser& opt, const MPICommEnsemble comm): opt_(opt), comm_(comm) {} 
    ~TestCavityFlow(){}

public:
    void  Flow();


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
