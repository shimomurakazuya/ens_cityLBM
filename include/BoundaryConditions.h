#pragma once
#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_


#include <iostream>
#include <string>
#include "definePrecision.h"
#include "defineLBM.h"
#include "Field.h"


class  BoundaryConditions {
private:
    int  rank_;

public:
    BoundaryConditions(){
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    ~BoundaryConditions(){}

public:
    void  InitializeByWRFData(
        const int           time_num,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  InitializeByGroundData(
        const int           time_num,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  ReadWRFData(
        const int           time_num,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

    void  ReadGroundData(
        const int           time_num,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues
        );

private:
    void  ReadOklahoma_uvwTE(float* u, float* v, float* w, float* T, float* E, int nx, int ny, int nz, int time_num);
    void  ReadOklahoma_GroundTb(float* Tb, int nx, int ny, int time_num);

    real  f_interpolate2d(const float* f, float xup, float yup, int i, int j, int nx, int ny);
    real  f_interpolate3d(const float* f, float xup, float yup, float zup, int i, int j, int k, int nx, int ny, int nz);

};


#endif
