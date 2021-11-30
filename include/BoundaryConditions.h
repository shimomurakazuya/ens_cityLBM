#pragma once
#ifndef BOUNDARYCONDITIONS_H_
#define BOUNDARYCONDITIONS_H_


#include <iostream>
#include <string>
#include "definePrecision.h"
#include "defineLBM.h"
#include "Field.h"
#include "MPICommEnsemble.h"
#include "ParticleFilterSt.h"


class  BoundaryConditions {
private:
    const MPICommEnsemble comm_;

public:
    BoundaryConditions() = delete;
    BoundaryConditions(const MPICommEnsemble comm): comm_(comm) {}
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

    void  InitializeHeatFluxObj(
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
              MeshValue*    meshValues,
        const float         coef_nuding
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
    void  ReadOklahoma_GroundTflux(float* Tflux, int nx, int ny, int time_num);

    void  ReadHeatFluxObjHeader(
        const std::string& fname,
        int&   nx_da,    int&   ny_da,
        float& x_min_da, float& y_min_da,
        float& pitch_x,  float& pitch_y
        );

    void  ReadHeatFluxObj(
        const std::string& fname,
        float* hflux,
        int nx_da, int ny_da
        );

    real  f_interpolate2d(const float* f, float xup, float yup, int i, int j, int nx, int ny);
    real  f_interpolate3d(const float* f, float xup, float yup, float zup, int i, int j, int k, int nx, int ny, int nz);

    real rand(int seed, real u0) const;
};


#endif
