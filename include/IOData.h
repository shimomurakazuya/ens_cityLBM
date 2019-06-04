#pragma once
#ifndef IODATA_H_
#define IODATA_H_


#include <iostream>
#include <cstdlib>
#include "Field.h"
#include "defineFilenames.h"
#include "defineAMR.h"
#include "defineLBM.h"
//#include "stHDF5Values.h"
#include "Array3D.h"
#include "ParaviewVTU.h"


class  IOData {
private:
    const int  rank_;
    const int  step_;

    const std::string  hdf5_grid_name_ = "hdf5-rank" + std::to_string(rank_) + "-grid";

public:
    IOData (const int rank, const int step) : rank_(rank), step_(step) {}
    ~IOData () {}

public:
//    void
//    writeHDF5Values (
//        const Grid*         grids,
//        const Tree&         tree,
//        const MeshValue*    meshValues
//        )
//    const;
//
//
//    void
//    readHDF5Values (
//        const Grid*         grids,
//        const Tree&         tree,
//              MeshValue*    meshValues
//        )
//    const;


#if 0
    void
    writeVTKFile (
        const Grid*         grids,
        const Tree&         tree,
        const MeshValue*    meshValues
        )
    const;
#endif

    void
    writeVTKFile (
        const Grid*          grids,
        const Tree&          tree,
        const Parameters&    parameters,
        const MeshValue*     meshValues,
        const VTKOutputScale vtkOutputScale = VTKOutputScale::Full
        )
    const;


private:
//    stHDF5Values
//    makeHDF5Values (
//        const int           nx,
//        const int           ny,
//        const int           nz,
//        const Array3D<int>& offsets,
//        const Array3D<int>& offsets_lbm,
//        const Node*         node,
//        const MeshValue&    meshValue
//        )
//    const;

};


#endif
