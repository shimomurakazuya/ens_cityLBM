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
#include "MPICommEnsemble.h"
#include "mpi_wrapper.hpp"


class  IOData {
private:
    const MPICommEnsemble comm_;
    const int  step_;

    const std::string  hdf5_grid_name_ = "hdf5-rank" + std::to_string(comm_.world().rank()) + "-grid";

public:
    IOData (const MPICommEnsemble comm, const int step) : comm_(comm), step_(step) {}
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
//

    void writeBinaries(
        const Grid*          grids,
        const Tree&          tree,
        const MeshValue*     meshValues0,
        const MeshValue*     meshValues1,
        const ValueTimeAverage* valueTimeAverage
    ) const;

    void readBinaries(
        const Grid*          grids,
        const Tree&          tree,
              MeshValue*     meshValues0,
              MeshValue*     meshValues1,
              ValueTimeAverage* valueTimeAverage
    ) const;


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
        const VTKOutputScale vtkOutputScale = VTKOutputScale::Full,
        const int            filter_bits = 8,
        const std::vector<real>& zSlices = std::vector<real>(),
        const std::vector<real>& ySlices = std::vector<real>(),
        const std::vector<real>& xSlices = std::vector<real>()
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
