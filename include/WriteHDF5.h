#if 0

#pragma once
#ifndef WRITEHDF5_H_
#define WRITEHDF5_H_


#include <iostream>
#include <fstream>
#include <string>
#include "HDFTypes.h"
#include "stHDF5Values.h"


namespace  WriteHDF5 {
using namespace HDFTypes;


void
write_hdf_datasets (
    const std::string                   folder,
    const std::string                   filename,
    const std::string                   filename_lbm,
    const std::string                   grid_name0,
    const std::vector<stHDF5Values>&    HDF5Values,
    const std::vector<int>&             lv_leaf,
    const std::vector<bool>&            flag_visualization,
    const int                           nx_cell[]
    );


void
write_hdf_dataset_impl (
    const hid_t&  file_id,
    const std::vector<WriteFormat>& vec
    );


void
write_hdf5_and_xdmf_datasets (
    const std::string               folder,
    const std::string               filename,
    const std::string               filename_lbm,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioLBMAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType,
    const std::vector<int>&         lv_leaf,
    const std::vector<bool>&        flag_visualization
    );


void
write_hdf5_dataset (
    const std::string               folder,
    const std::string               filename,
    const std::string               filename_lbm,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioLBMAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType
    );


void
write_hdf5_dataset_lbm (
    const std::string               folder,
    const std::string               filename_lbm,
    const std::vector<WriteFormat>& ioLBMAttributeDataType
    );


void
write_hdf5_dataset_for_visualization (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType
    );


void
write_xdmf_dataset (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioLBMAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType
    );


void
write_xdmf_dataset_for_visualization (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType,
    const std::vector<bool>&        flag_visualization
    );


void
write_xdmf_dataset_for_lv (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType,
    const std::vector<int>&         lv_leaf,
    const int                       lv
    );


};


#endif

#endif
