#if 0

#pragma once
#ifndef READHDF5_H_
#define READHDF5_H_


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <vector>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "defineLBM.h"
#include "HDFTypes.h"
#include "stHDF5Values.h"


namespace  ReadHDF5 {
using namespace HDFTypes;


template <typename T>
struct  ReadHDF5ptr {
    T*  ptr;
    int n;
};


template <typename T>
T*
read_hdf5_dataset(
    const int           ndim,
    const hid_t         file_id,
    const std::string   dset_name,
    const hid_t         type_id
    )
{
    hsize_t  dims[3];
    H5LTget_dataset_info(file_id, ("/" + dset_name).c_str(), dims, NULL, NULL);
    const int  nn = dims[0]*dims[1]*dims[2];

    T*  val = new T[nn];
    H5LTread_dataset(file_id, ("/" + dset_name).c_str(), type_id, val);

    return  val;
}


template <typename T>
T*
read_hdf5_dataset_series(
    const int           ndim,
    const hid_t         file_id,
    const std::vector<std::string>  dset_names,
    const hid_t         type_id
    )
{
    const int  nvec = dset_names.size();

    hsize_t  dims[3];
    H5LTget_dataset_info(file_id, ("/" + dset_names[0]).c_str(), dims, NULL, NULL);
    const int  nn = dims[0]*dims[1]*dims[2];

    T*  val = new T[nn*nvec];
    for (int i=0; i<nvec; i++) {
        const int  offset = nn*i;
        H5LTread_dataset(file_id, ("/" + dset_names[i]).c_str(), type_id, &val[offset]);
    }

    return  val;
}


stHDF5Values
read_hdf_datasets(
    const hid_t file_id,
    const hid_t file_id_lbm,
    const int   i
    );


};


#endif

#endif
