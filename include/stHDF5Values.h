#pragma once
#ifndef STHDF5VALUES_H_
#define STHDF5VALUES_H_


#include <iostream>
#include "HDFTypes.h"


struct  stHDF5Values {
    // coordinates xyz //
    std::unique_ptr<float[]>  xp;
    std::unique_ptr<float[]>  yp;
    std::unique_ptr<float[]>  zp;

    // values //
    std::unique_ptr<float[]>  lv_obj;
    std::unique_ptr<float[]>  rho_obj;
    std::unique_ptr<float[]>  u_obj;
    std::unique_ptr<float[]>  v_obj;
    std::unique_ptr<float[]>  w_obj;

    std::unique_ptr<float[]>  u;
    std::unique_ptr<float[]>  v;
    std::unique_ptr<float[]>  w;
    std::unique_ptr<float[]>  rho;

    std::unique_ptr<float[]>  scalar;


    // lbm values //
    std::unique_ptr<float[]>  f_lbm;

    // HDF5 grid information //
    std::unique_ptr<int[]>    amr_lv;
//    std::unique_ptr<int[]>    amr_posx;
//    std::unique_ptr<int[]>    amr_posy;
//    std::unique_ptr<int[]>    amr_posz;
    std::unique_ptr<int[]>    mpi_rank;
};


#endif
