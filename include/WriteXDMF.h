#pragma once
#ifndef WRITE_XDMF_H_
#define WRITE_XDMF_H_


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "HDFTypes.h"


namespace  WriteXDMF {
using namespace HDFTypes;


void  write_xdmf_header    (FILE*  xmf);
void  write_xdmf_footer    (FILE*  xmf);


void
write_xdmf_datas (
          FILE*                     xmf,
    const std::vector<WriteFormat>& geometries,
    const std::vector<WriteFormat>& attributes,
    const std::vector<WriteFormat>& lbm_attributes,
    const std::vector<WriteFormat>& grid_attributes
    );


void
write_xdmf_datas_for_visualization (
          FILE*                     xmf,
    const std::vector<WriteFormat>& geometries,
    const std::vector<WriteFormat>& attributes,
    const std::vector<WriteFormat>& grid_attributes,
    const std::vector<bool>&        flag_visualization
    );


void
write_xdmf_datas_for_lv (
          FILE*                     xmf,
    const std::vector<WriteFormat>& geometries,
    const std::vector<WriteFormat>& attributes,
    const std::vector<WriteFormat>& grid_attributes,
    const std::vector<int>&         lv_leaf,
    const int                       lv
    );


void
write_xdmf_geometry_x_y_z (
          FILE*                     xmf,
    const std::vector<WriteFormat>& vec,
    const int                       index,
    const int                       dim
    );


void
write_xdmf_attribute (
          FILE*                     xmf,
    const std::vector<WriteFormat>& vec,
    const int                       index,
    const int                       dim
    );


void
write_xdmf_grid_attribute (
          FILE*                     xmf,
    const std::vector<WriteFormat>& vec,
    const int                       index,
    const int                       dim
    );


};


#endif
