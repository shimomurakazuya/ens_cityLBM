#pragma once
#ifndef FUNCMAPDATA_H_
#define FUNCMAPDATA_H_


#include <fstream>
#include <string>
#include <cmath>
#include "definePrecision.h"


namespace  FuncMapData {


void  levelset_map(
          real* lv_obs,
    const real* x,
    const real* y,
    const real* z,
    const int   nx_leaf,
    const int   lv,
    const real  dx,
    const int*  height_map,
    const int   mx,
    const int   my,
    const real  dx_map,
    const real  map_offset_x,
    const real  map_offset_y,
    const real  domain_min_x,
    const real  domain_min_y,
    const real  domain_max_x,
    const real  domain_max_y
    );


double  get_levelset(
    const double xc,
    const double yc,
    const double zc,
    const int    lv,
    const double dx,
    const int*   height_map,
    const int    mx,
    const int    my,
    const double dx_map,
    const double map_offset_x,
    const double map_offset_y,
    const double domain_min_x,
    const double domain_min_y,
    const double domain_max_x,
    const double domain_max_y
    );


int  get_height(
    const int* height_map,
    const int mx,
    const int my,
    const int ix,
    const int iy
    );


bool  check_index_in_map(
    const int  ix,
    const int  iy,
    const int  nx,
    const int  ny
    );


void  convert_geometory_to_index(
          int& ix,
          int& iy,
    const double x,
    const double y,
    const double dx,
    const double map_offset_x,
    const double map_offset_y
    );


};


#endif
