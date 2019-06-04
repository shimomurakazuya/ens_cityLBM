#pragma once
#ifndef STIOGRID_H_
#define STIOGRID_H_


#include "definePrecision.h"


struct  stIOGrid {
    int  nx;
    int  ny;
    int  nz;

    real offset_x;
    real offset_y;
    real offset_z;

    real dx;
};


#endif
