#pragma once
#ifndef STNODEVALARRAY_H_
#define STNODEVALARRAY_H_


#include "Array3D.h"


struct  stNodeValArray {
    int*            id_parent_;
    int*            id_child_;

    int*            mesh_offsets_;
    int*            mesh_offsets3x3x3_;
    Array3D<int>*   mesh_offsets3d_;

    int*            node_pos_x_;
    int*            node_pos_y_;
    int*            node_pos_z_;
};


#endif
