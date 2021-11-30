#pragma once
#ifndef STNODEVALARRAY_H_
#define STNODEVALARRAY_H_


#include "Array3D.h"


struct  stNodeValArray {
    int*            id_parent_ = nullptr;
    int*            id_child_  = nullptr;

    int*            mesh_offsets_      = nullptr;
    int*            mesh_offsets3x3x3_ = nullptr;
    Array3D<int>*   mesh_offsets3d_    = nullptr;

    int*            node_pos_x_ = nullptr;
    int*            node_pos_y_ = nullptr;
    int*            node_pos_z_ = nullptr;
};


#endif
