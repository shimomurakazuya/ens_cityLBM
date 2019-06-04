#pragma once
#ifndef STNODEINFORMATIONS_H_
#define STNODEINFORMATIONS_H_


#include "defineAMR.h"


struct  NodeIndices {
    int  lv;
    int  i;
    int  j;
    int  k;

public:
    NodeIndices() {};

    NodeIndices(const int  _lv, const int  _i, const int  _j, const int  _k) :
        lv(_lv),
        i(_i), j(_j), k(_k)
    {};

    void  cout(const std::string  str) const { std::cout << str << lv << ", " << i << ", " << j << ", " << k << std::endl; }
};


namespace  NodeTypes {
    enum  Type {
            RootT               = -1, // does not have parent node //
            ExternalLeafT       = -2, // does not have child node //
            ExternalNeighborT   = -3, // does not have neighbor node //
            InvalidT            = -4, // is invalid node //
            LeafT               =  1, // defalut //
            };
};


struct  MeshInfo {
    int  mesh_offset; // offset is defined for each amr-level //
};


#endif
