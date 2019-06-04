#pragma once
#ifndef STIONODE_H_
#define STIONODE_H_


#include "stNodeInformations.h"


struct  IONode {
    int   index;

    // connections //
    int   neighbor[6];
    int   parent;
    int   child;

    // NodeInformations //
    NodeTypes::Type  node_type;
    NodeIndices      position;

    // NodeCalFlags //
    bool  nodeCalFlags[7];

    // MeshInfo //
    MeshInfo   meshInfo;
};


#endif
