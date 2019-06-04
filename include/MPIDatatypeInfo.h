#pragma once
#ifndef MPIDATATYPEINFO_H_
#define MPIPATATYPEINFO_H_


#include <mpi.h>
#include <iostream>
#include <vector>
#include "defineMemory.h"


struct  MPIDatatypeInfo {
    MPI_Datatype  mpi_datatype;

//    void*         base_ptr;
    int           rank;
    int           count;
    int           blength;
    MPI_Datatype  dtype;

    std::vector<int>  id; // mesh id //
};


struct MPIPutGetInfo {
    bool flag_put;
    int  num_vec_put;
    std::vector<int>  rank_put;
    std::vector<int>  id_put  ;
    std::vector<int>  offset_put  ;

    bool flag_get;
    int  rank_get;
    int  id_get; // node id //
    int  offset_get; // mesh_offset : NN_LEAF //
};


struct MPIPutInfo {
    int  id; // src : node id //
    int  rank_put;
    int  id_put  ; // dst : node id //
    int  offset_put; // dst : mesh_offset //
};


struct  MPIDatatypeInfoSrcDst {
    MPIDatatypeInfo  mpiDatatypeInfo_src;
    MPIDatatypeInfo  mpiDatatypeInfo_dst;
};


// pack unpack //
struct MPIPackUnpackCommInfo {
    int  rank;
    int  tag;
    int  nleaf;
    int  sum_leaf;

    std::vector<int>  id_leaf;
};


struct MPIPackUnpackInfo{
    std::vector<MPIPackUnpackCommInfo>  sendPackUnpackCommInfo;
    std::vector<MPIPackUnpackCommInfo>  recvPackUnpackCommInfo;


    // pack unpack on cpu/gpu //
    int   num_slist_leaf;
    int   num_rlist_leaf;
    std::vector<int>  slist_leaf;
    std::vector<int>  rlist_leaf;

    int   num_slist_val;
    int   num_rlist_val;
    int*  slist_val;
    int*  rlist_val;

    int   num_slist_lbm;
    int   num_rlist_lbm;
    int*  slist_lbm;
    int*  rlist_lbm;
};


#endif
