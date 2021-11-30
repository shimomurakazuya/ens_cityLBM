#pragma once
#ifndef MPIDATATYPEINFO_H_
#define MPIPATATYPEINFO_H_


#include <mpi.h>
#include <iostream>
#include <vector>
#include "defineMemory.h"


struct  v3x3x3 {
    bool v[3][3][3];

    void zero()
    {
        for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
        for (int iv=-1; iv<=1; iv++) {
                v[iv+1][jv+1][kv+1] = false;
        }
        }
        }
    }
};


struct  MPIDatatypeInfo {
    MPI_Datatype  mpi_datatype;

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

    std::vector<int>    id_leaf;
    std::vector<v3x3x3> leaf_type; // true : has neighbor send leaf
};


struct PackUnpackVal {
    std::vector<int>  num_slist_val_rank;
    std::vector<int>  num_rlist_val_rank;

    std::vector<int>  offset_slist_val_rank; // sum of num_slist_val_rank
    std::vector<int>  offset_rlist_val_rank;
    int   num_slist_val; 
    int   num_rlist_val;
    int*  slist_val = nullptr;
    int*  rlist_val = nullptr;
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
    int*  slist_val = nullptr;
    int*  rlist_val = nullptr;

    int   num_slist_lbm;
    int   num_rlist_lbm;
    int*  slist_lbm = nullptr;
    int*  rlist_lbm = nullptr;

    PackUnpackVal  all;
    PackUnpackVal  slice1;
    PackUnpackVal  slice2;
    PackUnpackVal  slice3;

    PackUnpackVal  all_lbm;
    PackUnpackVal  slice1_lbm;
    PackUnpackVal  slice2_lbm;
    PackUnpackVal  slice3_lbm;
};


#endif
