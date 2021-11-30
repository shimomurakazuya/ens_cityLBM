#pragma once
#ifndef INITTREE_FOR_MPI_H_
#define INITTREE_FOR_MPI_H_


#include "Tree.h"
#include "Grid.h"
#include "MPIDatatypeInfo.h"
#include <mpi.h>
#include <vector>
#include "MPICommEnsemble.h"


struct  LocalNode {
    // basic information //
    int  lv;
    int  mesh_offset;

    int  id;
    int  idg;
    int  cal_rank;
    bool cal_flag;
    bool allocate_flag;
    bool halo_flag;

    // connection information //
    bool flag_put; // @ cal_flag == true //
    std::vector<int>  rank_put;
    std::vector<int>  id_put  ;
    std::vector<int>  offset_put  ;

    bool flag_get; // @ cal_flag == false && allocate_flag == true //
    int  rank_get;
    int  id_get;
    int  offset_get;
};


struct  _MPIGetInfo {
    int  gindex; // global index //
    int  rank_dst;
};


class InitTree_for_MPI {
private:
    const MPICommEnsemble comm_;
    int     num_allocated_localNode_;
    std::vector<LocalNode>  localNode_;

public:
    InitTree_for_MPI() = delete;

    InitTree_for_MPI( const MPICommEnsemble comm): comm_(comm)
    {
    }
    ~InitTree_for_MPI() {}


public:
    void
    init_tree_uniform3d_div(
              Tree& tree,
        const Grid* grids,
        const int   lv_max
        ) = delete;

    void init_tree_mpi_with_map(
              Tree& tree,
        const Grid* grids,
        const int   lv_max,
        const MapData& map
        );
              


private:
    void
    set_LocalNode(
        const Tree& tree_global,
        const Grid* grids
        );

    void
    set_LocalTree(
              Tree& tree,
        const Tree& tree_global,
        const Grid* grids
        );

    void createMPIPutGetInfo(const Tree& tree_global, std::vector<MPIPutGetInfo>& mpiPutGetInfo) const;

private:
    void set_lv(std::vector<LocalNode>& localNode, const Tree& tree_global);
    void set_cal_rank_div(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div1d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div2d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div3d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div_nc3d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div_citybox(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div_channel(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids);
    void set_cal_rank_div_2d_block(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids, int bx=2, int by=2);
    void set_flags   (std::vector<LocalNode>& localNode, const Tree& tree_global);

    bool node_is_cal(const int cal_rank, const int rank);
    bool neighbor_is_cal(std::vector<LocalNode>& localNode, const Node* node);
    bool neighbor_is_cal2(std::vector<LocalNode>& localNode, const Node* node);

    void set_LocalNode_index(std::vector<LocalNode>& localNode, int& num_allocated_localNode);

    void set_MPIGetInfo(std::vector<LocalNode>& localNode);
    void set_MPIPutInfo(std::vector<LocalNode>& localNode, const Tree& tree_global);


    void
    commit_MPI_struct(
              MPI_Datatype&            mpi_Datatype,
              void*&                   base_ptr,
              std::vector<int>&        val,
        const std::vector<_MPIGetInfo>& mpiGetInfos,
        const int                      istart,
        const int                      count
        );

};


#endif
