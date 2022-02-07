#pragma once
#ifndef INITTREE_H_
#define INITTREE_H_


#include "Tree.h"
#include "Grid.h"
#include "MPICommEnsemble.h"


namespace  InitTree {


struct  id_amr_init {
    int  id;
    int  lv;
    bool flag_allocate;
    bool flag_cal;
//    int  flag_allocate;
//    int  flag_cal;
};


void
init_tree_uniform3d(
          Tree& tree,
    const Grid* grids,
    const int   lv_max
    );


void
init_tree_amr(
          Tree& tree,
    const Grid* grids,
    const int   lv_max,
    const MapData& map
    );


Tree
tree_all_region(
    const MPICommEnsemble comm,
    const Grid* grids,
    const int   lv_max
    );


std::vector<id_amr_init>
create_id_amr_init(
    const Tree& tree_org,
    const Grid* grids,
    const int   lv_max,
    const std::vector<int> cal_flags
    );


void update_basic_informations(
          Tree& tree,
    const Tree& tree_all,
    const Grid* grids,
    const int   lv_max,
    const std::vector<id_amr_init>&  id_input
    );


void update_connections(
          Tree& tree,
    const Tree& tree_all,
    const Grid* grids,
    const int   lv_max,
    const std::vector<id_amr_init>&  id_input
    );


void update_NodeCalFlags(Tree& tree);


bool check_neighbor(
    int id0, 
    int i, int j, int k,
    int imin, int jmin, int kmin, 
    int imax, int jmax, int kmax,
    int nx, int ny, int nz, 
    std::vector<id_amr_init>& ids
    );


};


#endif
