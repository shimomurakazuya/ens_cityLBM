#pragma once
#ifndef TREE_H_
#define TREE_H_


#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>

#include "defineAMR.h"
#include "defineMemory.h"
#include "Node.h"
#include "stIONode.h"
#include "stNodeValArray.h"
#include "Grid.h"
#include "Array3D.h"
#include "MPIDatatypeInfo.h"


class  Tree {
private:
    int  rank_;
    int  num_procs_;

    int  num_nodes_global_total_;
    int  num_nodes_global_max_;
    int  num_nodes_global_min_;
    int  num_nodes_lv_global_total_[DefAMR::LV_MAX];
    int  num_nodes_lv_global_max_  [DefAMR::LV_MAX];
    int  num_nodes_lv_global_min_  [DefAMR::LV_MAX];

    int  num_nodes_lv_[DefAMR::LV_MAX];

    MemType  memType_;
    std::vector<Node>  nodes_;

    // external nodes (flags) //
    Node  root_node_;
    Node  external_leaf_node_;
    Node  external_neighbor_node_;
    Node  invalid_node_;
    // causion : external node index is 0. but nodes(0) is different from external node(0). //


    // mesh val //
    stNodeValArray  stnodeValArray_;

    // MPI //
    std::vector<MPIPutGetInfo>  mpiPutGetInfo_; // mpiPutGetInfo_.size() = nodes_.size() //

    std::vector<MPIPutInfo>  mpiPutInfo_[DefAMR::LV_MAX]; // vector.size() = number of communications //

    std::vector<MPIDatatypeInfoSrcDst>  lbmPutInfoSrcDst_[DefAMR::LV_MAX]; // lbm : N^3 * 27 * node, vector.size() = number of dst ranks //
    std::vector<MPIDatatypeInfoSrcDst>  valPutInfoSrcDst_[DefAMR::LV_MAX]; // val : N^3 * node //

    MPIPackUnpackInfo        mpiPackUnpackInfo_[DefAMR::LV_MAX]; // lbm & val //

public:
    Tree ()
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &num_procs_);
        memType_ = MemType::Managed;

        allocate_stNodeValArray(0);
        set_default_external_nodes();
    }

    ~Tree () {
        delete_stNodeValArray();
    }

public:
    const std::vector<Node>& nodes() const { return nodes_; }

          Node*  nodes(const int  id)        { return  ( id>=0 && id<nodes_.size() ) ? &nodes_[id] : &invalid_node_; }
    const Node*  nodes(const int  id)  const { return  ( id>=0 && id<nodes_.size() ) ? &nodes_[id] : &invalid_node_; }

    // dammy nodes (flags) //
          Node*  root_node()                      { return  &root_node_; }
          Node*  external_leaf_node()             { return  &external_leaf_node_; }
          Node*  external_neighbor_node()         { return  &external_neighbor_node_; }
          Node*  invalid_node()                   { return  &invalid_node_; }

    const Node*  root_node()                const { return  &root_node_; }
    const Node*  external_leaf_node()       const { return  &external_leaf_node_; }
    const Node*  external_neighbor_node()   const { return  &external_neighbor_node_; }
    const Node*  invalid_node()             const { return  &invalid_node_; }

    // stNodeValArray //
    const stNodeValArray& stnodeValArray() const { return stnodeValArray_; }

    const int*          id_parent()         const { return stnodeValArray_.id_parent_; }
    const int*          id_child()          const { return stnodeValArray_.id_child_; }
    const int*          mesh_offsets()      const { return stnodeValArray_.mesh_offsets_; }
    const int*          mesh_offsets3x3x3() const { return stnodeValArray_.mesh_offsets3x3x3_; }
    const Array3D<int>* mesh_offsets3d()    const { return stnodeValArray_.mesh_offsets3d_; }

    const int*          node_pos_x()   const { return stnodeValArray_.node_pos_x_; }
    const int*          node_pos_y()   const { return stnodeValArray_.node_pos_y_; }
    const int*          node_pos_z()   const { return stnodeValArray_.node_pos_z_; }

          std::vector<MPIPutGetInfo>& mpiPutGetInfo()       { return mpiPutGetInfo_; }
    const std::vector<MPIPutGetInfo>& mpiPutGetInfo() const { return mpiPutGetInfo_; }

          std::vector<MPIPutInfo>& mpiPutInfo(const int lv)       { return mpiPutInfo_[lv]; }
    const std::vector<MPIPutInfo>& mpiPutInfo(const int lv) const { return mpiPutInfo_[lv]; }

          std::vector<MPIDatatypeInfoSrcDst>& lbmPutInfoSrcDst(const int lv)       { return lbmPutInfoSrcDst_[lv]; }
    const std::vector<MPIDatatypeInfoSrcDst>& lbmPutInfoSrcDst(const int lv) const { return lbmPutInfoSrcDst_[lv]; }

          std::vector<MPIDatatypeInfoSrcDst>& valPutInfoSrcDst(const int lv)       { return valPutInfoSrcDst_[lv]; }
    const std::vector<MPIDatatypeInfoSrcDst>& valPutInfoSrcDst(const int lv) const { return valPutInfoSrcDst_[lv]; }

          MPIPackUnpackInfo&  mpiPackUnpackInfo(const int lv)       { return mpiPackUnpackInfo_[lv]; }
    const MPIPackUnpackInfo&  mpiPackUnpackInfo(const int lv) const { return mpiPackUnpackInfo_[lv]; }

public:
    // maximum number of nodes //
    int   number_of_nodes_lv(const int lv)  const { return  num_nodes_lv_[lv]; }
    int   number_of_nodes_lv_global_total(const int lv)  const { return  num_nodes_lv_global_total_[lv]; }
    int   number_of_nodes_lv_global_max  (const int lv)  const { return  num_nodes_lv_global_max_  [lv]; }
    int   number_of_nodes_lv_global_min  (const int lv)  const { return  num_nodes_lv_global_min_  [lv]; }

    int   number_of_nodes()        const { return  nodes_.size(); }
//    int   number_of_nodes()        const { return  num_nodes_; }
    int   number_of_nodes_global_total() const { return  num_nodes_global_total_; }
    int   number_of_nodes_global_max()   const { return  num_nodes_global_max_; }
    int   number_of_nodes_global_min()   const { return  num_nodes_global_min_; }

    void  add_node  (const int  index) { nodes_.emplace_back(index); }
    void  resize_node(const int nsize) { nodes_.resize(nsize); }

    Node*  external_node(NodeTypes::Type  node_type)
    {
        return  ( node_type == NodeTypes::RootT             ) ? root_node() :
                ( node_type == NodeTypes::ExternalLeafT     ) ? external_leaf_node() :
                ( node_type == NodeTypes::ExternalNeighborT ) ? external_neighbor_node() :
                ( node_type == NodeTypes::InvalidT          ) ? invalid_node() :
                                                                NULL;
    }

private:
    void  set_default_external_nodes()
    {
        set_default_external_node( &root_node_,              NodeTypes::RootT );
        set_default_external_node( &external_leaf_node_,     NodeTypes::ExternalLeafT );
        set_default_external_node( &external_neighbor_node_, NodeTypes::ExternalNeighborT );
        set_default_external_node( &invalid_node_,           NodeTypes::InvalidT );
    }

    void  error_msg (const std::string  str)  const { std::cout << str << std::endl; }

public:
    void  preset_tree_data(const Grid*  grids);

    // initialize //
    void  init_tree_data(const Grid*  grids);

    // IO //
    std::vector<IONode>  make_vectorIONode()  const;

    void  write_stNodeValArray(std::ofstream& fout, const std::string  filename) const;
    void  read_stNodeValArray(std::ifstream& fin, const std::string  filename);

    void  write_node_structures(std::ofstream& fout, const std::string  filename, const std::vector<IONode>  ioNodes) const;
    void  read_node_structures(std::ifstream& fin, const std::string  filename);

    void  copy_node_structures(const Tree&  other, const bool  is_reset = false);

    void  write_MPIPutGetInfo(std::ofstream& fout, const std::string  filename) const;
    void  read_MPIPutGetInfo(std::ifstream& fin, const std::string  filename);

    IONode  makeIONode(const Node&  node) const;

    void  copyMPIPutGetInfo(const std::vector<MPIPutGetInfo>& mpiPutGetInfo);
    void  copy_stNodeValArray(int n, const stNodeValArray& stnodeValArray);

    std::string  filename(const int rank, const int step) const ;

    void createMPIPutInformations();

    void  set_num_nodes_lv();

private:
    void  init_stNodeValArray();

    void  allocate_stNodeValArray(const int n);
    void  reallocate_stNodeValArray(const int n);
    void  delete_stNodeValArray();

    void  set_default_external_node(Node*  node, NodeTypes::Type  node_type);

    void  set_node_i       (const int  i, const IONode&  data);
    void  check_node_index (const int  i, const IONode&  data)  const;

    void
    createMPIPutInfoSrcDst_lv(
              std::vector<MPIPutInfo>&    mpiPutInfo,
        const int                         lv,
        const std::vector<MPIPutGetInfo>& mpiPutGetInfo
        );

    void
    createMPIDataypeInfoSrcDst_lv(
              std::vector<MPIDatatypeInfoSrcDst>&   mpiDatatypeInfoSrcDst,
        const int                                   lv,
        const std::vector<MPIPutInfo>&              mpiPutInfo,
        const int                                   nQ
        );


    void
    createMPIPackUnpack_lv(
              MPIPackUnpackInfo&                    mpiPackUnpackInfo,
        const int                                   lv,
        const std::vector<MPIDatatypeInfoSrcDst>&   valDatatypeInfoSrcDst
        );


    void
    allocate_MPIPackUnpackInfo(
              MPIPackUnpackInfo&    mpiPackUnpackInfo,
        const int                   num_leaf_send,
        const int                   num_leaf_recv
        );


    void release_MPIPackUnpackInfo(MPIPackUnpackInfo&    mpiPackUnpackInfo);

};


#endif
