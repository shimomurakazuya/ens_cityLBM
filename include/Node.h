#pragma once
#ifndef NODE_H_
#define NODE_H_


#include <iostream>
#include <cstdlib>
#include <vector>

#include "defineLBM.h"
#include "stNodeInformations.h"
#include "NodeCalFlags.h"
#include "Array3D.h"


// oct-tree //
class  Node {
private:
    int     index_; // global index of node vectors //

    // connections //
    Node*   neighbor_[6];
    Node*   parent_;    // lower  level : lv - 1 //
    Node*   child_;     // higher level : lv + 1 //

    // NodeInformations //
    NodeTypes::Type  node_type_;
    NodeIndices      position_;

    // NodeCalFlags //
    NodeCalFlags     nodeCalFlags_;

    // MeshInfo //
    MeshInfo         meshInfo_;

public:
    Node () {}

    Node (const int  index) :
        index_(index),
        node_type_(NodeTypes::LeafT)
    {
        set_default_value (index);
    }

    ~Node () {}

public:
    int  index()        const { return  index_; }

    // connections //
          Node* node()       { return  this; }
    const Node* node() const { return  this; }

          Node* neighbor()                     { return  node(); }
    const Node* neighbor()               const { return  node(); }
    const Node* neighbor(const int  ii)  const { return  neighbor_[ii]; }

          Node* parent()        { return parent_; }
          Node* child ()        { return child_; }

    const Node* parent()  const { return parent_; }
    const Node* child ()  const { return child_; }

    // NodeInformations //
    NodeTypes::Type  node_type()  const { return node_type_; }

    int  level()        const { return  position_.lv; }
    int  depth()        const { return  position_.lv; }
    int  x()            const { return  position_.i; }
    int  y()            const { return  position_.j; }
    int  z()            const { return  position_.k; }

    const NodeIndices  position()   const { return  position_; }

    // NodeCalFlags //
          NodeCalFlags& nodeCalFlags()        { return  nodeCalFlags_; }
    const NodeCalFlags& nodeCalFlags()  const { return  nodeCalFlags_; }

    // MeshInfo //
          MeshInfo  meshInfo()        { return  meshInfo_; }
    const MeshInfo  meshInfo()  const { return  meshInfo_; }

    int  mesh_offset()      const { return  meshInfo().mesh_offset; }
    int  mesh_offset_lbm()  const { return  meshInfo().mesh_offset*LBM_velocity_model::nQ; }

    Array3D<int>  neighbor_mesh_offsets()      const;
    Array3D<int>  neighbor_lbm_mesh_offsets()  const;

public:
    // connections //
    Node*  neighbor(const int  i, const int  j, const int  k)
    {
        if      ( i ==  0 && j ==  0 && k ==  0 )   { return  node(); }

        else if ( i <= -1 ) { return neighbor_[ix_neighbor(-1,  0,  0)]->neighbor(i+1, j  , k  ); }
        else if ( i >=  1 ) { return neighbor_[ix_neighbor( 1,  0,  0)]->neighbor(i-1, j  , k  ); }
        else if ( j <= -1 ) { return neighbor_[ix_neighbor( 0, -1,  0)]->neighbor(i,   j+1, k  ); }
        else if ( j >=  1 ) { return neighbor_[ix_neighbor( 0,  1,  0)]->neighbor(i,   j-1, k  ); }
        else if ( k <= -1 ) { return neighbor_[ix_neighbor( 0,  0, -1)]->neighbor(i,   j  , k+1); }
        else if ( k >=  1 ) { return neighbor_[ix_neighbor( 0,  0,  1)]->neighbor(i,   j  , k-1); }

        else { error_msg(__PRETTY_FUNCTION__); exit(-1); }
    }

    const Node*  neighbor(const int  i, const int  j, const int  k)  const
    {
        if      ( i ==  0 && j ==  0 && k ==  0 )   { return  node(); }

        else if ( i <= -1 ) { return neighbor_[ix_neighbor(-1,  0,  0)]->neighbor(i+1, j  , k  ); }
        else if ( i >=  1 ) { return neighbor_[ix_neighbor( 1,  0,  0)]->neighbor(i-1, j  , k  ); }
        else if ( j <= -1 ) { return neighbor_[ix_neighbor( 0, -1,  0)]->neighbor(i,   j+1, k  ); }
        else if ( j >=  1 ) { return neighbor_[ix_neighbor( 0,  1,  0)]->neighbor(i,   j-1, k  ); }
        else if ( k <= -1 ) { return neighbor_[ix_neighbor( 0,  0, -1)]->neighbor(i,   j  , k+1); }
        else if ( k >=  1 ) { return neighbor_[ix_neighbor( 0,  0,  1)]->neighbor(i,   j  , k-1); }

        else { error_msg(__PRETTY_FUNCTION__); exit(-1); }
    }


    Node*  depth_search(const int  lv)
    {
        if      ( lv ==  0 ) { return  node(); }
        else if ( lv == -1 ) { return  parent(); }
        else if ( lv ==  1 ) { return  child(); }

        else if ( lv <= -2 ) { return  parent()->depth_search( lv+1 ); }
        else if ( lv >=  2 ) { return  child ()->depth_search( lv-1 ); }

        else { error_msg(__PRETTY_FUNCTION__); exit(-1); }
    }

    const Node*  depth_search(const int  lv)  const
    {
        if      ( lv ==  0 ) { return  node(); }
        else if ( lv == -1 ) { return  parent(); }
        else if ( lv ==  1 ) { return  child(); }

        else if ( lv <= -2 ) { return  parent()->depth_search( lv+1 ); }
        else if ( lv >=  2 ) { return  child ()->depth_search( lv-1 ); }

        else { error_msg(__PRETTY_FUNCTION__); exit(-1); }
    }

          Node* node(const int  lv, const int  i, const int  j, const int  k)        { return  depth_search(lv)->neighbor(i, j, k); }
    const Node* node(const int  lv, const int  i, const int  j, const int  k)  const { return  depth_search(lv)->neighbor(i, j, k); }

          Node* node(const NodeIndices& nodeIndices)        { return  node( nodeIndices.lv, nodeIndices.i, nodeIndices.j, nodeIndices.k ); }
    const Node* node(const NodeIndices& nodeIndices)  const { return  node( nodeIndices.lv, nodeIndices.i, nodeIndices.j, nodeIndices.k ); }


    // recursion function //
    const Node* node_search_impl()  const {  return  this; }
    const Node* node_search_impl(const NodeIndices&  nodeIndices)  const {  return  node( nodeIndices ); }

    template <class... Rest>
    const Node* node_search_impl(const NodeIndices& nodeIndices, const Rest&... rest)
    const
    {
        return  node(nodeIndices.lv, nodeIndices.i, nodeIndices.j, nodeIndices.k)->node_search_impl(rest...);
    }

    template <class... Ts>
    const Node* node_search(const Ts&... args) const { return  node_search_impl( args... ); }

public:
    void  set_neighbor(Node* node, const int  i, const int  j, const int  k)  { neighbor_[ix_neighbor(i, j, k)] = node; }
    void  set_neighbor(Node* node, const int  ix)                             { neighbor_[ix]                   = node; }

    void  set_parent (Node* node) { parent_ = node; }
    void  set_child  (Node* node) { child_  = node; }

    void  set_index(const int  index) { index_ = index; }

    void  set_node_type(const NodeTypes::Type  node_type)  { node_type_ = node_type; }

    void  set_level(const int  lv) { position_.lv = lv; }
    void  set_depth(const int  lv) { position_.lv = lv; }

public:
    void  set_default_value (const int  index);
    void  set_value_from_other (const int  index, const std::vector<Node>  other);

    void  set_position(const int  lv, const int  x, const int  y, const int  z);
    void  set_position(const NodeIndices  nodeIndices);

    void  set_mesh_offset(const int  mesh_offset)  { meshInfo_.mesh_offset = mesh_offset; }

    void  set_NodeCalFlags(const NodeCalFlags nodeCalFlags);
    void  set_MeshInfo(const MeshInfo  meshInfo);

    bool  is_same_position     (const Node&  node)  const;
    bool  has_same_connections (const Node&  node)  const;

private:
    int  ix_neighbor(const int  i, const int  j, const int  k)
    const
    {
        if      (i == -1 && j ==  0 && k ==  0) {  return  0; }
        else if (i ==  1 && j ==  0 && k ==  0) {  return  1; }
        else if (i ==  0 && j == -1 && k ==  0) {  return  2; }
        else if (i ==  0 && j ==  1 && k ==  0) {  return  3; }
        else if (i ==  0 && j ==  0 && k == -1) {  return  4; }
        else if (i ==  0 && j ==  0 && k ==  1) {  return  5; }
        else { error_msg(__PRETTY_FUNCTION__); exit(-1); }
    }

private:
    void  error_msg (const std::string  str)  const { std::cout << str << std::endl; }

};


#endif
