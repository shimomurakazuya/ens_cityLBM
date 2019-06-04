#include "Node.h"


Array3D<int>  Node::
neighbor_mesh_offsets()
const
{
    Array3D<int>  tmp;

    for (int k=-1; k<=1; k++) {
        for (int j=-1; j<=1; j++) {
            for (int i=-1; i<=1; i++) {
//                tmp.val[k+1][j+1][i+1] = neighbor(i, j, k)->mesh_offset();
                tmp.val[k+1][j+1][i+1] = (neighbor(i, j, k)->node_type() == NodeTypes::LeafT)
                                        ? neighbor(i, j, k)->mesh_offset() :
                                          neighbor(0, 0, 0)->mesh_offset() ;
            }
        }
    }

    return  tmp;
}


Array3D<int>  Node::
neighbor_lbm_mesh_offsets()
const
{
    Array3D<int>  tmp;

    for (int k=-1; k<=1; k++) {
        for (int j=-1; j<=1; j++) {
            for (int i=-1; i<=1; i++) {
//                tmp.val[k+1][j+1][i+1] = neighbor(i, j, k)->mesh_offset() * LBM_velocity_model::nQ;
                tmp.val[k+1][j+1][i+1] = (neighbor(i, j, k)->node_type() == NodeTypes::LeafT)
                                        ? neighbor(i, j, k)->mesh_offset() * LBM_velocity_model::nQ :
                                          neighbor(0, 0, 0)->mesh_offset() * LBM_velocity_model::nQ;
            }
        }
    }

    return  tmp;
}


void  Node::
set_default_value (const int  index)
{
    set_index(index);

    // connections //
    set_parent(this);
    set_child (this);

    set_neighbor (this, -1,  0,  0);
    set_neighbor (this,  1,  0,  0);
    set_neighbor (this,  0, -1,  0);
    set_neighbor (this,  0,  1,  0);
    set_neighbor (this,  0,  0, -1);
    set_neighbor (this,  0,  0,  1);

    // NodeInformations //
    set_node_type(NodeTypes::LeafT);
    set_position(0, 0, 0, 0);

    // NodeCalFlags //
    // defalut flag is false defined in constructor //

    // MeshInfo //
}


void  Node::
set_position(const int  lv, const int  x, const int  y, const int  z)
{
    position_.lv = lv;
    position_.i  = x;
    position_.j  = y;
    position_.k  = z;
}


void  Node::
set_position(const NodeIndices  nodeIndices)
{
    set_position(nodeIndices.lv, nodeIndices.i, nodeIndices.j, nodeIndices.k);
}


void  Node::
set_NodeCalFlags(const NodeCalFlags nodeCalFlags)
{
    nodeCalFlags_.set_flag_Cal   ( nodeCalFlags.Cal() );
    nodeCalFlags_.set_flag_L2C   ( nodeCalFlags.L2C() );
    nodeCalFlags_.set_flag_C2L   ( nodeCalFlags.C2L() );
    nodeCalFlags_.set_flag_L2F   ( nodeCalFlags.L2F() );
    nodeCalFlags_.set_flag_F2L   ( nodeCalFlags.F2L() );
    nodeCalFlags_.set_flag_putMPI( nodeCalFlags.putMPI() );
    nodeCalFlags_.set_flag_getMPI( nodeCalFlags.getMPI() );
}


void  Node::
set_MeshInfo(const MeshInfo  meshInfo)
{
    meshInfo_.mesh_offset = meshInfo.mesh_offset;
}


bool  Node::
is_same_position (const Node&  node)
const
{
    return  (  level() == node.level()
            && x()     == node.x()
            && y()     == node.y()
            && z()     == node.z() );
}


bool  Node::
has_same_connections (const Node&  node)
const
{
    return  (  parent()    == node.parent()
            && child()     == node.child()
            && neighbor(0) == node.neighbor(0)
            && neighbor(1) == node.neighbor(1)
            && neighbor(2) == node.neighbor(2)
            && neighbor(3) == node.neighbor(3)
            && neighbor(4) == node.neighbor(4)
            && neighbor(5) == node.neighbor(5) );
}
