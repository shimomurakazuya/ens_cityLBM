#include "Tree.h"
#include "InitTree.h"
#include "InitTree_for_MPI.h"
#include "defineFilenames.h"
#include "defineLBM.h"
#include "FuncAllocate.h"
#include "Index.h"

#include <algorithm>


// public //
void  Tree::
preset_tree_data(const Grid*  grids)
{
    InitTree_for_MPI  initTree_for_MPI;
    initTree_for_MPI.init_tree_uniform3d_div(*this, grids, DefAMR::LV_MAX);

    // offset //
    init_stNodeValArray();

    set_num_nodes_lv();

    // put get info //
//    std::cout << __PRETTY_FUNCTION__ << "--------------------------------------" << std::endl;
    createMPIPutInformations();
}


void  Tree::
init_tree_data(const Grid*  grids)
{
//    InitTree::init_tree_uniform3d(*this, grids, DefAMR::LV_MAX);
    InitTree::init_tree_amr(*this, grids, DefAMR::LV_MAX);

    // offset //
    init_stNodeValArray();

    set_num_nodes_lv();
}


// IO //
std::vector<IONode>  Tree::
make_vectorIONode()
const
{
//    std::cout << "nodes_.size() = " << nodes_.size() << std::endl;
    std::vector<IONode>  ioNodes;
    for (int i=0; i<(int)nodes_.size(); i++) { ioNodes.emplace_back(makeIONode(nodes_[i])); }

    return  ioNodes;
}


void  Tree::
write_stNodeValArray(std::ofstream& fout, const std::string  filename)
const
{
    constexpr int  nQ = LBM_velocity_model::nQ;
    const int  num_nodes = nodes_.size();

    fout.write(( char * ) &num_nodes, sizeof( int ));

    fout.write(( char * ) stnodeValArray_.id_parent_,      sizeof( int          )*num_nodes );
    fout.write(( char * ) stnodeValArray_.id_child_,       sizeof( int          )*num_nodes );

    fout.write(( char * ) stnodeValArray_.mesh_offsets_,      sizeof( int          )*num_nodes );
    fout.write(( char * ) stnodeValArray_.mesh_offsets3d_,    sizeof( Array3D<int> )*num_nodes );
    fout.write(( char * ) stnodeValArray_.mesh_offsets3x3x3_, sizeof( int          )*num_nodes*nQ );

    fout.write(( char * ) stnodeValArray_.node_pos_x_,     sizeof( int          )*num_nodes );
    fout.write(( char * ) stnodeValArray_.node_pos_y_,     sizeof( int          )*num_nodes );
    fout.write(( char * ) stnodeValArray_.node_pos_z_,     sizeof( int          )*num_nodes );
}


void  Tree::
read_stNodeValArray(std::ifstream& fin, const std::string  filename)
{
    constexpr int  nQ = LBM_velocity_model::nQ;
    int  num_nodes;

    fin.read(( char * ) &num_nodes, sizeof( int ));
    reallocate_stNodeValArray(num_nodes);

    fin.read(( char * ) stnodeValArray_.id_parent_,      sizeof( int          )*num_nodes );
    fin.read(( char * ) stnodeValArray_.id_child_,       sizeof( int          )*num_nodes );

    fin.read(( char * ) stnodeValArray_.mesh_offsets_,      sizeof( int          )*num_nodes );
    fin.read(( char * ) stnodeValArray_.mesh_offsets3d_,    sizeof( Array3D<int> )*num_nodes );
    fin.read(( char * ) stnodeValArray_.mesh_offsets3x3x3_, sizeof( int          )*num_nodes*nQ );

    fin.read(( char * ) stnodeValArray_.node_pos_x_,     sizeof( int          )*num_nodes );
    fin.read(( char * ) stnodeValArray_.node_pos_y_,     sizeof( int          )*num_nodes );
    fin.read(( char * ) stnodeValArray_.node_pos_z_,     sizeof( int          )*num_nodes );
}


void  Tree::
write_node_structures(std::ofstream& fout, const std::string  filename, const std::vector<IONode>  ioNodes)
const
{
    fout.write(( char * ) &DefAMR::LV_MAX, sizeof( int ));
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        fout.write(( char * ) &num_nodes_lv_[i], sizeof( int ));
    }

    const int  nodes_size = ioNodes.size();
    fout.write(( char * ) &nodes_size, sizeof( int ));

    for (int i=0; i<nodes_size; i++) {
        const IONode  data = ioNodes[i];

        fout.write(( char * ) &data, sizeof( data ));
    }
}


void  Tree::
read_node_structures(std::ifstream& fin, const std::string  filename)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    int  lv_max;
    fin.read(( char * ) &lv_max, sizeof( int ));
    if (lv_max != DefAMR::LV_MAX) { std::cout << __PRETTY_FUNCTION__ << " : " << "error lv_max" << std::endl; }

    for (int i=0; i<DefAMR::LV_MAX; i++) {
        fin.read(( char * ) &num_nodes_lv_[i], sizeof( int ));
//        std::cout << "num_nodes_lv_[" << i << "] = " << num_nodes_lv_[i] << std::endl;
    }


    int  nodes_size;
    fin.read(( char * ) &nodes_size, sizeof( int ));
    for (int i=0; i<nodes_size; i++) { nodes_.push_back(i); }


    for (int i=0; i<nodes_size; i++) {
        IONode  data;

        fin.read(( char * ) &data, sizeof( data ));

        set_node_i      (i, data);
        check_node_index(i, data);
    }

    set_num_nodes_lv();
}


void  Tree::
copy_node_structures(const Tree&  other, const bool  is_reset)
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    std::copy( other.num_nodes_lv_, other.num_nodes_lv_ + DefAMR::LV_MAX, num_nodes_lv_ );

    if (is_reset) {
        nodes_.resize( other.number_of_nodes() );
        nodes_.shrink_to_fit();
    }

    for (int i=0; i<other.number_of_nodes(); i++) {

        nodes(i)->set_index( other.nodes(i)->index() );

        // connections //
        if ( other.nodes(i)->parent()->node_type() == NodeTypes::LeafT ) {
            nodes(i)->set_parent( nodes( other.nodes(i)->parent()->index() ) );
        }
        else {
            nodes(i)->set_parent( external_node( other.nodes(i)->parent()->node_type() ) );
        }

        if ( other.nodes(i)->child()->node_type() == NodeTypes::LeafT ) {
            nodes(i)->set_child ( nodes( other.nodes(i)->child ()->index() ) );
        }
        else {
            nodes(i)->set_child( external_node( other.nodes(i)->child()->node_type() ) );
        }

        for (int ii=0; ii<6; ii++) {
            if ( other.nodes(i)->neighbor(ii)->node_type() == NodeTypes::LeafT ) {
                nodes(i)->set_neighbor ( nodes( other.nodes(i)->neighbor(ii)->index() ), ii );
            }
            else {
                nodes(i)->set_neighbor( external_node( other.nodes(i)->neighbor(ii)->node_type() ), ii );
            }
        }

        // NodeInformations //
        nodes(i)->set_node_type( other.nodes(i)->node_type() );
        nodes(i)->set_position ( other.nodes(i)->position() );

        // NodeCalFlags //
        nodes(i)->set_NodeCalFlags( other.nodes(i)->nodeCalFlags() );

        // MeshInfo //
        nodes(i)->set_MeshInfo( other.nodes(i)->meshInfo() );
    }

    set_num_nodes_lv();
}


void  Tree::
write_MPIPutGetInfo(std::ofstream& fout, const std::string  filename)
const
{
    const int  num_nodes = mpiPutGetInfo_.size();
    fout.write(( char * ) &num_nodes,  sizeof( int ));

    for (int i=0; i<num_nodes; i++) {
        fout.write(( char * ) &mpiPutGetInfo_[i].flag_put,    sizeof( bool ));
        fout.write(( char * ) &mpiPutGetInfo_[i].num_vec_put, sizeof( int ));
        for (int j=0; j<mpiPutGetInfo_[i].num_vec_put; j++) {
            fout.write(( char * ) &mpiPutGetInfo_[i].rank_put[j],   sizeof( int ));
            fout.write(( char * ) &mpiPutGetInfo_[i].id_put[j],     sizeof( int ));
            fout.write(( char * ) &mpiPutGetInfo_[i].offset_put[j], sizeof( int ));
        }

        fout.write(( char * ) &mpiPutGetInfo_[i].flag_get,   sizeof( bool ));
        fout.write(( char * ) &mpiPutGetInfo_[i].rank_get,   sizeof( int ));
        fout.write(( char * ) &mpiPutGetInfo_[i].id_get,     sizeof( int ));
        fout.write(( char * ) &mpiPutGetInfo_[i].offset_get, sizeof( int ));
    }
}


void  Tree::
read_MPIPutGetInfo(std::ifstream& fin, const std::string  filename)
{
    int  num_nodes;
    fin.read(( char * ) &num_nodes,  sizeof( int ));
    mpiPutGetInfo_.resize(num_nodes);
    mpiPutGetInfo_.shrink_to_fit();

    for (int i=0; i<num_nodes; i++) {
        fin.read(( char * ) &mpiPutGetInfo_[i].flag_put,    sizeof( bool ));

        fin.read(( char * ) &mpiPutGetInfo_[i].num_vec_put, sizeof( int ));
        mpiPutGetInfo_[i].rank_put  .resize( mpiPutGetInfo_[i].num_vec_put );
        mpiPutGetInfo_[i].id_put    .resize( mpiPutGetInfo_[i].num_vec_put );
        mpiPutGetInfo_[i].offset_put.resize( mpiPutGetInfo_[i].num_vec_put );
        for (int j=0; j<mpiPutGetInfo_[i].num_vec_put; j++) {
            fin.read(( char * ) &mpiPutGetInfo_[i].rank_put[j],   sizeof( int ));
            fin.read(( char * ) &mpiPutGetInfo_[i].id_put[j],     sizeof( int ));
            fin.read(( char * ) &mpiPutGetInfo_[i].offset_put[j], sizeof( int ));
        }

        fin.read(( char * ) &mpiPutGetInfo_[i].flag_get,   sizeof( bool ));
        fin.read(( char * ) &mpiPutGetInfo_[i].rank_get,   sizeof( int ));
        fin.read(( char * ) &mpiPutGetInfo_[i].id_get,     sizeof( int ));
        fin.read(( char * ) &mpiPutGetInfo_[i].offset_get, sizeof( int ));
    }

    set_num_nodes_lv();
}


IONode  Tree::
makeIONode(const Node&  node)
const
{
    // convert neighbor address to index //
    IONode  ioNode  {
                node.index(),
                // connections //
                node.neighbor(0)->index(), node.neighbor(1)->index(),
                node.neighbor(2)->index(), node.neighbor(3)->index(),
                node.neighbor(4)->index(), node.neighbor(5)->index(),
                node.parent()->index(),
                node.child()->index(),
                // NodeInformations //
                node.node_type(),
                node.position(),
                // NodeCalFlags //
                node.nodeCalFlags().Cal(),
                node.nodeCalFlags().L2C(),
                node.nodeCalFlags().C2L(),
                node.nodeCalFlags().L2F(),
                node.nodeCalFlags().F2L(),
                node.nodeCalFlags().putMPI(),
                node.nodeCalFlags().getMPI(),
                // MeshInfo //
                node.meshInfo(),
                };

    return  ioNode;
}


void Tree::
copyMPIPutGetInfo(const std::vector<MPIPutGetInfo>& mpiPutGetInfo)
{
    const int nsize = mpiPutGetInfo.size();
    mpiPutGetInfo_.resize(nsize);
    mpiPutGetInfo_.shrink_to_fit();

    for (int i=0; i<nsize; i++) {
        mpiPutGetInfo_[i].flag_put    = mpiPutGetInfo[i].flag_put;
        mpiPutGetInfo_[i].num_vec_put = mpiPutGetInfo[i].num_vec_put;
        mpiPutGetInfo_[i].rank_put    = mpiPutGetInfo[i].rank_put;
        mpiPutGetInfo_[i].id_put      = mpiPutGetInfo[i].id_put;
        mpiPutGetInfo_[i].offset_put  = mpiPutGetInfo[i].offset_put;

        mpiPutGetInfo_[i].flag_get    = mpiPutGetInfo[i].flag_get;
        mpiPutGetInfo_[i].rank_get    = mpiPutGetInfo[i].rank_get;
        mpiPutGetInfo_[i].id_get      = mpiPutGetInfo[i].id_get;
        mpiPutGetInfo_[i].offset_get  = mpiPutGetInfo[i].offset_get;


//        std::cout << "flag_put, flag_get = " << mpiPutGetInfo_[i].flag_put << ", " << mpiPutGetInfo_[i].flag_get << std::endl;
    }
}


void  Tree::
copy_stNodeValArray(const int n, const stNodeValArray& stnodeValArray)
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    reallocate_stNodeValArray(n);

    for (int i=0; i<n; i++) {
        stnodeValArray_.id_parent_     [i] = stnodeValArray.id_parent_     [i];
        stnodeValArray_.id_child_      [i] = stnodeValArray.id_child_      [i];

        stnodeValArray_.mesh_offsets_     [i] = stnodeValArray.mesh_offsets_  [i];
        stnodeValArray_.mesh_offsets3d_   [i] = stnodeValArray.mesh_offsets3d_[i];
        for (int idv=0; idv<nQ; idv++) {
            stnodeValArray_.mesh_offsets3x3x3_[i*nQ + idv] = stnodeValArray.mesh_offsets3x3x3_[i*nQ + idv];
        }

        stnodeValArray_.node_pos_x_    [i] = stnodeValArray.node_pos_x_    [i];
        stnodeValArray_.node_pos_y_    [i] = stnodeValArray.node_pos_y_    [i];
        stnodeValArray_.node_pos_z_    [i] = stnodeValArray.node_pos_z_    [i];
    }
}


std::string  Tree::
filename(const int rank, const int step)
const
{
    return    Foldernames::io_folder + "/"
            + Filenames::tree_name0
            + "-rank" + std::to_string(rank)
            + "-step" + std::to_string(step)
            + ".dat";
}


// private //
void  Tree::
init_stNodeValArray()
{
//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    constexpr int  nQ = LBM_velocity_model::nQ;
    const int n = nodes_.size();

    reallocate_stNodeValArray(n);

    for (int i=0; i<n; i++) {
        if ( i != nodes(i)->index() ) { std::cout << __PRETTY_FUNCTION__ << " : " << "i != nodes(i)->index()" << std::endl; }

        const int          id_parent = (nodes(i)->node_type() == NodeTypes::LeafT) ? nodes(i)->parent()->index() : -1;
        const int          id_child  = (nodes(i)->node_type() == NodeTypes::LeafT) ? nodes(i)->child ()->index() : -1;

        const Array3D<int> offset3d  = nodes(i)->neighbor_mesh_offsets();
        const int          offset    = offset3d.offset0();

        const int          node_pos_x = nodes(i)->x();
        const int          node_pos_y = nodes(i)->y();
        const int          node_pos_z = nodes(i)->z();

        stnodeValArray_.id_parent_     [i] = id_parent;
        stnodeValArray_.id_child_      [i] = id_child;

        stnodeValArray_.mesh_offsets_  [i] = offset;
        stnodeValArray_.mesh_offsets3d_[i] = offset3d;

        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv = Index::idv(iv, jv, kv);

                    stnodeValArray_.mesh_offsets3x3x3_[i*nQ + idv] = offset3d.offset(iv,jv,kv);
                }
            }
        }

        stnodeValArray_.node_pos_x_    [i] = node_pos_x;
        stnodeValArray_.node_pos_y_    [i] = node_pos_y;
        stnodeValArray_.node_pos_z_    [i] = node_pos_z;
    }
}


void  Tree::
allocate_stNodeValArray(const int n)
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    FuncAllocate::allocate_value<int>( &stnodeValArray_.id_parent_, n, memType_ );
    FuncAllocate::allocate_value<int>( &stnodeValArray_.id_child_,  n, memType_ );

    FuncAllocate::allocate_value<int>( &stnodeValArray_.mesh_offsets_,       n,    memType_ );
    stnodeValArray_.mesh_offsets3d_    = new Array3D<int>[n];
    FuncAllocate::allocate_value<int>( &stnodeValArray_.mesh_offsets3x3x3_,  n*nQ, memType_ );

    FuncAllocate::allocate_value<int>( &stnodeValArray_.node_pos_x_,  n, memType_ );
    FuncAllocate::allocate_value<int>( &stnodeValArray_.node_pos_y_,  n, memType_ );
    FuncAllocate::allocate_value<int>( &stnodeValArray_.node_pos_z_,  n, memType_ );

//    stnodeValArray_.id_parent_      = new int[n];
//    stnodeValArray_.id_child_       = new int[n];
//
//    stnodeValArray_.mesh_offsets_      = new int[n];
//    stnodeValArray_.mesh_offsets3d_    = new Array3D<int>[n];
//    stnodeValArray_.mesh_offsets3x3x3_ = new int[n*nQ];
//
//    stnodeValArray_.node_pos_x_     = new int[n];
//    stnodeValArray_.node_pos_y_     = new int[n];
//    stnodeValArray_.node_pos_z_     = new int[n];
}


void  Tree::
reallocate_stNodeValArray(const int n)
{
    delete_stNodeValArray();
    allocate_stNodeValArray(n);
}


void  Tree::
delete_stNodeValArray()
{
    FuncAllocate::release_value<int>(stnodeValArray_.id_parent_, memType_);
    FuncAllocate::release_value<int>(stnodeValArray_.id_child_,  memType_);

    FuncAllocate::release_value<int>(stnodeValArray_.mesh_offsets_,       memType_);
    delete [] stnodeValArray_.mesh_offsets3d_;
    FuncAllocate::release_value<int>(stnodeValArray_.mesh_offsets3x3x3_,  memType_);

    FuncAllocate::release_value<int>(stnodeValArray_.node_pos_x_,  memType_);
    FuncAllocate::release_value<int>(stnodeValArray_.node_pos_y_,  memType_);
    FuncAllocate::release_value<int>(stnodeValArray_.node_pos_z_,  memType_);

//    delete [] stnodeValArray_.id_parent_;
//    delete [] stnodeValArray_.id_child_;
//
//    delete [] stnodeValArray_.mesh_offsets_;
//    delete [] stnodeValArray_.mesh_offsets3d_;
//    delete [] stnodeValArray_.mesh_offsets3x3x3_;
//
//    delete [] stnodeValArray_.node_pos_x_;
//    delete [] stnodeValArray_.node_pos_y_;
//    delete [] stnodeValArray_.node_pos_z_;
}


void  Tree::
set_default_external_node(Node*  node, const NodeTypes::Type  node_type)
{
    const int  index0 = 0;
    const int  depth  = 0;

    node->set_index(index0);
    node->set_depth(depth);
    node->set_node_type(node_type);

    node->set_parent(&root_node_);
    node->set_child (&root_node_);

    node->set_neighbor (&root_node_, -1,  0,  0);
    node->set_neighbor (&root_node_,  1,  0,  0);
    node->set_neighbor (&root_node_,  0, -1,  0);
    node->set_neighbor (&root_node_,  0,  1,  0);
    node->set_neighbor (&root_node_,  0,  0, -1);
    node->set_neighbor (&root_node_,  0,  0,  1);
}


void  Tree::
set_num_nodes_lv()
{
    int  num_nodes = 0;
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        auto  is_lv = [lv](const Node& node) { return  ( node.level() == lv ); };

        num_nodes_lv_[lv] = std::count_if( nodes_.begin(), nodes_.end(), is_lv );

        num_nodes += num_nodes_lv_[lv];
    }
    MPI_Barrier(MPI_COMM_WORLD);

//    MPI_Allreduce(&num_nodes, &num_nodes_global_total_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&num_nodes, &num_nodes_global_max_  , 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//    MPI_Allreduce(&num_nodes, &num_nodes_global_min_  , 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
//        MPI_Allreduce(&num_nodes_lv_[lv], &num_nodes_lv_global_total_[lv],   1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce(&num_nodes_lv_[lv], &num_nodes_lv_global_max_  [lv],   1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//        MPI_Allreduce(&num_nodes_lv_[lv], &num_nodes_lv_global_min_  [lv],   1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
//    }

    MPI_Request requ[3*(1+DefAMR::LV_MAX)];

    int i_iallreduces = 0;
    MPI_Iallreduce(&num_nodes, &num_nodes_global_total_, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &requ[i_iallreduces]);  i_iallreduces++;
    MPI_Iallreduce(&num_nodes, &num_nodes_global_max_  , 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &requ[i_iallreduces]);  i_iallreduces++;
    MPI_Iallreduce(&num_nodes, &num_nodes_global_min_  , 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &requ[i_iallreduces]);  i_iallreduces++;
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        MPI_Iallreduce(&num_nodes_lv_[lv], &num_nodes_lv_global_total_[lv],   1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &requ[i_iallreduces]);  i_iallreduces++;
        MPI_Iallreduce(&num_nodes_lv_[lv], &num_nodes_lv_global_max_  [lv],   1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &requ[i_iallreduces]);  i_iallreduces++;
        MPI_Iallreduce(&num_nodes_lv_[lv], &num_nodes_lv_global_min_  [lv],   1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, &requ[i_iallreduces]);  i_iallreduces++;
    }
    MPI_Waitall(i_iallreduces, requ, MPI_STATUSES_IGNORE);
}


void  Tree::
set_node_i (const int  i, const IONode&  data)
{
    nodes_[i].set_index(data.index);

    // connections //
    for (int ix=0; ix<6; ix++) { nodes_[i].set_neighbor(&nodes_[data.neighbor[ix]], ix); }
    nodes_[i].set_parent(&nodes_[data.parent]);
    nodes_[i].set_child (&nodes_[data.child ]);

    // NodeInformations //
    nodes_[i].set_node_type(data.node_type);
    nodes_[i].set_position(data.position);

    // NodeCalFlags //
    nodes_[i].nodeCalFlags().set_flag_Cal   ( data.nodeCalFlags[0] );
    nodes_[i].nodeCalFlags().set_flag_L2C   ( data.nodeCalFlags[1] );
    nodes_[i].nodeCalFlags().set_flag_C2L   ( data.nodeCalFlags[2] );
    nodes_[i].nodeCalFlags().set_flag_L2F   ( data.nodeCalFlags[3] );
    nodes_[i].nodeCalFlags().set_flag_F2L   ( data.nodeCalFlags[4] );
    nodes_[i].nodeCalFlags().set_flag_putMPI( data.nodeCalFlags[5] );
    nodes_[i].nodeCalFlags().set_flag_getMPI( data.nodeCalFlags[6] );

    // MeshInfo //
    nodes_[i].set_MeshInfo(data.meshInfo);
}


void  Tree::
check_node_index (const int  i, const IONode&  data)
const
{
    if (i != data.index) {
        std::cout << "----------" << std::endl;
        std::cout << "error" << std::endl;
        std::cout << __PRETTY_FUNCTION__ << " : " << i << ", " << data.index << std::endl;
        std::cout << "----------" << std::endl;
        exit(-1);
    }
}


void Tree::
createMPIPutInformations()
{
    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        if (rank == 0) { std::cout << "lv = " << lv << " : createMPIPutInfoSrcDst_lv" << std::endl; }
        createMPIPutInfoSrcDst_lv(
            mpiPutInfo_[lv],
            lv,
            mpiPutGetInfo_
            );
    }
    MPI_Barrier(MPI_COMM_WORLD);


    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        // lbm //
        if (rank == 0) { std::cout << "lv = " << lv << " : createMPIDataypeInfoSrcDst_lv (lbm)" << std::endl; }
        createMPIDataypeInfoSrcDst_lv(
            lbmPutInfoSrcDst_[lv],
            lv,
            mpiPutInfo_[lv],
            LBM_velocity_model::nQ
            );

        // val //
        if (rank == 0) { std::cout << "lv = " << lv << " : createMPIDataypeInfoSrcDst_lv (val)" << std::endl; }
        createMPIDataypeInfoSrcDst_lv(
            valPutInfoSrcDst_[lv],
            lv,
            mpiPutInfo_[lv],
            1
            );

        // pack unpack //
        if (rank == 0) { std::cout << "lv = " << lv << " : createMPIPackUnpack_lv" << std::endl; }
        createMPIPackUnpack_lv(
            mpiPackUnpackInfo_[lv],
            lv,
            valPutInfoSrcDst_[lv]
            );

    }
    MPI_Barrier(MPI_COMM_WORLD);
}


void  Tree::
createMPIPutInfoSrcDst_lv(
          std::vector<MPIPutInfo>&    mpiPutInfo,
    const int                         lv,
    const std::vector<MPIPutGetInfo>& mpiPutGetInfo
    )
{
    const int nsize = mpiPutGetInfo.size();

    // mpiPutInfo //
    mpiPutInfo.resize(0);
    mpiPutInfo.shrink_to_fit();
    for (int i=0; i<nsize; i++) {
        if ( nodes(i)->level() != lv  ||  !mpiPutGetInfo[i].flag_put ) { continue; }

        for (int j=0; j<mpiPutGetInfo[i].num_vec_put; j++) {
            if (mpiPutGetInfo[i].id_put[j] == -1) { std::cout << __PRETTY_FUNCTION__ << " : error id_put is -1" << std::endl; }

            mpiPutInfo.push_back( MPIPutInfo{ i, mpiPutGetInfo[i].rank_put[j], mpiPutGetInfo[i].id_put[j], mpiPutGetInfo[i].offset_put[j] } );
        }
    }
    std::sort( mpiPutInfo.begin(), mpiPutInfo.end(), [](MPIPutInfo a, MPIPutInfo b) { return a.rank_put < b.rank_put; } );

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << __PRETTY_FUNCTION__ << " : rank, lv, num_mpiPutInfo = " << rank << ", " << lv << ", " << mpiPutInfo.size() << std::endl;

    mpiPutInfo.shrink_to_fit();
}


void  Tree::
createMPIDataypeInfoSrcDst_lv(
          std::vector<MPIDatatypeInfoSrcDst>&   mpiDatatypeInfoSrcDst,
    const int                                   lv,
    const std::vector<MPIPutInfo>&              mpiPutInfo,
    const int                                   nq
    )
{
    std::vector<MPIPutInfo>  tmp_mpiPutInfo;
    for (int i=0; i<(int)mpiPutInfo.size(); i++) {
        tmp_mpiPutInfo.emplace_back( mpiPutInfo[i] );
    }
    std::sort( tmp_mpiPutInfo.begin(), tmp_mpiPutInfo.end(), [](MPIPutInfo a, MPIPutInfo b) { return a.rank_put < b.rank_put; } );
    tmp_mpiPutInfo.shrink_to_fit();


    // find //
    std::vector<MPIPutInfo>::iterator it_begin  = tmp_mpiPutInfo.begin();
    std::vector<MPIPutInfo>::iterator it_end    = tmp_mpiPutInfo.end();
    std::vector<MPIPutInfo>::iterator it_find;

    while (it_begin != it_end) {
        const int target = it_begin->rank_put;
        it_find = std::find_if(it_begin, it_end, [target](MPIPutInfo& a) { return (a.rank_put != target); } ); // find : not target //
        std::sort( it_begin, it_find, [](MPIPutInfo  a, MPIPutInfo  b) { return a.id < b.id; } ); // sort : optional //

        std::vector<int>  id_src;
        std::vector<int>  id_dst;
        std::vector<MPIPutInfo>::iterator it_tmp = it_begin;
        while (it_tmp != it_find) {
            id_src.push_back( nodes(it_tmp->id)->mesh_offset() * nq );
            id_dst.push_back( it_tmp->offset_put * nq );

            it_tmp++;
        }

        const int rank_src = rank_;
        const int rank_dst = target;

        const int count   = id_src.size();
        const int blength = DefAMR::NN_LEAF * nq;

        // src & dst //
        MPIDatatypeInfo  mpiDatatypeInfo_src;
        MPIDatatypeInfo  mpiDatatypeInfo_dst;

        // src //
        mpiDatatypeInfo_src.rank    = rank_src;
        mpiDatatypeInfo_src.count   = count;
        mpiDatatypeInfo_src.blength = blength;
        mpiDatatypeInfo_src.dtype   = MFLOAT;
        MPI_Type_create_indexed_block( mpiDatatypeInfo_src.count, mpiDatatypeInfo_src.blength, id_src.data(), mpiDatatypeInfo_src.dtype, &mpiDatatypeInfo_src.mpi_datatype );

        // dst //
        mpiDatatypeInfo_dst.rank    = rank_dst;
        mpiDatatypeInfo_dst.count   = count;
        mpiDatatypeInfo_dst.blength = blength;
        mpiDatatypeInfo_dst.dtype   = MFLOAT;
        MPI_Type_create_indexed_block( mpiDatatypeInfo_dst.count, mpiDatatypeInfo_dst.blength, id_dst.data(), mpiDatatypeInfo_dst.dtype, &mpiDatatypeInfo_dst.mpi_datatype );

        // commit //
        MPI_Type_commit( &mpiDatatypeInfo_src.mpi_datatype );
        MPI_Type_commit( &mpiDatatypeInfo_dst.mpi_datatype );

        // debug //
        for (int i=0; i<(int)id_src.size(); i++) {
            mpiDatatypeInfo_src.id.push_back( id_src[i] );
            mpiDatatypeInfo_dst.id.push_back( id_dst[i] );
        }

        // registration //
        mpiDatatypeInfoSrcDst.emplace_back( MPIDatatypeInfoSrcDst{  mpiDatatypeInfo_src, mpiDatatypeInfo_dst } );

#if 0 // check //
        std::cout << "rank_src, dst = " <<  rank_src << ", " << rank_dst << " : count = " << count << std::endl;
#endif

        // iterator //
        it_begin = it_find;
    }

//    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    std::cout << __PRETTY_FUNCTION__ << " : rank, lv, num_mpiDatatypeInfoSrcDst = " << rank << ", " << lv << ", " << mpiDatatypeInfoSrcDst.size() << std::endl;
}


void Tree::
createMPIPackUnpack_lv(
          MPIPackUnpackInfo&                    mpiPackUnpackInfo,
    const int                                   lv,
    const std::vector<MPIDatatypeInfoSrcDst>&   valDatatypeInfoSrcDst
    )
{
    std::cout << __PRETTY_FUNCTION__ << " : lv = " << lv << std::endl;
    const int  rank      = rank_;
    const int  num_procs = num_procs_;

    /// info send ///
    const int  nprocs_send = valDatatypeInfoSrcDst.size();

    mpiPackUnpackInfo.sendPackUnpackCommInfo.resize( nprocs_send );

    int _nleaf_max = 0;
    int sum_leaf = 0;
    for (int i=0; i<nprocs_send; i++) {
        mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank     = valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.rank;
        mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag      = rank;
        mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf    = valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.count;
        mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf = sum_leaf;

        mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.resize(mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf);
        for (int j=0; j<mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf; j++) {
//            mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.push_back( valDatatypeInfoSrcDst[i].mpiDatatypeInfo_src.id[j] / DefAMR::NN_LEAF );
            mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf[j] = valDatatypeInfoSrcDst[i].mpiDatatypeInfo_src.id[j] / DefAMR::NN_LEAF;
        }

        sum_leaf += mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf;
    
        _nleaf_max = std::max(_nleaf_max, mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf);
    }
    int nleaf_global_max;
    MPI_Allreduce(&_nleaf_max, &nleaf_global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//    std::cout << "nleaf_global_max = " << nleaf_global_max << std::endl;
//    exit(-1);


    /// communication ///
    constexpr int nelem = 5;
    auto id_tmp_pack = [nelem](int i){ return nelem*i; };

    int  mpiPackInfo_send[num_procs*nelem];
    int  mpiPackInfo_recv[num_procs*nelem];
    for (int i=0; i<num_procs; i++) {
        mpiPackInfo_send[id_tmp_pack(i)+0] =  0;     mpiPackInfo_recv[id_tmp_pack(i)+0] =  0; // flag  //
        mpiPackInfo_send[id_tmp_pack(i)+1] = -1;     mpiPackInfo_recv[id_tmp_pack(i)+1] = -1; // rank  //
        mpiPackInfo_send[id_tmp_pack(i)+2] = -1;     mpiPackInfo_recv[id_tmp_pack(i)+2] = -1; // tag   //
        mpiPackInfo_send[id_tmp_pack(i)+3] =  0;     mpiPackInfo_recv[id_tmp_pack(i)+3] =  0; // nleaf //
        mpiPackInfo_send[id_tmp_pack(i)+4] = -1;     mpiPackInfo_recv[id_tmp_pack(i)+4] = -1; // pack index //
    }

    for (int i=0; i<nprocs_send; i++) {
        const int rank_send = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;

        mpiPackInfo_send[id_tmp_pack(rank_send)+0] = 1;                                                    // flag  //
        mpiPackInfo_send[id_tmp_pack(rank_send)+1] = rank;                                                 // rank  //
        mpiPackInfo_send[id_tmp_pack(rank_send)+2] = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag;      // tag   //
        mpiPackInfo_send[id_tmp_pack(rank_send)+3] = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf;    // nleaf //
        mpiPackInfo_send[id_tmp_pack(rank_send)+4] = i;                                                    // pack index //
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Alltoall(mpiPackInfo_send, nelem, MPI_INT, mpiPackInfo_recv, nelem, MPI_INT, MPI_COMM_WORLD);


    int* id_sleaf_all = new int[num_procs*nleaf_global_max];
    int* id_rleaf_all = new int[num_procs*nleaf_global_max];
    for (int i=0; i<num_procs*nleaf_global_max; i++) { id_sleaf_all[i] = -1; id_rleaf_all[i] = -1; }

    for (int i=0; i<nprocs_send; i++) {
        const int  _target_rank = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;
        for (int j=0; j<mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf; j++) {
            id_sleaf_all[_target_rank*nleaf_global_max + j] = valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.id[j] / DefAMR::NN_LEAF;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Alltoall(id_sleaf_all, nleaf_global_max, MPI_INT, id_rleaf_all, nleaf_global_max, MPI_INT, MPI_COMM_WORLD);


    /// info recv ///
    int  nprocs_recv = 0;
    {
        for (int i=0; i<num_procs; i++) { if ( mpiPackInfo_recv[id_tmp_pack(i) + 0] == 1 ) { nprocs_recv++; } } // flag //
        mpiPackUnpackInfo.recvPackUnpackCommInfo.resize( nprocs_recv );

        int  idr = 0; int  sum_leafr = 0;
        for (int i=0; i<num_procs; i++) {
            if ( mpiPackInfo_recv[id_tmp_pack(i) + 0] == 1 ) {
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].rank     = mpiPackInfo_recv[id_tmp_pack(i) + 1]; // rank  //
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].tag      = mpiPackInfo_recv[id_tmp_pack(i) + 2]; // tag   //
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf    = mpiPackInfo_recv[id_tmp_pack(i) + 3]; // nleaf //
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].sum_leaf = sum_leafr;

                const int  _target_rank = mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].rank;
                const int  _target_tag  = mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].tag;
                if (_target_rank != i) { std::cout << "error : _target_rank != i\n"; exit(-1); }
                if (_target_tag  != i) { std::cout << "error : _target_tag  != i\n"; exit(-1); }

                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].id_leaf.resize(mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf);
                for (int j=0; j<mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf; j++) {
                    mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].id_leaf[j] = id_rleaf_all[ i*nleaf_global_max + j ];

                    if ( id_rleaf_all[ _target_rank*nleaf_global_max + j ] == -1 ) { std::cout << "error id_rleaf\n"; exit(-1); }
                }

                sum_leafr += mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf;
                idr++;
            }
        }

    }
    MPI_Barrier(MPI_COMM_WORLD);

    delete [] id_sleaf_all;
    delete [] id_rleaf_all;


    /// slist & rlist ///
    constexpr int nQ = LBM_velocity_model::nQ;
    {
        // leaf //
        mpiPackUnpackInfo.num_slist_leaf = 0;
        for (int i=0; i<nprocs_send; i++) {
            const int tmp_nleaf = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.size();
            mpiPackUnpackInfo.num_slist_leaf += tmp_nleaf;

            for (int j=0; j<tmp_nleaf; j++) {
                const int tmp_id = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf[j];
                mpiPackUnpackInfo.slist_leaf.push_back( tmp_id );
            }
        }

        mpiPackUnpackInfo.num_rlist_leaf = 0;
        for (int i=0; i<nprocs_recv; i++) {
            const int tmp_nleaf = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].id_leaf.size();
            mpiPackUnpackInfo.num_rlist_leaf += tmp_nleaf;

            for (int j=0; j<tmp_nleaf; j++) {
                const int tmp_id = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].id_leaf[j];
                mpiPackUnpackInfo.rlist_leaf.push_back( tmp_id );
            }
        }
        // check //
        if (mpiPackUnpackInfo.num_slist_leaf != mpiPackUnpackInfo.slist_leaf.size()) { std::cout << "error num_slist_leaf\n"; }
        if (mpiPackUnpackInfo.num_rlist_leaf != mpiPackUnpackInfo.rlist_leaf.size()) { std::cout << "error num_rlist_leaf\n"; }


        // val & lbm //
        const int num_slist_leaf = mpiPackUnpackInfo.num_slist_leaf;
        const int num_rlist_leaf = mpiPackUnpackInfo.num_rlist_leaf;
        allocate_MPIPackUnpackInfo(mpiPackUnpackInfo, num_slist_leaf, num_rlist_leaf);

        std::cout << "rank, lv = " << rank << ", " << lv << " : num_slist_leaf, num_rlist_leaf = " << num_slist_leaf << ", " << num_rlist_leaf << std::endl;


        // lambda //
        auto set_leaf = [](int num_leaf, int m, const int* list_leaf, int* list_val)
        {
            int tmp_id_list = 0;
            for (int i=0; i<num_leaf; i++) {
                const int tmp_offset_leaf = list_leaf[i]*m;
                for (int j=0; j<m; j++) {
                    list_val[tmp_id_list] = tmp_offset_leaf + j;
                    tmp_id_list++;
                }
            }
        };


        // val //
        mpiPackUnpackInfo.num_slist_val = num_slist_leaf * DefAMR::NN_LEAF;
        mpiPackUnpackInfo.num_rlist_val = num_rlist_leaf * DefAMR::NN_LEAF;

        set_leaf(num_slist_leaf, DefAMR::NN_LEAF, mpiPackUnpackInfo.slist_leaf.data(), mpiPackUnpackInfo.slist_val);
        set_leaf(num_rlist_leaf, DefAMR::NN_LEAF, mpiPackUnpackInfo.rlist_leaf.data(), mpiPackUnpackInfo.rlist_val);

        // lbm //
        mpiPackUnpackInfo.num_slist_lbm = num_slist_leaf * DefAMR::NN_LEAF*nQ;
        mpiPackUnpackInfo.num_rlist_lbm = num_rlist_leaf * DefAMR::NN_LEAF*nQ;

        set_leaf(num_slist_leaf, DefAMR::NN_LEAF*nQ, mpiPackUnpackInfo.slist_leaf.data(), mpiPackUnpackInfo.slist_lbm);
        set_leaf(num_rlist_leaf, DefAMR::NN_LEAF*nQ, mpiPackUnpackInfo.rlist_leaf.data(), mpiPackUnpackInfo.rlist_lbm);
    }

#if 0
//    for (int i=0; i<nprocs_send; i++) {
//        std::cout << "send : rank_src, dst = " << rank << ", " << mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank << " : offset = " << mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf 
//                  << " : slist               = ";
//        for (int j=0; j<mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf; j++) {
//            std::cout << mpiPackUnpackInfo.slist_leaf[mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf + j]  << ", ";
//        }
//        std::cout << std::endl << std::flush;
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    std::cout << std::endl << std::flush;
//    MPI_Barrier(MPI_COMM_WORLD);

//    for (int i=0; i<nprocs_send; i++) {
//        std::cout << "send : rank_src, dst = " << rank << ", " << mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank << " : offset = " << mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf 
//                  << " : mpiDatatypeInfo_dst = ";
//        for (int j=0; j<mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf; j++) {
//            std::cout << valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.id[j] / DefAMR::NN_LEAF  << ", ";
//        }
//        std::cout << std::endl << std::flush;
//    }
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    std::cout << std::endl << std::flush;
//    MPI_Barrier(MPI_COMM_WORLD);

//    int _jj;
//    _jj = 0;
    for (int i=0; i<nprocs_recv; i++) {
        std::cout << "recv : rank_src, dst = " << rank << ", " << mpiPackUnpackInfo.recvPackUnpackCommInfo[i].rank << " : offset = " << mpiPackUnpackInfo.recvPackUnpackCommInfo[i].sum_leaf 
                  << " : rlist               = ";
        for (int j=0; j<mpiPackUnpackInfo.recvPackUnpackCommInfo[i].nleaf; j++) {
            std::cout << mpiPackUnpackInfo.rlist_leaf[mpiPackUnpackInfo.recvPackUnpackCommInfo[i].sum_leaf + j]  << ", ";
//            std::cout << mpiPackUnpackInfo.rlist_leaf[mpiPackUnpackInfo.recvPackUnpackCommInfo[i].sum_leaf + j]  << "/ " << mpiPackUnpackInfo.rlist_leaf[_jj] << ", ";
//            _jj++;
        }
        std::cout << std::endl << std::flush;
    }


    exit(-1);
#endif
    MPI_Barrier(MPI_COMM_WORLD);

}


#if 0
void Tree::
_createMPIPackUnpack_lv(
          MPIPackUnpackInfo&                    mpiPackUnpackInfo,
    const int                                   lv,
    const std::vector<MPIDatatypeInfoSrcDst>&   valDatatypeInfoSrcDst
    )
{
    std::cout << __PRETTY_FUNCTION__ << " : lv = " << lv << std::endl;
    const int  rank      = rank_;
    const int  num_procs = num_procs_;

    /// info send ///
    const int  nprocs_send = valDatatypeInfoSrcDst.size();

    mpiPackUnpackInfo.sendPackUnpackCommInfo.resize( nprocs_send );

    {
        int sum_leaf = 0;
        for (int i=0; i<nprocs_send; i++) {
            mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank     = valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.rank;
            mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag      = rank;
            mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf    = valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.count;
            mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf = sum_leaf;     sum_leaf += mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf;

            for (int j=0; j<mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf; j++) {
                mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.push_back( valDatatypeInfoSrcDst[i].mpiDatatypeInfo_src.id[j] / DefAMR::NN_LEAF );
            }
        }
    }


    /// communication ///
    constexpr int nelem = 4;
    auto id_tmp_pack = [nelem](int i){ return nelem*i; };

    int  mpiPackInfo_send[num_procs*nelem];
    int  mpiPackInfo_recv[num_procs*nelem];

    for (int i=0; i<num_procs; i++) {
        mpiPackInfo_send[id_tmp_pack(i)+0] =  0;     mpiPackInfo_recv[id_tmp_pack(i)+0] =  0; // flag  //
        mpiPackInfo_send[id_tmp_pack(i)+1] = -1;     mpiPackInfo_recv[id_tmp_pack(i)+1] = -1; // rank  //
        mpiPackInfo_send[id_tmp_pack(i)+2] = -1;     mpiPackInfo_recv[id_tmp_pack(i)+2] = -1; // tag   //
        mpiPackInfo_send[id_tmp_pack(i)+3] =  0;     mpiPackInfo_recv[id_tmp_pack(i)+3] =  0; // nleaf //
    }

    for (int i=0; i<nprocs_send; i++) {
        const int rank_send = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;

        mpiPackInfo_send[id_tmp_pack(rank_send)+0] = 1;                                                    // flag  //
        mpiPackInfo_send[id_tmp_pack(rank_send)+1] = rank;                                                 // rank  //
        mpiPackInfo_send[id_tmp_pack(rank_send)+2] = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag;      // tag   //
        mpiPackInfo_send[id_tmp_pack(rank_send)+3] = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf;    // nleaf //
    }
    MPI_Barrier(MPI_COMM_WORLD);


    // check //
#if 0
    if (rank == 0) { std::cout << "befor put\n"; }
    for (int i=0; i<num_procs; i++) {
        std::cout << "send : rank, i =  " << rank << ", " << i << " : flag, rank, tag, nleaf = " << mpiPackInfo_send[id_tmp_pack(i)] << ", " << mpiPackInfo_send[id_tmp_pack(i)+1] << ", " << mpiPackInfo_send[id_tmp_pack(i)+2] << ", " << mpiPackInfo_send[id_tmp_pack(i)+3] << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << std::endl;
#endif


    MPI_Alltoall(mpiPackInfo_send, nelem, MPI_INT, mpiPackInfo_recv, nelem, MPI_INT, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);


#if 0
    if (rank == 0) { std::cout << "after put\n"; }
    // check //
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0; i<num_procs; i++) {
        std::cout << "recv : rank, i =  " << rank << ", " << i << " : flag, rank, tag, nleaf = " << mpiPackInfo_recv[id_tmp_pack(i)] << ", " << mpiPackInfo_recv[id_tmp_pack(i)+1] << ", " << mpiPackInfo_recv[id_tmp_pack(i)+2] << ", " << mpiPackInfo_recv[id_tmp_pack(i)+3] << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
#endif


    /// id leaf dst ///
    std::vector<MPIPackUnpackCommInfo>  tmpPackUnpackCommInfo;
    tmpPackUnpackCommInfo.resize( nprocs_send );
    for (int i=0; i<nprocs_send; i++) {
        for (int j=0; j<mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf; j++) {
            tmpPackUnpackCommInfo[i].id_leaf.push_back( valDatatypeInfoSrcDst[i].mpiDatatypeInfo_dst.id[j] / DefAMR::NN_LEAF );
        }
    }


    /// info recv ///
    int  nprocs_recv = 0;
    {
        for (int i=0; i<num_procs; i++) {
            if ( mpiPackInfo_recv[id_tmp_pack(i) + 0] == 1 ) { nprocs_recv++; }

#if 0 // check //
            if ( mpiPackInfo_recv[id_tmp_pack(i) + 0] == 1 && i != mpiPackInfo_recv[id_tmp_pack(i) + 1] ) {
                std::cout << "error : info recv : rank = " << rank  << " : i, tmpPackInfo_recv[] = " << i << ", " << mpiPackInfo_recv[id_tmp_pack(i) + 1] << std::endl;
                exit(-1);
            }
#endif
        }
        mpiPackUnpackInfo.recvPackUnpackCommInfo.resize( nprocs_recv );


        int  idr = 0; int  sum_leafr = 0;
        for (int i=0; i<num_procs; i++) {
            if ( mpiPackInfo_recv[id_tmp_pack(i) + 0] == 1 ) {
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].rank     = mpiPackInfo_recv[id_tmp_pack(i) + 1]; // rank  //
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].tag      = mpiPackInfo_recv[id_tmp_pack(i) + 2]; // tag   //
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf    = mpiPackInfo_recv[id_tmp_pack(i) + 3]; // nleaf //
                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].sum_leaf = sum_leafr;

                mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].id_leaf.resize( mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf );

                sum_leafr += mpiPackUnpackInfo.recvPackUnpackCommInfo[idr].nleaf;
                idr++;
            }
        }


        // id_leaf (recv) //
        MPI_Request  request_s[nprocs_send], request_r[nprocs_recv];
        MPI_Status   status_s[nprocs_send],  status_r[nprocs_recv];

        for (int i=0; i<nprocs_send; i++) {
//            const int*  sbuff = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.data();
            const int*  sbuff =                    tmpPackUnpackCommInfo[i].id_leaf.data();
            const int   snum  = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.size();
            const int   srank = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;
            const int   stag  = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag;

            MPI_Isend( sbuff, snum, MPI_INT, srank, stag, MPI_COMM_WORLD, request_s );
        }

        for (int i=0; i<nprocs_recv; i++) {
                  int*  rbuff = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].id_leaf.data();
            const int   rnum  = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].id_leaf.size();
            const int   rrank = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].rank;
            const int   rtag  = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].tag;

            MPI_Irecv( rbuff, rnum, MPI_INT, rrank, rtag, MPI_COMM_WORLD, request_r );
        }
        
        for (int i=0; i<nprocs_recv; i++) { MPI_Wait(request_r, status_r); }
        for (int i=0; i<nprocs_send; i++) { MPI_Wait(request_s, status_s); }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    /// slist & rlist ///
    constexpr int nQ = LBM_velocity_model::nQ;
    {
        // leaf //
        mpiPackUnpackInfo.num_slist_leaf = 0;
        for (int i=0; i<nprocs_send; i++) {
            const int tmp_nleaf = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf.size();
            mpiPackUnpackInfo.num_slist_leaf += tmp_nleaf;

            for (int j=0; j<tmp_nleaf; j++) {
                const int tmp_id = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].id_leaf[j];
                mpiPackUnpackInfo.slist_leaf.push_back( tmp_id );
            }
        }


        mpiPackUnpackInfo.num_rlist_leaf = 0;
        for (int i=0; i<nprocs_recv; i++) {
            const int tmp_nleaf = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].id_leaf.size();
            mpiPackUnpackInfo.num_rlist_leaf += tmp_nleaf;

//            std::cout << "recv : rank, i, tmp_nleaf " << rank << ", " << i << ", " << tmp_nleaf << " : "; // check //
            for (int j=0; j<tmp_nleaf; j++) {
                const int tmp_id = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].id_leaf[j];
                mpiPackUnpackInfo.rlist_leaf.push_back( tmp_id );

//                std::cout << tmp_id << ", "; // check //
            }
//            std::cout << std::endl; // check //
        }


        // val & lbm //
        const int num_slist_leaf = mpiPackUnpackInfo.num_slist_leaf;
        const int num_rlist_leaf = mpiPackUnpackInfo.num_rlist_leaf;
        allocate_MPIPackUnpackInfo(mpiPackUnpackInfo, num_slist_leaf, num_rlist_leaf);

        std::cout << "rank, lv = " << rank << ", " << lv << " : num_slist_leaf, num_rlist_leaf = " << num_slist_leaf << ", " << num_rlist_leaf << std::endl;


        // lambda //
        auto set_leaf = [](int num_leaf, int m, const int* list_leaf, int* list_val)
        {
            int tmp_id_list = 0;
            for (int i=0; i<num_leaf; i++) {
                const int tmp_offset_leaf = list_leaf[i]*m;
                for (int j=0; j<m; j++) {
                    list_val[tmp_id_list] = tmp_offset_leaf + j;
                    tmp_id_list++;
                }
            }
        };


        // val //
        mpiPackUnpackInfo.num_slist_val = num_slist_leaf * DefAMR::NN_LEAF;
        mpiPackUnpackInfo.num_rlist_val = num_rlist_leaf * DefAMR::NN_LEAF;

        set_leaf(num_slist_leaf, DefAMR::NN_LEAF, mpiPackUnpackInfo.slist_leaf.data(), mpiPackUnpackInfo.slist_val);
        set_leaf(num_rlist_leaf, DefAMR::NN_LEAF, mpiPackUnpackInfo.rlist_leaf.data(), mpiPackUnpackInfo.rlist_val);

        // lbm //
        mpiPackUnpackInfo.num_slist_lbm = num_slist_leaf * DefAMR::NN_LEAF*nQ;
        mpiPackUnpackInfo.num_rlist_lbm = num_rlist_leaf * DefAMR::NN_LEAF*nQ;

        set_leaf(num_slist_leaf, DefAMR::NN_LEAF*nQ, mpiPackUnpackInfo.slist_leaf.data(), mpiPackUnpackInfo.slist_lbm);
        set_leaf(num_rlist_leaf, DefAMR::NN_LEAF*nQ, mpiPackUnpackInfo.rlist_leaf.data(), mpiPackUnpackInfo.rlist_lbm);

    }

}
#endif


void Tree::
allocate_MPIPackUnpackInfo(
          MPIPackUnpackInfo&    mpiPackUnpackInfo,
    const int                   num_leaf_send,
    const int                   num_leaf_recv
    )
{
    constexpr int nQ   = LBM_velocity_model::nQ;

    FuncAllocate::allocate_value<int>(&mpiPackUnpackInfo.slist_val, num_leaf_send * DefAMR::NN_LEAF,      MemType::Managed);
    FuncAllocate::allocate_value<int>(&mpiPackUnpackInfo.slist_lbm, num_leaf_send * DefAMR::NN_LEAF * nQ, MemType::Managed);

    FuncAllocate::allocate_value<int>(&mpiPackUnpackInfo.rlist_val, num_leaf_recv * DefAMR::NN_LEAF,      MemType::Managed);
    FuncAllocate::allocate_value<int>(&mpiPackUnpackInfo.rlist_lbm, num_leaf_recv * DefAMR::NN_LEAF * nQ, MemType::Managed);


    FuncAllocate::fill_value<int>(mpiPackUnpackInfo.slist_val, 0.0, num_leaf_send * DefAMR::NN_LEAF,      MemType::Managed);
    FuncAllocate::fill_value<int>(mpiPackUnpackInfo.slist_lbm, 0.0, num_leaf_send * DefAMR::NN_LEAF * nQ, MemType::Managed);

    FuncAllocate::fill_value<int>(mpiPackUnpackInfo.rlist_val, 0.0, num_leaf_recv * DefAMR::NN_LEAF,      MemType::Managed);
    FuncAllocate::fill_value<int>(mpiPackUnpackInfo.rlist_lbm, 0.0, num_leaf_recv * DefAMR::NN_LEAF * nQ, MemType::Managed);
}


void Tree::
release_MPIPackUnpackInfo(MPIPackUnpackInfo&    mpiPackUnpackInfo)
{
    FuncAllocate::release_value<int>(mpiPackUnpackInfo.slist_val, MemType::Managed);
    FuncAllocate::release_value<int>(mpiPackUnpackInfo.slist_lbm, MemType::Managed);

    FuncAllocate::release_value<int>(mpiPackUnpackInfo.rlist_val, MemType::Managed);
    FuncAllocate::release_value<int>(mpiPackUnpackInfo.rlist_lbm, MemType::Managed);
}
