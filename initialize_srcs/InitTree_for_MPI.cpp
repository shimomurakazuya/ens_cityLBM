#include "InitTree_for_MPI.h"
#include "InitTree.h"
#include "defineAMR.h"
#include <algorithm>
#include <array>
#include <functional> // std::ref
#include <map>
#include <tuple>
#include "mpi_wrapper.hpp"
#include "divisor.hpp"
#include "range.hpp"
#include "runtime_error.hpp"


//void  InitTree_for_MPI::
//init_tree_uniform3d_div(
//          Tree& tree,
//    const Grid* grids,
//    const int   lv_max
//    )
//{
//    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    Tree  tree_global;
//    InitTree::init_tree_amr(tree_global, grids, lv_max);
////    InitTree::init_tree_cavity2d(tree_global, grids, lv_max);
////    InitTree::init_tree_flow_around_cube(tree_global, grids, lv_max);
//    MPI_Barrier(MPI_COMM_WORLD);
//
//
//    // allocate //
//    const int  nsize = tree_global.nodes().size();
////    std::cout << rank << " : " << nsize << std::endl;
//    localNode_.resize(nsize);
//
//    // set rank //
//    set_LocalNode(tree_global, grids);
//
//    // tree //
//    set_LocalTree(tree, tree_global, grids); // not check yet
//
//
//    // MPI information //
//    std::vector<MPIPutGetInfo>  tmp_mpiPutGetInfo;
//    createMPIPutGetInfo(tree_global, tmp_mpiPutGetInfo);
//    tree.copyMPIPutGetInfo(tmp_mpiPutGetInfo, true);
//}

void  InitTree_for_MPI::
init_tree_mpi_with_map(
          Tree& tree,
    const Grid* grids,
    const int   lv_max,
    const MapData& map
    )
{
    Tree  tree_global(tree.comm());
    InitTree::init_tree_amr_with_map(tree_global, grids, lv_max, map);
    MPI_Barrier(MPI_COMM_WORLD);


    // allocate //
    const int  nsize = tree_global.nodes().size();
    localNode_.resize(nsize);

    // set rank //
    set_LocalNode(tree_global, grids);

    // tree //
    set_LocalTree(tree, tree_global, grids); // not check yet


    // MPI information //
    std::vector<MPIPutGetInfo>  tmp_mpiPutGetInfo;
    createMPIPutGetInfo(tree_global, tmp_mpiPutGetInfo);
    tree.copyMPIPutGetInfo(tmp_mpiPutGetInfo, true);
}


void  InitTree_for_MPI::
set_LocalNode(
    const Tree& tree_global,
    const Grid* grids
    )
{
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    // basic information //
    set_lv          (localNode_, tree_global);
    set_cal_rank_div(localNode_, tree_global, grids);
    set_flags       (localNode_, tree_global);

    set_LocalNode_index(localNode_, num_allocated_localNode_); // write localNode_.id  and  num_allocated_localNode_ //

    // connection information //
    set_MPIGetInfo(localNode_);
    MPI_Barrier(MPI_COMM_WORLD);

    set_MPIPutInfo(localNode_, tree_global);
    MPI_Barrier(MPI_COMM_WORLD);
}


void  InitTree_for_MPI::
set_LocalTree(
          Tree& tree,
    const Tree& tree_global,
    const Grid* grids
    )
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int nsize = tree_global.nodes().size();
    const int num_allocated_localNode = num_allocated_localNode_;

    int  itmp = 0;
    int  mesh_offset[DefAMR::LV_MAX];
    for (int l=0; l<DefAMR::LV_MAX; l++) {
        mesh_offset[l] = 0;
    }

    tree.resize_node(num_allocated_localNode);
    for (int i=0; i<nsize; i++) {
        if ( !localNode_[i].allocate_flag ) { continue; }
        int  lid = localNode_[i].id;
        if (lid != itmp) { std::cout << __PRETTY_FUNCTION__ << " : error" << std::endl; }

        // index //
        tree.nodes(lid)->set_index(lid);

        // this //
        // position //
        tree.nodes(lid)->set_node_type( tree_global.nodes(i)->node_type() );
        tree.nodes(lid)->set_position( tree_global.nodes(i)->level(), tree_global.nodes(i)->x(), tree_global.nodes(i)->y(), tree_global.nodes(i)->z() );

        // mesh_offset //
        const int lv = tree_global.nodes(i)->level();
        tree.nodes(lid)->set_mesh_offset( mesh_offset[lv] );

        // connections //
        // parent //
        if ( tree_global.nodes(i)->parent()->node_type() == NodeTypes::LeafT ) {
            tree.nodes(lid)->set_parent( tree.nodes( localNode_[tree_global.nodes(i)->parent()->index()].id ) );
        }
        else {
            tree.nodes(lid)->set_parent( tree.external_node( tree_global.nodes(i)->parent()->node_type() ) );
        }

        // child //
        if ( tree_global.nodes(i)->child()->node_type() == NodeTypes::LeafT ) {
            tree.nodes(lid)->set_child( tree.nodes( localNode_[tree_global.nodes(i)->child()->index()].id ) );
        }
        else {
            tree.nodes(lid)->set_child( tree.external_node( tree_global.nodes(i)->child()->node_type() ) );
        }

        // neighbor //
        std::vector<std::array<int, 3>>  ineighbors;
        ineighbors.push_back({-1,  0,  0}); ineighbors.push_back({ 1,  0,  0});
        ineighbors.push_back({ 0, -1,  0}); ineighbors.push_back({ 0,  1,  0});
        ineighbors.push_back({ 0,  0, -1}); ineighbors.push_back({ 0,  0,  1});
        for (auto& in : ineighbors) {
            if ( tree_global.nodes(i)->neighbor( in[0],in[1],in[2] )->node_type() == NodeTypes::LeafT ) {
                tree.nodes(lid)->set_neighbor( tree.nodes( localNode_[tree_global.nodes(i)->neighbor( in[0],in[1],in[2] )->index()].id ), in[0],in[1],in[2] );
            }
            else {
                tree.nodes(lid)->set_neighbor( tree.external_node( tree_global.nodes(i)->neighbor( in[0],in[1],in[2] )->node_type() ), in[0],in[1],in[2] );
            }
        }

        // node flags //
        tree.nodes(lid)->nodeCalFlags().set_flag_Cal( tree_global.nodes(i)->nodeCalFlags().Cal() );
        if ( !localNode_[i].cal_flag ) { tree.nodes(lid)->nodeCalFlags().set_flag_Cal( false ); }


        // MPI //
        tree.nodes(lid)->nodeCalFlags().set_flag_putMPI( localNode_[i].flag_put );
        tree.nodes(lid)->nodeCalFlags().set_flag_getMPI( localNode_[i].flag_get );

        // iteration //
        itmp++; // check //
        mesh_offset[lv] += DefAMR::NN_LEAF;
    }


    // update node flags //
    const int  n_leaf  = tree.number_of_nodes();
    for (int i=0; i<n_leaf; i++) {
        Node*  node = tree.nodes(i);

        // put : L2X //
        if ( node->nodeCalFlags().Cal() ) {
            // parent : L2C //
            if ( node->parent()->node_type() == NodeTypes::LeafT ) {
                if ( node->parent()->nodeCalFlags().Cal() ) {
                    std::cout << "error : parent flag is calculate" << std::endl;
                }
                else {
                    node->nodeCalFlags().set_flag_L2C(true);
                }
            }
            // child : L2F //
            if ( node->child()->node_type() == NodeTypes::LeafT ) {
                if (   node->child()->neighbor(0,0,0)->nodeCalFlags().Cal()
                    || node->child()->neighbor(1,0,0)->nodeCalFlags().Cal()
                    || node->child()->neighbor(0,1,0)->nodeCalFlags().Cal()
                    || node->child()->neighbor(1,1,0)->nodeCalFlags().Cal()
                    || node->child()->neighbor(0,0,1)->nodeCalFlags().Cal()
                    || node->child()->neighbor(1,0,1)->nodeCalFlags().Cal()
                    || node->child()->neighbor(0,1,1)->nodeCalFlags().Cal()
                    || node->child()->neighbor(1,1,1)->nodeCalFlags().Cal() ) {
                    std::cout << "error : child flag is calculate" << std::endl;
                }
                else {
                    node->nodeCalFlags().set_flag_L2F(true);
                }
            }
        }

        // get : X2L //
        if ( !node->nodeCalFlags().Cal() && node->node_type() == NodeTypes::LeafT ) {
            // parent : C2L //
            if ( node->parent()->node_type() == NodeTypes::LeafT ) {
                node->nodeCalFlags().set_flag_C2L(true);
            }
            // child : F2L //
            if ( node->child()->node_type() == NodeTypes::LeafT ) {
                node->nodeCalFlags().set_flag_F2L(true);
            }
        }
    }

}


void  InitTree_for_MPI::
createMPIPutGetInfo(const Tree& tree_global, std::vector<MPIPutGetInfo>& mpiPutGetInfo)
const
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    mpiPutGetInfo.resize(num_allocated_localNode_);
    mpiPutGetInfo.shrink_to_fit();

    for (int i=0; i<(int)tree_global.nodes().size(); i++) {
        if ( !localNode_[i].allocate_flag ) { continue; }
        const int  lid = localNode_[i].id;


        mpiPutGetInfo[lid].flag_put    = localNode_[i].flag_put;
        mpiPutGetInfo[lid].num_vec_put = localNode_[i].rank_put.size();
        mpiPutGetInfo[lid].rank_put    = localNode_[i].rank_put;
        mpiPutGetInfo[lid].id_put      = localNode_[i].id_put;
        mpiPutGetInfo[lid].offset_put  = localNode_[i].offset_put;

        mpiPutGetInfo[lid].flag_get    = localNode_[i].flag_get;
        mpiPutGetInfo[lid].rank_get    = localNode_[i].rank_get;
        mpiPutGetInfo[lid].id_get      = localNode_[i].id_get;
        mpiPutGetInfo[lid].offset_get  = localNode_[i].offset_get;
    }
}


void  InitTree_for_MPI::
set_lv(std::vector<LocalNode>& localNode, const Tree& tree_global)
{
    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        localNode[i].lv = node->level();
    }
}


void  InitTree_for_MPI::
set_cal_rank_div(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
//    set_cal_rank_div1d(localNode, tree_global, grids);
//    set_cal_rank_div2d(localNode, tree_global, grids);
//    set_cal_rank_div3d(localNode, tree_global, grids);

//    set_cal_rank_div_nc3d(localNode, tree_global, grids);
//    set_cal_rank_div_channel(localNode, tree_global, grids);
//    set_cal_rank_div_citybox(localNode, tree_global, grids);

//    set_cal_rank_div_2d_block(localNode, tree_global, grids, 1, 1);
    set_cal_rank_div_2d_block(localNode, tree_global, grids, 2, 2);
}


#if 0
void  InitTree_for_MPI::
set_cal_rank_div1d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
    const auto ncpu = comm_.col_vector().size();


#if 1 // x //
    int  mpi_nx[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        mpi_nx[i] = grids[i].nx().x() / ncpu;
    }

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };

        const int  cal_rank = pos[0] / mpi_nx[lv];

        localNode[i].cal_rank = cal_rank;
    }
#endif

#if 0 // z //
    int  mpi_nz[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        mpi_nz[i] = grids[i].nx().z() / ncpu;
    }

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };

        const int  cal_rank = pos[2] / mpi_nz[lv];

        localNode[i].cal_rank = cal_rank;
    }
#endif
}


void  InitTree_for_MPI::
set_cal_rank_div2d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
    const auto ncpu = comm_.col_vector().size();

    const int  ncpux = (int)sqrt(ncpu);
    const int  ncpuy = ncpu/ncpux;

    int  mpi_nx[DefAMR::LV_MAX];
    int  mpi_ny[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        mpi_nx[i] = grids[i].nx().x() / ncpux;
        mpi_ny[i] = grids[i].nx().y() / ncpuy;
    }

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };

        const int  _cal_rank_x = pos[0] / mpi_nx[lv];
        const int  _cal_rank_y = pos[1] / mpi_ny[lv];
        const int  cal_rank_x = (_cal_rank_x >= ncpux) ? ncpux - 1 : _cal_rank_x;
        const int  cal_rank_y = (_cal_rank_y >= ncpuy) ? ncpuy - 1 : _cal_rank_y;
        const int  cal_rank = cal_rank_x + ncpux*cal_rank_y;

//        std::cout << "cal_rank, ncpu = " << cal_rank << ", " << ncpu << std::endl;
        if ( !(cal_rank >= 0 && cal_rank <= ncpu-1) ) {
            std::cout << __PRETTY_FUNCTION__ 
                      << " : error cal_rank : ncpux, ncpuy, mpi_nx, mpi_ny, cal_rank_x, cal_rank_y, cal_rank = " 
                      << ncpux << ", " << ncpuy << ", " << mpi_nx[lv] << ", " << mpi_ny[lv] << ", " << cal_rank_x << ", " << cal_rank_y << ", " << cal_rank << std::endl;
        }

        localNode[i].cal_rank = cal_rank;
    }

}
#endif

void  InitTree_for_MPI::
set_cal_rank_div_2d_block(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids, int bx, int by)
{
    const auto ncpu = comm_.col_vector().size();

    /// decide ncpu_x, ncpu_y
#if 0 
    //// a brute-force map for BDEC (2021/07/16)
    const std::map<int, std::tuple<int, int>> ncpu_xy_map {
        // ss 40 -- 360
        {40, {10, 4}},
        {80, {10, 8}},
        {120, {12, 10}},
        {240, {20, 12}},
        {360, {20, 18}},
        // ss 36 -- 324
        {36, {6, 6}},
        {64, {8, 8}},
        {144, {12, 12}},
        {256, {16, 16}},
        {324, {18, 18}}
    };
    const auto [ncpu_x, ncpu_y] = ncpu_xy_map.at(ncpu);
#else
    int ncpu_x = 1, ncpu_y = 1;
    for(int res = ncpu; res > 1;) {
        for(const auto& m: util::divisor(res)) {
            if(m > 1) {
                int& nxy = (ncpu_x < ncpu_y ? std::ref(ncpu_x) : std::ref(ncpu_y));
                nxy *= m;
                res /= m;
                break;
            }
        }
    }
#endif
    const int nb_x = ncpu_x / bx;
    const int nb_y = ncpu_y / by;

    if(util::mpi(MPI_COMM_WORLD).rank() == 0) {
        std::cout << "cal_rank_div = " << ncpu_x << ", " << ncpu_y << std::endl;
        std::cout << " with block =  " << bx << ", " << by << std::endl;
    }

    runtime_assert(ncpu == nb_x * bx * nb_y * by, "InvalidConfiguration: ncpu is not divisible by such blocking");

    std::map<int, int> rank_map;
    for(const auto j_outer: util::irange(nb_y)) {
    for(const auto i_outer: util::irange(nb_x)) {
        for(const auto j_inner: util::irange(by)) {
        for(const auto i_inner: util::irange(bx)) {
            const auto i = i_inner + i_outer*bx;
            const auto j = j_inner + j_outer*by;
            const auto ij = i + j * ncpu_x;
            const auto id_block = i_inner + j_inner*bx + by*bx*(i_outer + j_outer*nb_x);
            rank_map[ij] = id_block;
        }
        }
    }
    }


    const int  nsize = tree_global.nodes().size();
    int imax = 0, jmax = 0;
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };
        const int  pos_nx = grids[lv].nx().x();
        const int  pos_ny = grids[lv].nx().y();

        const int  mpi_i = pos[0] * ncpu_x / pos_nx;
        const int  mpi_j = pos[1] * ncpu_y / pos_ny;
        const int  cal_rank = rank_map.at(mpi_i + mpi_j * ncpu_x);
        runtime_assert(0<=cal_rank && cal_rank<ncpu, "InternalError: invalid rank mapping");

        localNode[i].cal_rank = cal_rank;
        /// test
        imax = std::max(mpi_i, imax);
        jmax = std::max(mpi_j, jmax);
    }

}

void InitTree_for_MPI::
set_cal_rank_div_citybox(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
    const auto ncpu = comm_.col_vector().size();
    /// decide ncpu_x, ncpu_y
    int ncpu_x = 1, ncpu_y = 1;
    for(int res = ncpu; res > 1;) {
        for(const auto& m: util::divisor(res)) {
            if(m > 1) {
                int& nxy = (ncpu_x < ncpu_y ? std::ref(ncpu_x) : std::ref(ncpu_y));
                nxy *= m;
                res /= m;
                break;
            }
        }
    }
    if(util::mpi(MPI_COMM_WORLD).rank() == 0) {
        std::cout << "cal_rank_div = " << ncpu_x << ", " << ncpu_y << std::endl;
    }

    /// solve mesh_finest.{west,east,south,north}
    double west = NAN, east = NAN, south = NAN, north = NAN;
    auto&& min = [](double a, double b) { return (a<b) ? a : b; }; // allow nan
    auto&& max = [](double a, double b) { return (a>b) ? a : b; }; // allow nan
    for(const auto& node: tree_global.nodes()) {
        if(node.level() == DefAMR::LV_MAX-1) {
            west  = min(west , node.x());
            east  = max(east , node.y());
            south = min(south, node.x());
            north = max(north, node.y());
        }
    }

    /// domain decomposition
    auto&& floorceil = [](int rank_tmp, int rank_begin, int rank_end) -> int {
        int ret = rank_tmp;
        if(ret < rank_begin) { ret = rank_begin; }
        if(ret >= rank_end) { ret = rank_end - 1; }
        return ret;
    };
    for(std::intptr_t i=0; i<tree_global.nodes().size(); i++) {
        const auto& node = tree_global.nodes().at(i);
        const int lv = node.level();
        const int n_lv = std::pow(2, DefAMR::LV_MAX-1-lv);
        const double coo_x = n_lv * node.x();
        const double coo_y = n_lv * node.y();
        const int rank_x_tmp = (ncpu_x * (coo_x - west ))/(east  - west );
        const int rank_y_tmp = (ncpu_y * (coo_y - south))/(north - south);
        const int rank_x = floorceil(rank_x_tmp, 0, ncpu_x);
        const int rank_y = floorceil(rank_y_tmp, 0, ncpu_y);
        const int cal_rank = rank_x + ncpu_x * rank_y;
        localNode.at(i).cal_rank = cal_rank;
    }
}


void  InitTree_for_MPI::
set_cal_rank_div3d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
    const auto ncpu = comm_.col_vector().size();

    const int  ncpuz  = cbrt((double)ncpu);
    const int  ncpuxy = ncpu/ncpuz;

    const int  ncpux = (int)sqrt(ncpuxy);
    const int  ncpuy = ncpuxy/ncpux;

    int  mpi_nx[DefAMR::LV_MAX];
    int  mpi_ny[DefAMR::LV_MAX];
    int  mpi_nz[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        mpi_nx[i] = grids[i].nx().x() / ncpux;
        mpi_ny[i] = grids[i].nx().y() / ncpuy;
        mpi_nz[i] = grids[i].nx().z() / ncpuz;
    }

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };

        const int  _cal_rank_x = pos[0] / mpi_nx[lv];
        const int  _cal_rank_y = pos[1] / mpi_ny[lv];
        const int  _cal_rank_z = pos[2] / mpi_nz[lv];

        const int  cal_rank_x = (_cal_rank_x >= ncpux) ? ncpux - 1 : _cal_rank_x;
        const int  cal_rank_y = (_cal_rank_y >= ncpuy) ? ncpuy - 1 : _cal_rank_y;
        const int  cal_rank_z = (_cal_rank_z >= ncpuz) ? ncpuz - 1 : _cal_rank_z;
        const int  cal_rank = cal_rank_x + ncpux*cal_rank_y + ncpux*ncpuy*cal_rank_z;

//        std::cout << "cal_rank, ncpu = " << cal_rank << ", " << ncpu << std::endl;
        if ( !(cal_rank >= 0 && cal_rank <= ncpu-1) ) {
            std::cout << __PRETTY_FUNCTION__ 
                      << " : error cal_rank : ncpux, ncpuy, mpi_nx, mpi_ny, cal_rank_x, cal_rank_y, cal_rank = " 
                      << ncpux << ", " << ncpuy << ", " << mpi_nx[lv] << ", " << mpi_ny[lv] << ", " << cal_rank_x << ", " << cal_rank_y << ", " << cal_rank << std::endl;
        }

        localNode[i].cal_rank = cal_rank;
    }

}


void  InitTree_for_MPI::
set_cal_rank_div_nc3d(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
    const auto ncpu = comm_.col_vector().size();

    const int  ncpuz  = 4;
    const int  ncpuxy = ncpu/ncpuz;

    const int  ncpux = (int)sqrt(ncpuxy);
    const int  ncpuy = ncpuxy/ncpux;

    int  mpi_nx[DefAMR::LV_MAX];
    int  mpi_ny[DefAMR::LV_MAX];
    int  mpi_nz[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        mpi_nx[i] = grids[i].nx().x() / ncpux;
        mpi_ny[i] = grids[i].nx().y() / ncpuy;
        mpi_nz[i] = grids[i].nx().z() / ncpuz;
    }

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };

        const int  _cal_rank_x = pos[0] / mpi_nx[lv];
        const int  _cal_rank_y = pos[1] / mpi_ny[lv];
        const int  _cal_rank_z = pos[2] / mpi_nz[lv];

        const int  cal_rank_x = (_cal_rank_x >= ncpux) ? ncpux - 1 : _cal_rank_x;
        const int  cal_rank_y = (_cal_rank_y >= ncpuy) ? ncpuy - 1 : _cal_rank_y;
        const int  cal_rank_z = (_cal_rank_z >= ncpuz) ? ncpuz - 1 : _cal_rank_z;
        const int  cal_rank = cal_rank_x + ncpux*cal_rank_y + ncpux*ncpuy*cal_rank_z;

//        std::cout << "cal_rank, ncpu = " << cal_rank << ", " << ncpu << std::endl;
        if ( !(cal_rank >= 0 && cal_rank <= ncpu-1) ) {
            std::cout << __PRETTY_FUNCTION__ 
                      << " : error cal_rank : ncpux, ncpuy, mpi_nx, mpi_ny, cal_rank_x, cal_rank_y, cal_rank = " 
                      << ncpux << ", " << ncpuy << ", " << mpi_nx[lv] << ", " << mpi_ny[lv] << ", " << cal_rank_x << ", " << cal_rank_y << ", " << cal_rank << std::endl;
        }

        localNode[i].cal_rank = cal_rank;
    }

}


void  InitTree_for_MPI::
set_cal_rank_div_channel(std::vector<LocalNode>& localNode, const Tree& tree_global, const Grid* grids)
{
    const auto ncpu = comm_.col_vector().size();

    const int  _ncpu_tmp  = 1;
    const int  ncpuxy = ncpu/_ncpu_tmp;

    const int  _ncpux = (int)sqrt(ncpuxy);
    const int  ncpuy = ncpuxy/_ncpux;
    const int  ncpux = _ncpux*_ncpu_tmp;
    const int  ncpuz = 1;

    int  mpi_nx[DefAMR::LV_MAX];
    int  mpi_ny[DefAMR::LV_MAX];
    int  mpi_nz[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        mpi_nx[i] = grids[i].nx().x() / ncpux;
        mpi_ny[i] = grids[i].nx().y() / ncpuy;
        mpi_nz[i] = grids[i].nx().z() / ncpuz;
    }

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        const Node* node = tree_global.nodes(i);

        const int  lv    = node->level();
        const int  pos[] = { node->x(), node->y(), node->z() };

        const int  _cal_rank_x = pos[0] / mpi_nx[lv];
        const int  _cal_rank_y = pos[1] / mpi_ny[lv];
        const int  _cal_rank_z = pos[2] / mpi_nz[lv];

        const int  cal_rank_x = (_cal_rank_x >= ncpux) ? ncpux - 1 : _cal_rank_x;
        const int  cal_rank_y = (_cal_rank_y >= ncpuy) ? ncpuy - 1 : _cal_rank_y;
        const int  cal_rank_z = (_cal_rank_z >= ncpuz) ? ncpuz - 1 : _cal_rank_z;
        const int  cal_rank = cal_rank_x + ncpux*cal_rank_y + ncpux*ncpuy*cal_rank_z;

//        std::cout << "cal_rank, ncpu = " << cal_rank << ", " << ncpu << std::endl;
        if ( !(cal_rank >= 0 && cal_rank <= ncpu-1) ) {
            std::cout << __PRETTY_FUNCTION__ 
                      << " : error cal_rank : ncpux, ncpuy, mpi_nx, mpi_ny, cal_rank_x, cal_rank_y, cal_rank = " 
                      << ncpux << ", " << ncpuy << ", " << mpi_nx[lv] << ", " << mpi_ny[lv] << ", " << cal_rank_x << ", " << cal_rank_y << ", " << cal_rank << std::endl;
        }

        localNode[i].cal_rank = cal_rank;
    }
}


void  InitTree_for_MPI::
set_flags(std::vector<LocalNode>& localNode, const Tree& tree_global)
{
    const auto rank = comm_.col_vector().rank();

    const int  nsize = tree_global.nodes().size();
    for (int i=0; i<nsize; i++) {
        localNode[i].cal_flag = node_is_cal(localNode[i].cal_rank, rank);
    }
    for (int i=0; i<nsize; i++) {
        localNode[i].halo_flag = ( neighbor_is_cal(localNode, tree_global.nodes(i)) && !localNode[i].cal_flag );
    }
    for (int i=0; i<nsize; i++) {
        localNode[i].allocate_flag = neighbor_is_cal2(localNode, tree_global.nodes(i));
    }
}


bool  InitTree_for_MPI::
node_is_cal(const int cal_rank, const int rank)
{
    return  (cal_rank == rank) ? true : false;
}


bool  InitTree_for_MPI::
neighbor_is_cal(std::vector<LocalNode>& localNode, const Node* node)
{
    const int  imin = -1;
    const int  jmin = -1;
    const int  kmin = -1;
    const int  imax =  1;
    const int  jmax =  1;
    const int  kmax =  1;


    for (int kk=kmin; kk<=kmax; kk++) {
    for (int jj=jmin; jj<=jmax; jj++) {
    for (int ii=imin; ii<=imax; ii++) {
        if ( !node->neighbor(ii, jj, kk)->nodeCalFlags().Cal() ) { continue; }

        if (localNode[ node->neighbor(ii, jj, kk)->index() ].cal_flag == true) {
            return true;
        }
    }
    }
    }

    return false;
}


bool  InitTree_for_MPI::
neighbor_is_cal2(std::vector<LocalNode>& localNode, const Node* node)
{
    const int  lv = node->level();
    const int  pos[] = { node->x(), node->y(), node->z() };

    const int  imin = ( lv != 0 && pos[0]%2 == 1 ) ? -2 : -1;
    const int  jmin = ( lv != 0 && pos[1]%2 == 1 ) ? -2 : -1;
    const int  kmin = ( lv != 0 && pos[2]%2 == 1 ) ? -2 : -1;

    const int  imax = ( lv != 0 && pos[0]%2 == 0 ) ?  2 :  1;
    const int  jmax = ( lv != 0 && pos[1]%2 == 0 ) ?  2 :  1;
    const int  kmax = ( lv != 0 && pos[2]%2 == 0 ) ?  2 :  1;


    for (int kk=kmin; kk<=kmax; kk++) {
    for (int jj=jmin; jj<=jmax; jj++) {
    for (int ii=imin; ii<=imax; ii++) {
        if ( !node->neighbor(ii, jj, kk)->nodeCalFlags().Cal() ) { continue; }

        if (localNode[ node->neighbor(ii, jj, kk)->index() ].cal_flag == true) {
            return true;
        }
    }
    }
    }

    return false;
}


void  InitTree_for_MPI::
set_LocalNode_index(std::vector<LocalNode>& localNode, int& num_allocated_localNode)
{
    int  lid = 0;
    int  mesh_offset[DefAMR::LV_MAX];
    for (int l=0; l<DefAMR::LV_MAX; l++) { mesh_offset[l] = 0; }


    const int  nsize = localNode.size();
    for (int i=0; i<nsize; i++) {
        localNode[i].idg = i;

        const int lv = localNode[i].lv;
        if (localNode[i].allocate_flag) {
            localNode[i].id = lid;
            lid++;

            localNode[i].mesh_offset = mesh_offset[lv];
            mesh_offset[lv] += DefAMR::NN_LEAF;
        }
        else {
            localNode[i].id = -1;
            localNode[i].mesh_offset = -1;
        }
    }

    num_allocated_localNode = lid;
}


void  InitTree_for_MPI::
set_MPIGetInfo(std::vector<LocalNode>& localNode)
{
    const auto comm_ens = comm_.col_vector();
    const auto rank = comm_ens.rank();
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  _nsize = localNode.size();
    const int nsize = comm_.col_vector().reduce_max(_nsize);

    int* id_local     = new int[nsize];
    int* offset_local = new int[nsize];
    for (int i=0; i<nsize; i++) {
        id_local    [i] = localNode[i].id;
        offset_local[i] = localNode[i].mesh_offset;
    }

    int* lid_dst    = new int[nsize];
    int* offset_dst = new int[nsize];
    for (int i=0; i<nsize; i++) {
        lid_dst[i]    = id_local[i];
        offset_dst[i] = offset_local[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);


    // lid_dst //
    if (rank == 0) { std::cout << "win_create : lid_dst" << std::endl; }

    MPI_Win  win_id_get;
    MPI_Win_create(id_local, sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, comm_ens.comm(), &win_id_get);
    MPI_Win_fence(0, win_id_get);

    for (int i=0; i<nsize; i++) { // @ cal_flag == false && allocate_flag == true //
        if ( !localNode[i].halo_flag ) { continue; }

        const int  rank_dst = localNode[i].cal_rank;
        if (rank_dst == rank) { continue; }

        MPI_Get(&lid_dst[i], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_id_get);
    }
    MPI_Win_fence(0, win_id_get);
    MPI_Win_free(&win_id_get);
    // lid_dst //


    // offset_dst //
    if (rank == 0) { std::cout << "win_create : offset_dst" << std::endl; }

    MPI_Win  win_offset_get;
    MPI_Win_create(offset_local, sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, comm_ens.comm(), &win_offset_get);
    MPI_Win_fence(0, win_offset_get);

    for (int i=0; i<nsize; i++) { // @ cal_flag == false && allocate_flag == true //
        if ( !localNode[i].halo_flag ) { continue; }

        const int  rank_dst = localNode[i].cal_rank;
        if (rank_dst == rank) { continue; }

        MPI_Get(&offset_dst[i], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_offset_get);
    }
    MPI_Win_fence(0, win_offset_get);
    MPI_Win_free(&win_offset_get);
    // offset_dst //

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) { std::cout << "update" << std::endl; }

    for (int i=0; i<nsize; i++) {
        if ( !localNode[i].halo_flag ) {
            localNode[i].flag_get   = false;
            localNode[i].id_get     = -1;
            localNode[i].offset_get = -1;
            localNode[i].rank_get   = rank;
            continue;
        }
        const int  rank_dst = localNode[i].cal_rank;
        if (rank_dst == rank) {
            localNode[i].flag_get   = false;
            localNode[i].id_get     = -1;
            localNode[i].offset_get = -1;
            localNode[i].rank_get   = rank;
            continue;
        }

        localNode[i].flag_get   = true;
        localNode[i].id_get     = lid_dst[i];
        localNode[i].offset_get = offset_dst[i];
        localNode[i].rank_get   = rank_dst;
    }

    delete [] id_local;
    delete [] offset_local;

    delete [] lid_dst;
    delete [] offset_dst;
}


//void  InitTree_for_MPI::
//_set_MPIGetInfo(std::vector<LocalNode>& localNode)
//{
//    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//
//    const int  _nsize = localNode.size();
//    int  nsize;  MPI_Allreduce(&_nsize, &nsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//
//    std::vector<int>  id_local(nsize);
//    std::vector<int>  offset_local(nsize);
//    for (int i=0; i<nsize; i++) {
//        id_local    [i] = localNode[i].id;
//        offset_local[i] = localNode[i].mesh_offset;
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if (rank == 0) { std::cout << "win_create" << std::endl; }
//
//    MPI_Win  win_id_get;
//    MPI_Win_create(id_local.    data(), sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_id_get);
//    MPI_Win_fence(0, win_id_get);
//
//    MPI_Win  win_offset_get;
//    MPI_Win_create(offset_local.data(), sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_offset_get);
//    MPI_Win_fence(0, win_offset_get);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if (rank == 0) { std::cout << "lid_dst, offset_dst" << std::endl; }
//
//    std::vector<int>  lid_dst(nsize);
//    std::vector<int>  offset_dst(nsize);
//    for (int i=0; i<nsize; i++) {
//        lid_dst[i]    = id_local[i];
//        offset_dst[i] = offset_local[i];
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//
//    if (rank == 0) { std::cout << "mpi_get" << std::endl; }
//
//
//    // @ cal_flag == false && allocate_flag == true //
//    for (int i=0; i<nsize; i++) {
//        if ( !localNode[i].halo_flag ) { continue; }
//
//        const int  rank_dst = localNode[i].cal_rank;
//        if (rank_dst == rank) { continue; }
//
//        MPI_Get(&lid_dst   [i], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_id_get);
//        MPI_Get(&offset_dst[i], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_offset_get);
//    }
//    MPI_Win_fence(0, win_id_get);
//    MPI_Win_fence(0, win_offset_get);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    MPI_Win_free(&win_id_get);
//    MPI_Win_free(&win_offset_get);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    if (rank == 0) { std::cout << "update" << std::endl; }
//
//    for (int i=0; i<nsize; i++) {
//        if ( !localNode[i].halo_flag ) {
//            localNode[i].flag_get   = false;
//            localNode[i].id_get     = -1;
//            localNode[i].offset_get = -1;
//            localNode[i].rank_get   = rank;
//            continue;
//        }
//        const int  rank_dst = localNode[i].cal_rank;
//        if (rank_dst == rank) {
//            localNode[i].flag_get   = false;
//            localNode[i].id_get     = -1;
//            localNode[i].offset_get = -1;
//            localNode[i].rank_get   = rank;
//            continue;
//        }
//
//        localNode[i].flag_get   = true;
//        localNode[i].id_get     = lid_dst[i];
//        localNode[i].offset_get = offset_dst[i];
//        localNode[i].rank_get   = rank_dst;
//    }
//
//}


void  InitTree_for_MPI::
set_MPIPutInfo(std::vector<LocalNode>& localNode, const Tree& tree_global)
{
    const auto comm_ens = comm_.col_vector();
    const auto rank = comm_ens.rank();
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  _nsize = localNode.size();
    const int nsize = comm_.col_vector().reduce_max(_nsize);

    int* id_local     = new int[nsize];
    int* offset_local = new int[nsize];
    for (int i=0; i<nsize; i++) {
        id_local[i] = localNode[i].id;
        offset_local[i] = localNode[i].mesh_offset;
    }
    MPI_Barrier(MPI_COMM_WORLD);


    auto check_duplicate_list = [](std::vector<int> list, int val)
    {
        for (auto& elem : list) { if (elem == val) { return true; } }
        return  false;
    };


    // id_local //
    MPI_Win  win_id_get;
    MPI_Win_create(id_local, sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, comm_ens.comm(), &win_id_get);
    MPI_Win_fence(0, win_id_get);

    // offset_local //
    MPI_Win  win_offset_get;
    MPI_Win_create(offset_local, sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, comm_ens.comm(), &win_offset_get);
    MPI_Win_fence(0, win_offset_get);

    MPI_Barrier(MPI_COMM_WORLD);

    // @ cal_flag == true //
    int* lid_dst    = new int[nsize*27];
    int* offset_dst = new int[nsize*27];

    int  icount;
    icount = 0;
    for (int i=0; i<nsize; i++) {
        if ( !localNode[i].cal_flag ) {
            localNode[i].flag_put = false;
            continue;
        }
        const Node* node = tree_global.nodes(i);

        // remove_rank //
        std::vector<int>  remove_rank;  remove_rank.push_back(rank);
        for(int kk=-1; kk<=1; kk++) {
        for(int jj=-1; jj<=1; jj++) {
        for(int ii=-1; ii<=1; ii++) {
            if ( !node->neighbor(ii, jj, kk)->nodeCalFlags().Cal() ) { continue; }

            const int  id_neighbor = node->neighbor(ii, jj, kk)->index();
            const int  rank_dst = localNode[id_neighbor].cal_rank;

            if ( check_duplicate_list(remove_rank, rank_dst) ) { continue; }
            remove_rank.push_back(rank_dst);

            // MPI_Get //
            MPI_Get(&lid_dst   [icount], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_id_get);
            MPI_Get(&offset_dst[icount], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_offset_get);
            icount++;
        }
        }
        }
    }
    MPI_Win_fence(0, win_id_get);
    MPI_Win_fence(0, win_offset_get);

    MPI_Win_free(&win_id_get);
    MPI_Win_free(&win_offset_get);

    MPI_Barrier(MPI_COMM_WORLD);

    icount = 0;
    for (int i=0; i<nsize; i++) {
        if ( !localNode[i].cal_flag ) {
            localNode[i].flag_put = false;
            continue;
        }
        const Node* node = tree_global.nodes(i);

        // remove_rank //
        std::vector<int>  remove_rank;  remove_rank.push_back(rank);
        for(int kk=-1; kk<=1; kk++) {
        for(int jj=-1; jj<=1; jj++) {
        for(int ii=-1; ii<=1; ii++) {
            if ( !node->neighbor(ii, jj, kk)->nodeCalFlags().Cal() ) { continue; }

            const int  id_neighbor = node->neighbor(ii, jj, kk)->index();
            const int  rank_dst = localNode[id_neighbor].cal_rank;

            if ( check_duplicate_list(remove_rank, rank_dst) ) { continue; }
            remove_rank.push_back(rank_dst);

            // update //
            localNode[i].flag_put = true;
            localNode[i].id_put.    push_back( lid_dst[icount] );
            localNode[i].offset_put.push_back( offset_dst[icount] );
            localNode[i].rank_put.  push_back( rank_dst );

            icount++;
        }
        }
        }
    }

    delete [] id_local;
    delete [] offset_local;

    delete [] lid_dst;
    delete [] offset_dst;
}


//void  InitTree_for_MPI::
//_set_MPIPutInfo(std::vector<LocalNode>& localNode, const Tree& tree_global)
//{
//    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
//
//    const int  _nsize = localNode.size();
//    int  nsize;  MPI_Allreduce(&_nsize, &nsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//
//    std::vector<int>  id_local(nsize);
//    std::vector<int>  offset_local(nsize);
//    for (int i=0; i<nsize; i++) {
//        id_local[i] = localNode[i].id;
//        offset_local[i] = localNode[i].mesh_offset;
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//
//
//    auto check_duplicate_list = [](std::vector<int> list, int val)
//    {
//        for (auto& elem : list) { if (elem == val) { return true; } }
//        return  false;
//    };
//
//
//    MPI_Win  win_id_get;
//    MPI_Win_create(id_local.    data(), sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_id_get);
//
//    MPI_Win  win_offset_get;
//    MPI_Win_create(offset_local.data(), sizeof(int)*nsize, sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win_offset_get);
//
//    MPI_Win_fence(0, win_id_get);
//    MPI_Win_fence(0, win_offset_get);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    // @ cal_flag == true //
//    std::vector<int>  lid_dst(nsize*27);
//    std::vector<int>  offset_dst(nsize*27);
//
//    int  icount;
//    icount = 0;
//    for (int i=0; i<nsize; i++) {
//        if ( !localNode[i].cal_flag ) {
//            localNode[i].flag_put = false;
//            continue;
//        }
//        const Node* node = tree_global.nodes(i);
//
//        // remove_rank //
//        std::vector<int>  remove_rank;  remove_rank.push_back(rank);
//        for(int kk=-1; kk<=1; kk++) {
//        for(int jj=-1; jj<=1; jj++) {
//        for(int ii=-1; ii<=1; ii++) {
//            if ( !node->neighbor(ii, jj, kk)->nodeCalFlags().Cal() ) { continue; }
//
//            const int  id_neighbor = node->neighbor(ii, jj, kk)->index();
//            const int  rank_dst = localNode[id_neighbor].cal_rank;
//
//            if ( check_duplicate_list(remove_rank, rank_dst) ) { continue; }
//            remove_rank.push_back(rank_dst);
//
//            // MPI_Get //
//            MPI_Get(&lid_dst   [icount], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_id_get);
//            MPI_Get(&offset_dst[icount], 1, MPI_INT, rank_dst, i, 1, MPI_INT, win_offset_get);
//            icount++;
//        }
//        }
//        }
//    }
//    MPI_Win_fence(0, win_id_get);
//    MPI_Win_fence(0, win_offset_get);
//
//    MPI_Win_free(&win_id_get);
//    MPI_Win_free(&win_offset_get);
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    icount = 0;
//    for (int i=0; i<nsize; i++) {
//        if ( !localNode[i].cal_flag ) {
//            localNode[i].flag_put = false;
//            continue;
//        }
//        const Node* node = tree_global.nodes(i);
//
//        // remove_rank //
//        std::vector<int>  remove_rank;  remove_rank.push_back(rank);
//        for(int kk=-1; kk<=1; kk++) {
//        for(int jj=-1; jj<=1; jj++) {
//        for(int ii=-1; ii<=1; ii++) {
//            if ( !node->neighbor(ii, jj, kk)->nodeCalFlags().Cal() ) { continue; }
//
//            const int  id_neighbor = node->neighbor(ii, jj, kk)->index();
//            const int  rank_dst = localNode[id_neighbor].cal_rank;
//
//            if ( check_duplicate_list(remove_rank, rank_dst) ) { continue; }
//            remove_rank.push_back(rank_dst);
//
//            // update //
//            localNode[i].flag_put = true;
//            localNode[i].id_put.    push_back( lid_dst[icount] );
//            localNode[i].offset_put.push_back( offset_dst[icount] );
//            localNode[i].rank_put.  push_back( rank_dst );
//
//            icount++;
//        }
//        }
//        }
//    }
//
//}


void  InitTree_for_MPI::
commit_MPI_struct(
          MPI_Datatype&            mpi_Datatype,
          void*&                   base_ptr,
          std::vector<int>&        val,
    const std::vector<_MPIGetInfo>& mpiGetInfos,
    const int                      istart,
    const int                      count
    )
{
    MPI_Aint  addr[count];
    for (int i=0; i<count; i++) {
        const int j = istart + i;

        const int gindex = mpiGetInfos[j].gindex;
        MPI_Get_address( &val.data()[gindex], &addr[i] );
    }
    MPI_Aint  base = addr[0];

    MPI_Aint  displ[count];
    for (int i=0; i<count; i++) {
        displ[i] = addr[i] - base;
    }

    MPI_Datatype dtype[count];
    int  blength[count];
    for (int i=0; i<count; i++) {
        dtype  [i] = MPI_INT; // local index //
        blength[i] = 1;
    }

    // submit //
    base_ptr = reinterpret_cast<void*>( base );
    MPI_Type_create_struct( count, blength, displ, dtype, &mpi_Datatype );
    MPI_Type_commit( &mpi_Datatype );
}
