#include "InitTree.h"
#include "defineCal.h"
#include "FuncLoop.h"
#include "Vector3d.h"
#include "InitMonitorLevelset.h"
#include "mpi_wrapper.hpp"


namespace  InitTree {


void
init_tree_uniform3d(
          Tree& tree,
    const Grid* grids,
    const int   lv_max
    )
{
    Tree tree_all = tree_all_region(tree.comm(), grids, lv_max);

    const bool  is_reset = true;
    tree.copy_node_structures(tree_all, is_reset);
}


//void
//init_tree_amr(
//          Tree& tree,
//    const Grid* grids,
//    const int   lv_max
//    )
//{
////    std::cout << __PRETTY_FUNCTION__ << std::endl;
//
//    Tree tree_all    = tree_all_region(grids, lv_max);
//
//    // refinement condition //
//    InitMonitorLevelset  initMonitorLevelset;
//    std::vector<int>  cal_flags = initMonitorLevelset.create_id_flags(grids);
//    MPI_Barrier(MPI_COMM_WORLD);
//    // refinement condition //
//
//
//    std::vector<id_amr_init>  id_init = create_id_amr_init(tree, grids, lv_max, cal_flags);
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    Tree tree_new;
//    update_basic_informations(tree_new, tree_all, grids, lv_max, id_init); MPI_Barrier(MPI_COMM_WORLD);
//    update_connections       (tree_new, tree_all, grids, lv_max, id_init); MPI_Barrier(MPI_COMM_WORLD);
//    update_NodeCalFlags      (tree_new);
//
//    // copy //
//    const bool  is_reset = true;
//    tree.copy_node_structures(tree_new, is_reset);
//}

void
init_tree_amr_with_map(
          Tree& tree,
    const Grid* grids,
    const int   lv_max,
    const MapData& map
    )
{
    const bool output_flag = tree.comm().is_rank0();
    if(output_flag) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    Tree tree_all    = tree_all_region(tree.comm(), grids, lv_max);

    // refinement condition //
    InitMonitorLevelset  initMonitorLevelset;
    std::vector<int>  cal_flags = initMonitorLevelset.create_id_flags_with_map(grids, map);
    MPI_Barrier(MPI_COMM_WORLD);
    // refinement condition //

    MPI_Barrier(MPI_COMM_WORLD);
    if (output_flag){ std::cout << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; }


    std::vector<id_amr_init>  id_init = create_id_amr_init(tree, grids, lv_max, cal_flags);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; }

    auto&& tree_new = tree.fork();
    update_basic_informations(tree_new, tree_all, grids, lv_max, id_init); MPI_Barrier(MPI_COMM_WORLD);
    update_connections       (tree_new, tree_all, grids, lv_max, id_init); MPI_Barrier(MPI_COMM_WORLD);
    update_NodeCalFlags      (tree_new);

    MPI_Barrier(MPI_COMM_WORLD);
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; }


    // copy //
    const bool  is_reset = true;
    tree.copy_node_structures(tree_new, is_reset);

    if(output_flag) std::cout << __PRETTY_FUNCTION__ << " finished" << std::endl;
}


Tree
tree_all_region(
    const MPICommEnsemble comm,
    const Grid* grids,
    const int   lv_max
    )
{
    MPI_Barrier(MPI_COMM_WORLD);
    const bool output_flag = comm.is_rank0();
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    int  nx[lv_max];
    int  ny[lv_max];
    int  nz[lv_max];
    for (int i=0; i<lv_max; i++) {
        nx[i] = grids[i].nx().x();
        ny[i] = grids[i].nx().y();
        nz[i] = grids[i].nx().z();
    }

    Tree  tree(comm);

    // allocate //
    int  offset = 0;
    for (int lv=0; lv<lv_max; lv++) {
        FuncLoop::Loop1d  loop1d_all(0, nx[lv]*ny[lv]*nz[lv]);
        loop1d_all.for_each( [&offset, &tree](const int i){ tree.add_node(offset+i); } );

        offset += (nx[lv]*ny[lv]*nz[lv]);
    }

    // informations //
    int  id_lv[lv_max];
    id_lv[0] = 0;
    for (int lv=0; lv<lv_max-1; lv++) { id_lv[lv+1] = id_lv[lv] + (nx[lv]*ny[lv]*nz[lv]); }

    for (int lv=0; lv<lv_max; lv++) {

        int  mesh_offset = 0;
        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
                    const int  id  = i + nx[lv]*j + nx[lv]*ny[lv]*k;
                    const int  idg = id + id_lv[lv];

                    // this //
                    tree.nodes(idg)->set_level(lv);
                    tree.nodes(idg)->set_node_type(NodeTypes::LeafT);
                    tree.nodes(idg)->set_position(lv, i,j,k);

                    // nodeCalFlags //
                    if (lv == lv_max-1) { tree.nodes(idg)->nodeCalFlags().set_flag_Cal(true);  }
                    else                { tree.nodes(idg)->nodeCalFlags().set_flag_Cal(false); }


                    // mesh_offset //
                    tree.nodes(idg)->set_mesh_offset(mesh_offset);
                    mesh_offset += pow(DefAMR::NX_LEAF, 3);


                    // parent //
                    if (lv == 0) {
                        tree.nodes(idg)->set_parent(tree.root_node());
                    }
                    else {
                        const int  i_p = i/2;
                        const int  j_p = j/2;
                        const int  k_p = k/2;
                        const int  id_p = id_lv[lv-1] + i_p + nx[lv-1]*j_p + nx[lv-1]*ny[lv-1]*k_p;

                        tree.nodes(idg)->set_parent(tree.nodes(id_p)->node());
                    }


                    // child //
                    if (lv == lv_max-1) {
                        tree.nodes(idg)->set_child(tree.external_leaf_node());
                    }
                    else {
                        const int  i_c = i*2;
                        const int  j_c = j*2;
                        const int  k_c = k*2;
                        const int  id_c = id_lv[lv+1] + i_c + nx[lv+1]*j_c + nx[lv+1]*ny[lv+1]*k_c;

                        tree.nodes(idg)->set_child(tree.nodes(id_c)->node());
                    }


                    // neighbor //
                    // default : periodic boundary //
                    const int  idw = ((i-1+nx[lv])%nx[lv]) + nx[lv]*j + nx[lv]*ny[lv]*k + id_lv[lv];
                    const int  ide = ((i+1+nx[lv])%nx[lv]) + nx[lv]*j + nx[lv]*ny[lv]*k + id_lv[lv];

                    const int  ids = i + nx[lv]*((j-1+ny[lv])%ny[lv]) + nx[lv]*ny[lv]*k + id_lv[lv];
                    const int  idn = i + nx[lv]*((j+1+ny[lv])%ny[lv]) + nx[lv]*ny[lv]*k + id_lv[lv];

                    const int  idb = i + nx[lv]*j + nx[lv]*ny[lv]*((k-1+nz[lv])%nz[lv]) + id_lv[lv];
                    const int  idt = i + nx[lv]*j + nx[lv]*ny[lv]*((k+1+nz[lv])%nz[lv]) + id_lv[lv];

                    tree.nodes(idg)->set_neighbor(tree.nodes(idw)->node(), -1,  0,  0);
                    tree.nodes(idg)->set_neighbor(tree.nodes(ide)->node(),  1,  0,  0);

                    tree.nodes(idg)->set_neighbor(tree.nodes(ids)->node(),  0, -1,  0);
                    tree.nodes(idg)->set_neighbor(tree.nodes(idn)->node(),  0,  1,  0);

                    tree.nodes(idg)->set_neighbor(tree.nodes(idb)->node(),  0,  0, -1);
                    tree.nodes(idg)->set_neighbor(tree.nodes(idt)->node(),  0,  0,  1);
                    if ( !OuterBoundaryConditions::OuterPeriodicCommZ ) { // defineCal //
//                        tree.nodes(idg)->set_neighbor(tree.nodes(idb)->node(),  0,  0, 0); // error
//                        tree.nodes(idg)->set_neighbor(tree.nodes(idt)->node(),  0,  0, 0); // error
                        if ( k == 0        ) { tree.nodes(idg)->set_neighbor(tree.external_neighbor_node(),  0,  0, -1); }
                        if ( k == nz[lv]-1 ) { tree.nodes(idg)->set_neighbor(tree.external_neighbor_node(),  0,  0,  1); }
                    }

                } // i //
            } // j //
        } // k //

    } // lv //

    MPI_Barrier(MPI_COMM_WORLD);
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << ": " << __LINE__ << std::endl; }

    return  tree;
}


std::vector<id_amr_init>
create_id_amr_init(
    const Tree& tree_org,
    const Grid* grids,
    const int   lv_max,
    const std::vector<int> cal_flags
    )
{
    const bool output_flag = tree_org.comm().is_rank0();
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    std::vector<id_amr_init>  id_init;

    int  nx[lv_max];
    int  ny[lv_max];
    int  nz[lv_max];
    for (int i=0; i<lv_max; i++) {
        nx[i] = grids[i].nx().x();
        ny[i] = grids[i].nx().y();
        nz[i] = grids[i].nx().z();
    }


    // index //
    const int  cal_true  = 1;
    const int  cal_false = 0;

    int  offset = 0;
    for (int lv=0; lv<lv_max; lv++) {
        FuncLoop::Loop1d  loop1d_all(0, nx[lv]*ny[lv]*nz[lv]);
        loop1d_all.for_each( [&offset, &id_init, cal_true, cal_false](const int i){ id_init.push_back( { offset+i, cal_true } ); } );

        offset += (nx[lv]*ny[lv]*nz[lv]);
    }


    int  id_tmp;
    int  id_filtered;

    // calculate flag //
    if (output_flag) { std::cout << "calculate flag" << std::endl; }
    id_tmp      = 0;
    for (int lv=0; lv<lv_max; lv++) {
        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
                    id_init[id_tmp].flag_cal = cal_flags[id_tmp];

                    // update indices //
                    id_tmp++;
                }
            }
        }
    }

    // allocate flag //
    if (output_flag) { std::cout << "allocate flag" << std::endl; }
    id_tmp      = 0;
    for (int lv=0; lv<lv_max; lv++) {
        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
                    int  imin = (i%2 == 1) ? -2 : -1;
                    int  jmin = (j%2 == 1) ? -2 : -1;
                    int  kmin = (k%2 == 1) ? -2 : -1;

                    int  imax = (i%2 == 0) ?  2 :  1;
                    int  jmax = (j%2 == 0) ?  2 :  1;
                    int  kmax = (k%2 == 0) ?  2 :  1;

                    bool  allocate_flag = check_neighbor(id_tmp, i,j,k, imin,jmin,kmin, imax,jmax,kmax, nx[lv],ny[lv],nz[lv], id_init);
//                    bool  cal_flag = id_init[id_tmp].flag_cal;

                    id_init[id_tmp].flag_allocate = allocate_flag;

                    // update indices //
                    id_tmp++;
                }
            }
        }
    }



    if (output_flag) { std::cout << "id_init" << std::endl; }
    id_tmp      = 0;
    id_filtered = 0;
    for (int lv=0; lv<lv_max; lv++) {
        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
                    id_init[id_tmp].id  = id_filtered;
                    id_init[id_tmp].lv  = lv;

                    // update indices //
                    if ( id_init[id_tmp].flag_allocate ) { id_filtered++; }
                    id_tmp++;
                }
            }
        }
    }

#if 0
    // check //
    id_tmp      = 0;
    id_filtered = 0;
    for (int lv=0; lv<lv_max; lv++) {
        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
//                    if (id_init[id_tmp].flag_allocate) {
                    if (id_init[id_tmp].flag_cal) {
                        std::cout << " o " <<  " ";
//                        std::cout << id_init[id_tmp].id <<  " ";
                    }
                    else {
                        std::cout << "---" << " ";
                    }
//                    std::cout << id_init[id_tmp].flag_cal << "," << id_init[id_tmp].flag_allocate << " ";
                    // update indices //
                    id_tmp++;
                    if ( id_init[id_tmp].flag_allocate ) { id_filtered++; }
                }
                    std::cout << std::endl;
            }
                    std::cout << std::endl;
        }
                    std::cout << std::endl;
    }
#endif

    return id_init;
}


void
update_basic_informations(
          Tree& tree,
    const Tree& tree_all,
    const Grid* grids,
    const int   lv_max,
    const std::vector<id_amr_init>&  id_input
    )
{
    const bool output_flag = tree.comm().is_rank0();
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    int  nx[lv_max];
    int  ny[lv_max];
    int  nz[lv_max];
    for (int i=0; i<lv_max; i++) {
        nx[i] = grids[i].nx().x();
        ny[i] = grids[i].nx().y();
        nz[i] = grids[i].nx().z();
    }


    // index //
    int  id_tmp;
    int  id_filtered;

    id_tmp      = 0;
    id_filtered = 0;
    for (int lv=0; lv<lv_max; lv++) {
        int  mesh_offset = 0;

        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
                    if ( (id_input[id_tmp].flag_allocate) ) {
                        // index //
                        tree.add_node(id_filtered);

                        // this //
                        // position //
                        tree.nodes(id_filtered)->set_level(lv);
                        tree.nodes(id_filtered)->set_node_type(NodeTypes::LeafT);
                        tree.nodes(id_filtered)->set_position(lv, i,j,k);

                        // periodic boundary //
                        tree.nodes(id_filtered)->nodeCalFlags().reset();
                        tree.nodes(id_filtered)->nodeCalFlags().set_flag_Cal( id_input[id_tmp].flag_cal );

                        // mesh_offset //
                        tree.nodes(id_filtered)->set_mesh_offset(mesh_offset);
                        mesh_offset += pow(DefAMR::NX_LEAF, 3);
                    }

                    // increment //
                    if ( (id_input[id_tmp].flag_allocate) ) { id_filtered++; }
                    id_tmp++;
                }
            }
        }
    }

    if (output_flag) { std::cout << "set_num_nodes_lv" << std::endl; }
    tree.set_num_nodes_lv();

//    for (int lv=0; lv<lv_max; lv++) {
//        std::cout << __PRETTY_FUNCTION__ << " : lv, number_of_nodes_lv " << lv << ", " << tree.number_of_nodes_lv(lv) << std::endl;
//    }
}


void
update_connections(
          Tree& tree,
    const Tree& tree_all,
    const Grid* grids,
    const int   lv_max,
    const std::vector<id_amr_init>&  id_input
    )
{
    const bool output_flag = tree.comm().is_rank0();
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    int  nx[lv_max];
    int  ny[lv_max];
    int  nz[lv_max];
    for (int i=0; i<lv_max; i++) {
        nx[i] = grids[i].nx().x();
        ny[i] = grids[i].nx().y();
        nz[i] = grids[i].nx().z();
    }


    // index //
    int  id_tmp;
    int  id_filtered;

    // parent, child, neighbor //
    id_tmp      = 0;
    id_filtered = 0;
    for (int lv=0; lv<lv_max; lv++) {

        for (int k=0; k<nz[lv]; k++) {
            for (int j=0; j<ny[lv]; j++) {
                for (int i=0; i<nx[lv]; i++) {
                    if ( id_input[id_tmp].flag_allocate ) {
                        // parent //
                        if ( tree_all.nodes(id_tmp)->parent()->node_type() == NodeTypes::LeafT ) {
                            int  _id_p = tree_all.nodes(id_tmp)->parent()->index();
                            if ( id_input[_id_p].flag_allocate ) {
                                int  id_p = id_input[_id_p].id;
                                tree.nodes(id_filtered)->set_parent( tree.nodes(id_p)->node() );
                            }
                            else {
                                tree.nodes(id_filtered)->set_parent( tree.external_node( NodeTypes::RootT ) );
                            }
//                            else { std::cout << __PRETTY_FUNCTION__ << " : parent : " << "error" << std::endl;  exit(-1); } // comment out ? //

                        }
                        else { tree.nodes(id_filtered)->set_parent( tree.external_node( tree_all.nodes(id_tmp)->parent()->node_type() ) ); }
//                        std::cout << lv << " : tree : parent = " << tree.nodes(id_filtered)->parent()->node_type() << std::endl;


                        // child //
                        if ( tree_all.nodes(id_tmp)->child()->node_type() == NodeTypes::LeafT ) {
                            int  _id_c = tree_all.nodes(id_tmp)->child()->index();
                            if ( id_input[_id_c].flag_allocate ) {
                                int  id_c = id_input[_id_c].id;
                                tree.nodes(id_filtered)->set_child( tree.nodes(id_c)->node() );
                            }
                            else {
                                tree.nodes(id_filtered)->set_child( tree.external_node( NodeTypes::ExternalLeafT ) );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_child( tree.external_node( tree_all.nodes(id_tmp)->child()->node_type() ) ); }
//                        std::cout << "tree : child = " << tree.nodes(id_filtered)->child()->node_type() << std::endl;


                        // neighbor //
                        // default : periodic boundary //
                        if ( tree_all.nodes(id_tmp)->neighbor(-1, 0, 0)->node_type() == NodeTypes::LeafT ) {
                            int _idw = tree_all.nodes(id_tmp)->neighbor(-1, 0, 0)->index();
                            if ( id_input[_idw].flag_allocate ) {
                                int idw = id_input[_idw].id;
                                tree.nodes(id_filtered)->set_neighbor( tree.nodes(idw)->node(), -1, 0, 0 );
                            }
                            else {
                                tree.nodes(id_filtered)->set_neighbor( tree.external_node( NodeTypes::ExternalNeighborT ), -1, 0, 0 );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_neighbor( tree.external_node( tree_all.nodes(id_tmp)->neighbor(-1, 0, 0)->node_type() ), -1, 0, 0 ); }


                        if ( tree_all.nodes(id_tmp)->neighbor(1, 0, 0)->node_type() == NodeTypes::LeafT ) {
                            int _ide = tree_all.nodes(id_tmp)->neighbor(1, 0, 0)->index();
                            if ( id_input[_ide].flag_allocate ) {
                                int ide = id_input[_ide].id;
                                tree.nodes(id_filtered)->set_neighbor( tree.nodes(ide)->node(),  1, 0, 0 );
                            }
                            else {
                                tree.nodes(id_filtered)->set_neighbor( tree.external_node( NodeTypes::ExternalNeighborT ), 1, 0, 0 );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_neighbor( tree.external_node( tree_all.nodes(id_tmp)->neighbor(1, 0, 0)->node_type() ),  1, 0, 0 ); }


                        if ( tree_all.nodes(id_tmp)->neighbor(0, -1, 0)->node_type() == NodeTypes::LeafT ) {
                            int _ids = tree_all.nodes(id_tmp)->neighbor(0, -1, 0)->index();
                            if ( id_input[_ids].flag_allocate ) {
                                int ids = id_input[_ids].id;
                                tree.nodes(id_filtered)->set_neighbor( tree.nodes(ids)->node(), 0, -1, 0 );
                            }
                            else {
                                tree.nodes(id_filtered)->set_neighbor( tree.external_node( NodeTypes::ExternalNeighborT ), 0, -1, 0 );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_neighbor( tree.external_node( tree_all.nodes(id_tmp)->neighbor(0, -1, 0)->node_type() ),  0, -1, 0 ); }


                        if ( tree_all.nodes(id_tmp)->neighbor(0, 1, 0)->node_type() == NodeTypes::LeafT ) {
                            int _idn = tree_all.nodes(id_tmp)->neighbor(0, 1, 0)->index();
                            if ( id_input[_idn].flag_allocate ) {
                                int idn = id_input[_idn].id;
                                tree.nodes(id_filtered)->set_neighbor( tree.nodes(idn)->node(),  0, 1, 0 );
                            }
                            else {
                                tree.nodes(id_filtered)->set_neighbor( tree.external_node( NodeTypes::ExternalNeighborT ), 0, 1, 0 );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_neighbor( tree.external_node( tree_all.nodes(id_tmp)->neighbor(0, 1, 0)->node_type() ),  0, 1, 0 ); }


                        if ( tree_all.nodes(id_tmp)->neighbor(0, 0, -1)->node_type() == NodeTypes::LeafT ) {
                            int _idb = tree_all.nodes(id_tmp)->neighbor(0, 0, -1)->index();
                            if ( id_input[_idb].flag_allocate ) {
                                int idb = id_input[_idb].id;
                                tree.nodes(id_filtered)->set_neighbor( tree.nodes(idb)->node(), 0, 0, -1 );
                            }
                            else {
                                tree.nodes(id_filtered)->set_neighbor( tree.external_node( NodeTypes::ExternalNeighborT ), 0, 0, -1 );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_neighbor( tree.external_node( tree_all.nodes(id_tmp)->neighbor(0, 0, -1)->node_type() ),  0, 0, -1 ); }


                        if ( tree_all.nodes(id_tmp)->neighbor(0, 0, 1)->node_type() == NodeTypes::LeafT ) {
                            int _idt = tree_all.nodes(id_tmp)->neighbor(0, 0, 1)->index();
                            if ( id_input[_idt].flag_allocate ) {
                                int idt = id_input[_idt].id;
                                tree.nodes(id_filtered)->set_neighbor( tree.nodes(idt)->node(),  0, 0, 1 );
                            }
                            else {
                                tree.nodes(id_filtered)->set_neighbor( tree.external_node( NodeTypes::ExternalNeighborT ), 0, 0, 1 );
                            }
                        }
                        else { tree.nodes(id_filtered)->set_neighbor( tree.external_node( tree_all.nodes(id_tmp)->neighbor(0, 0, 1)->node_type() ),  0, 0, 1 ); }

                    }

                    // increment //
                    if ( id_input[id_tmp].flag_allocate ) { id_filtered++; }
                    id_tmp++;
                }
            }
        }
    }
}


void update_NodeCalFlags(Tree& tree)
{
    const bool output_flag = tree.comm().is_rank0();
    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();

    for (int i=0; i<n_leaf; i++) {
        Node*  node = tree.nodes(i);

        // put : L2X //
        if ( node->nodeCalFlags().Cal() ) {
            // parent : L2C //
            if ( node->parent()->node_type() == NodeTypes::LeafT ) {
                if ( node->parent()->nodeCalFlags().Cal() ) {
                    std::cout << "error : parent flag is calculate" << std::endl;
                    exit(-1);
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
                    exit(-1);
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

    if (output_flag) { std::cout << __PRETTY_FUNCTION__ << " finished" << std::endl; }
}


bool check_neighbor(
    int id0,
    int i, int j, int k,
    int imin, int jmin, int kmin,
    int imax, int jmax, int kmax,
    int nx, int ny, int nz,
    std::vector<id_amr_init>& ids
    )
{
    auto ii_cast = [](int i, int ii, int nx) {
        return (i+ii < 0   ) ? 0 :
               (i+ii > nx-1) ? 0:
                               ii;
    };


    bool  is_cal = false;

    for (int kk=kmin; kk<=kmax; kk++) {
    for (int jj=jmin; jj<=jmax; jj++) {
    for (int ii=imin; ii<=imax; ii++) {
        const int  strides[] = { 1, nx, nx*ny };


        const int  ix = ii_cast(i, ii, nx);
        const int  jx = ii_cast(j, jj, ny);
        const int  kx = ii_cast(k, kk, nz);

        const int  id = id0 + strides[0]*ix + strides[1]*jx + strides[2]*kx;

        if ( ids[id].flag_cal == true ) {
            is_cal = true;
            return  is_cal;
        }
    }
    }
    }

    return  is_cal;
}


};
