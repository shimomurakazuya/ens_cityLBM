#pragma once
#ifndef TASKID_H_
#define TASKID_H_


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <array>

#include "managed_vector.hpp"
#include "defineAMR.h"
#include "defineFilenames.h"
#include "Tree.h"
#include "MeshValue.h"

#include "FuncLoop.h"
#include "FuncMath.h"
#include "IOData.h"
#include "Index.h"
#include "IndexLBM.h"
#include "FuncObj.h"


class  TaskID {
public:
    using vector_type = util::managed_vector<int>;

private:
    std::array < vector_type, DefAMR::LV_MAX >
        id_tasks_all_,      
        id_tasks_lbm_,      
        id_tasks_L2F_,      
        id_tasks_F2L_,      
        id_tasks_with_halo_, // tb //
        // optimized for cal //
        // wo halo //
        id_tasks_lbm_w_bc_,
        id_tasks_lbm_wo_bc_,
        // tb //
        id_tasks_with_halo_w_bc_,
        id_tasks_with_halo_wo_bc_;

public:
    TaskID () {}

    ~TaskID () {}

public: // accessors //
    const auto& id_tasks_all      (int lv) const { return  id_tasks_all_      [lv]; }
    const auto& id_tasks_lbm      (int lv) const { return  id_tasks_lbm_      [lv]; }
    const auto& id_tasks_L2F      (int lv) const { return  id_tasks_L2F_      [lv]; }
    const auto& id_tasks_F2L      (int lv) const { return  id_tasks_F2L_      [lv]; }
    const auto& id_tasks_with_halo(int lv) const { return  id_tasks_with_halo_[lv]; }

    const auto* id_task_array_all      (int lv) const { return id_tasks_all      (lv).data();  }
    const auto* id_task_array_lbm      (int lv) const { return id_tasks_lbm      (lv).data();  }
    const auto* id_task_array_L2F      (int lv) const { return id_tasks_L2F      (lv).data();  }
    const auto* id_task_array_F2L      (int lv) const { return id_tasks_F2L      (lv).data();  }
    const auto* id_task_array_with_halo(int lv) const { return id_tasks_with_halo(lv).data();  }

    auto num_task_array_all      (int lv) const { return  id_tasks_all      (lv).size(); }
    auto num_task_array_lbm      (int lv) const { return  id_tasks_lbm      (lv).size(); }
    auto num_task_array_L2F      (int lv) const { return  id_tasks_L2F      (lv).size(); }
    auto num_task_array_F2L      (int lv) const { return  id_tasks_F2L      (lv).size(); }
    auto num_task_array_with_halo(int lv) const { return  id_tasks_with_halo(lv).size(); }

    // opt //
    // without halo //
    const auto& id_tasks_lbm_w_bc (int lv) const { return  id_tasks_lbm_w_bc_ [lv]; }
    const auto& id_tasks_lbm_wo_bc(int lv) const { return  id_tasks_lbm_wo_bc_[lv]; }

    const auto* id_task_array_lbm_w_bc (int lv) const { return  id_tasks_lbm_w_bc (lv).data(); }
    const auto* id_task_array_lbm_wo_bc(int lv) const { return  id_tasks_lbm_wo_bc(lv).data(); }

    auto num_task_array_lbm_w_bc (int lv) const { return  id_tasks_lbm_w_bc (lv).size(); }
    auto num_task_array_lbm_wo_bc(int lv) const { return  id_tasks_lbm_wo_bc(lv).size(); }

    // with halo //
    const auto& id_tasks_with_halo_w_bc (int lv) const { return  id_tasks_with_halo_w_bc_ [lv]; }
    const auto& id_tasks_with_halo_wo_bc(int lv) const { return  id_tasks_with_halo_wo_bc_[lv]; }

    const auto* id_task_array_with_halo_w_bc (int lv) const { return  id_tasks_with_halo_w_bc (lv).data(); }
    const auto* id_task_array_with_halo_wo_bc(int lv) const { return  id_tasks_with_halo_wo_bc(lv).data(); }

    auto  num_task_array_with_halo_w_bc (int lv) const { return  id_tasks_with_halo_w_bc (lv).size(); }
    auto  num_task_array_with_halo_wo_bc(int lv) const { return  id_tasks_with_halo_wo_bc(lv).size(); }

    // mlups //
    int  num_task_lbm_per_step() const {
        int  tmp = 0;
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            tmp += num_task_array_lbm(lv) * pow(2, lv);
        }
        return tmp;
    }

public:
    void  set_taskID(const Tree&  tree) {
        clear_taskID();
        init_taskID_vector(tree);
    }

    void  set_taskID_opt(const Tree&  tree, const MeshValue* meshValues) {
        clear_taskID_opt();
        init_taskID_vector_opt(tree, meshValues);
    }

    void  writeTaskID(const std::string filename) const {
        // write //
        std::ofstream  fout;
        fout.open(filename, std::ios::binary);

        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            write_vector( fout, id_tasks_all_[lv] );
            write_vector( fout, id_tasks_lbm_[lv] );
            write_vector( fout, id_tasks_L2F_[lv] );
            write_vector( fout, id_tasks_F2L_[lv] );

            write_vector( fout, id_tasks_with_halo_[lv] );

            // without halo //
            write_vector( fout, id_tasks_lbm_w_bc_ [lv] );
            write_vector( fout, id_tasks_lbm_wo_bc_[lv] );

            // with halo //
            write_vector( fout, id_tasks_with_halo_w_bc_ [lv] );
            write_vector( fout, id_tasks_with_halo_wo_bc_[lv] );
        }

        fout.close();
    }

    void  readTaskID(const std::string filename) {
        clear_taskID();
        clear_taskID_opt();

        // write //
        std::ifstream  fin;
        fin.open(filename, std::ios::binary);

        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            read_vector( fin, id_tasks_all_[lv] );
            read_vector( fin, id_tasks_lbm_[lv] );
            read_vector( fin, id_tasks_L2F_[lv] );
            read_vector( fin, id_tasks_F2L_[lv] );

            read_vector( fin, id_tasks_with_halo_[lv] );

            // without halo //
            read_vector( fin, id_tasks_lbm_w_bc_ [lv] );
            read_vector( fin, id_tasks_lbm_wo_bc_[lv] );

            // with halo //
            read_vector( fin, id_tasks_with_halo_w_bc_ [lv] );
            read_vector( fin, id_tasks_with_halo_wo_bc_[lv] );
        }

        fin.close();

    }

    void  copyTaskID(const TaskID& taskID, bool  is_reset) {
        if (is_reset) {
            clear_taskID();
            clear_taskID_opt();
        }

        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            id_tasks_all_      [lv] = taskID.id_tasks_all_      [lv];
            id_tasks_lbm_      [lv] = taskID.id_tasks_lbm_      [lv];
            id_tasks_L2F_      [lv] = taskID.id_tasks_L2F_      [lv];
            id_tasks_F2L_      [lv] = taskID.id_tasks_F2L_      [lv];
            id_tasks_with_halo_[lv] = taskID.id_tasks_with_halo_[lv];

            // without halo //
            id_tasks_lbm_w_bc_ [lv] = taskID.id_tasks_lbm_w_bc_ [lv];
            id_tasks_lbm_wo_bc_[lv] = taskID.id_tasks_lbm_wo_bc_[lv];

            // with halo //
            id_tasks_with_halo_w_bc_ [lv] = taskID.id_tasks_with_halo_w_bc_ [lv];
            id_tasks_with_halo_wo_bc_[lv] = taskID.id_tasks_with_halo_wo_bc_[lv];
        }
    }

    std::string filename(const int rank, const int step) const {
        return    Foldernames::io_folder + "/"
                + Filenames::taskID_name0
                + "-rank" + std::to_string(rank)
                + "-step" + std::to_string(step)
                + ".dat";
    }

private:
    void  clear_taskID() {
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            id_tasks_all_      [lv].clear();
            id_tasks_lbm_      [lv].clear();
            id_tasks_L2F_      [lv].clear();
            id_tasks_F2L_      [lv].clear();
            id_tasks_with_halo_[lv].clear();

            id_tasks_all_      [lv].shrink_to_fit();
            id_tasks_lbm_      [lv].shrink_to_fit();
            id_tasks_L2F_      [lv].shrink_to_fit();
            id_tasks_F2L_      [lv].shrink_to_fit();
            id_tasks_with_halo_[lv].shrink_to_fit();

        }
    }

    void  clear_taskID_opt() {
        // without halo //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            id_tasks_lbm_w_bc_ [lv].clear();
            id_tasks_lbm_wo_bc_[lv].clear();

            id_tasks_lbm_w_bc_ [lv].shrink_to_fit();
            id_tasks_lbm_wo_bc_[lv].shrink_to_fit();
        }


        // with halo //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            id_tasks_with_halo_w_bc_ [lv].clear();
            id_tasks_with_halo_wo_bc_[lv].clear();

            id_tasks_with_halo_w_bc_ [lv].shrink_to_fit();
            id_tasks_with_halo_wo_bc_[lv].shrink_to_fit();
        }
    }

    void  init_taskID_vector(const Tree& tree) {
        const int  n_leaf = tree.number_of_nodes();

        // id_tasks //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            for (int i=0; i<n_leaf; i++) {
                const Node*  node = tree.nodes(i);

                if ( node->level() == lv ) {
                    id_tasks_all_[lv].push_back(i);
                    if ( node->nodeCalFlags().Cal() ) { id_tasks_lbm_[lv].push_back(i); }

                    // L2X //
                    if ( node->nodeCalFlags().L2F() ) { id_tasks_L2F_[lv].push_back(i); }

                    // X2L //
                    if ( node->nodeCalFlags().F2L() ) {
                        for (int kk=0; kk<=1; kk++) {
                        for (int jj=0; jj<=1; jj++) {
                        for (int ii=0; ii<=1; ii++) {
                            if ( node->child()->neighbor(ii,jj,kk)->nodeCalFlags().L2C() ) {
                                id_tasks_F2L_[lv].push_back( node->child()->neighbor(ii,jj,kk)->index() );
                            }
                        }
                        }
                        }
                    }

                    // with halo //
                    bool flag_with_halo = false;
                    for (int kk=-1; kk<=1; kk++) {
                    for (int jj=-1; jj<=1; jj++) {
                    for (int ii=-1; ii<=1; ii++) {
                        if ( node->neighbor(ii,jj,kk)->nodeCalFlags().Cal() ) {
                            flag_with_halo = true;
                            ii = 2; jj = 2; kk = 2;
                        }
                    }
                    }
                    }
                    if (flag_with_halo) {
                        id_tasks_with_halo_[lv].push_back(i);
                    }
                }

            }
        }
    } // void init_taskID_vector()

    void  init_taskID_vector_opt(const Tree&  tree, const MeshValue*  meshValue) {
        const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

        // without halo //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            const int*   id_tasks = id_task_array_lbm(lv);
            const real*  lv_obj   = meshValue[lv].valueObjLS().lv_obj();

            for (int l=0; l<num_task_array_lbm(lv); l++) {
                bool has_object = leaf_has_object( id_tasks[l], mesh_offsets3x3x3, lv_obj );

                if (has_object) { id_tasks_lbm_w_bc_ [lv].push_back( id_tasks[l] ); }
                else            { id_tasks_lbm_wo_bc_[lv].push_back( id_tasks[l] ); }
            }
        }

        // with halo //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            const int*   id_tasks = id_task_array_with_halo(lv);
            const real*  lv_obj   = meshValue[lv].valueObjLS().lv_obj();

            for (int l=0; l<num_task_array_with_halo(lv); l++) {
                bool has_object = leaf_has_object( id_tasks[l], mesh_offsets3x3x3, lv_obj );

                if (has_object) { id_tasks_with_halo_w_bc_ [lv].push_back( id_tasks[l] ); }
                else            { id_tasks_with_halo_wo_bc_[lv].push_back( id_tasks[l] ); }
            }
        }
    } // void init_taskID_vector_opt()

    bool leaf_has_object(
        const int    idl,
        const int*   mesh_offsets3x3x3,
        const real*  lv_obj
    ) {
        constexpr int  nQ = LBM_velocity_model::nQ;

        int  offset3d[27];
        for (int idv=0; idv<nQ; idv++) {
            offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
        }

        for (int k=0; k<DefAMR::NX_LEAF; k++) {
        for (int j=0; j<DefAMR::NX_LEAF; j++) {
        for (int i=0; i<DefAMR::NX_LEAF; i++) {
            int  ids[27];
            for (int kv=-1; kv<=1; kv++) {
                for (int jv=-1; jv<=1; jv++) {
                    for (int iv=-1; iv<=1; iv++) {
                        const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                        ids[idv_leaf] = Index::id(i+iv,j+jv,k+kv, offset3d);
                    }
                }
            }


            for (int kv=-1; kv<=1; kv++) {
                for (int jv=-1; jv<=1; jv++) {
                    for (int iv=-1; iv<=1; iv++) {
                        const int  id_up = Index::id(ids, -iv,-jv,-kv);

                        if ( FuncObj::is_obj( lv_obj[id_up] ) )  { return true; }
                    }
                }
            }
        }
        }
        }
        return  false;
    } // bool leaf_has_object()

    void  write_vector(std::ofstream& fout, const vector_type& vec) const {
        const int count = vec.size();
        runtime_assert(fout.write( ( char * ) &count, sizeof( int ) ), "FileIOError");
        runtime_assert(fout.write( ( char * ) vec.data(), sizeof( int ) * count ), "FileIOError");
    }

    void  read_vector(std::ifstream& fin, vector_type& vec) {
        int count = -1;
        runtime_assert(fin.read( ( char * ) &count, sizeof( int ) ), "FileIOError");
        vec.resize(count);
        runtime_assert(fin.read( ( char * ) vec.data(), sizeof( int ) * count ), "FileIOError");
    }

};


#endif
