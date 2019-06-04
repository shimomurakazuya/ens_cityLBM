#pragma once
#ifndef TASKID_H_
#define TASKID_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "defineAMR.h"
#include "Tree.h"
#include "MeshValue.h"


class  TaskID {
private:
    MemType  memType_;

    bool  taskID_is_initialized_;
    bool  taskIDopt_is_initialized_;

    std::vector<int>  id_tasks_all_      [DefAMR::LV_MAX];
    std::vector<int>  id_tasks_lbm_      [DefAMR::LV_MAX];
    std::vector<int>  id_tasks_L2F_      [DefAMR::LV_MAX];
    std::vector<int>  id_tasks_F2L_      [DefAMR::LV_MAX];
    std::vector<int>  id_tasks_with_halo_[DefAMR::LV_MAX]; // tb //

    int*  id_task_array_all_      [DefAMR::LV_MAX];
    int*  id_task_array_lbm_      [DefAMR::LV_MAX]; // calculation reagion //
    int*  id_task_array_L2F_      [DefAMR::LV_MAX];
    int*  id_task_array_F2L_      [DefAMR::LV_MAX];
    int*  id_task_array_with_halo_[DefAMR::LV_MAX];

    int   num_task_array_all_      [DefAMR::LV_MAX];
    int   num_task_array_lbm_      [DefAMR::LV_MAX];
    int   num_task_array_L2F_      [DefAMR::LV_MAX];
    int   num_task_array_F2L_      [DefAMR::LV_MAX];
    int   num_task_array_with_halo_[DefAMR::LV_MAX];

    bool  flag_allocated_array_all_      [DefAMR::LV_MAX];
    bool  flag_allocated_array_lbm_      [DefAMR::LV_MAX];
    bool  flag_allocated_array_L2F_      [DefAMR::LV_MAX];
    bool  flag_allocated_array_F2L_      [DefAMR::LV_MAX];
    bool  flag_allocated_array_with_halo_[DefAMR::LV_MAX];

    // optimize //
    // without halo (calculation reagion) //
    std::vector<int>  id_tasks_lbm_w_bc_ [DefAMR::LV_MAX];
    std::vector<int>  id_tasks_lbm_wo_bc_[DefAMR::LV_MAX];

    int*  id_task_array_lbm_w_bc_ [DefAMR::LV_MAX];
    int*  id_task_array_lbm_wo_bc_[DefAMR::LV_MAX];

    int   num_task_array_lbm_w_bc_ [DefAMR::LV_MAX];
    int   num_task_array_lbm_wo_bc_[DefAMR::LV_MAX];

    bool  flag_allocated_array_lbm_w_bc_ [DefAMR::LV_MAX];
    bool  flag_allocated_array_lbm_wo_bc_[DefAMR::LV_MAX];


    // with halo //
    std::vector<int>  id_tasks_with_halo_w_bc_ [DefAMR::LV_MAX];
    std::vector<int>  id_tasks_with_halo_wo_bc_[DefAMR::LV_MAX];

    int*  id_task_array_with_halo_w_bc_ [DefAMR::LV_MAX];
    int*  id_task_array_with_halo_wo_bc_[DefAMR::LV_MAX];

    int   num_task_array_with_halo_w_bc_ [DefAMR::LV_MAX];
    int   num_task_array_with_halo_wo_bc_[DefAMR::LV_MAX];

    bool  flag_allocated_array_with_halo_w_bc_ [DefAMR::LV_MAX];
    bool  flag_allocated_array_with_halo_wo_bc_[DefAMR::LV_MAX];

public:
    TaskID () {
        taskID_is_initialized_    = false;
        taskIDopt_is_initialized_ = false;

        memType_ = MemType::Host;
#ifdef GPU_CALCULATION__
        memType_ = MemType::Managed;
#endif

        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            flag_allocated_array_all_      [lv] = false;
            flag_allocated_array_lbm_      [lv] = false;
            flag_allocated_array_L2F_      [lv] = false;
            flag_allocated_array_F2L_      [lv] = false;
            flag_allocated_array_with_halo_[lv] = false;
        }
    }

    ~TaskID () {}

public:
    const std::vector<int>& id_tasks_all      (int lv) const { return  id_tasks_all_      [lv]; }
    const std::vector<int>& id_tasks_lbm      (int lv) const { return  id_tasks_lbm_      [lv]; }
    const std::vector<int>& id_tasks_L2F      (int lv) const { return  id_tasks_L2F_      [lv]; }
    const std::vector<int>& id_tasks_F2L      (int lv) const { return  id_tasks_F2L_      [lv]; }
    const std::vector<int>& id_tasks_with_halo(int lv) const { return  id_tasks_with_halo_[lv]; }

    const int* id_task_array_all      (int lv) const { return  id_task_array_all_      [lv]; }
    const int* id_task_array_lbm      (int lv) const { return  id_task_array_lbm_      [lv]; }
    const int* id_task_array_L2F      (int lv) const { return  id_task_array_L2F_      [lv]; }
    const int* id_task_array_F2L      (int lv) const { return  id_task_array_F2L_      [lv]; }
    const int* id_task_array_with_halo(int lv) const { return  id_task_array_with_halo_[lv]; }

    int  num_task_array_all      (int lv) const { return  num_task_array_all_      [lv]; }
    int  num_task_array_lbm      (int lv) const { return  num_task_array_lbm_      [lv]; }
    int  num_task_array_L2F      (int lv) const { return  num_task_array_L2F_      [lv]; }
    int  num_task_array_F2L      (int lv) const { return  num_task_array_F2L_      [lv]; }
    int  num_task_array_with_halo(int lv) const { return  num_task_array_with_halo_[lv]; }

    // opt //
    // without halo //
    const std::vector<int>& id_tasks_lbm_w_bc (int lv) const { return  id_tasks_lbm_w_bc_ [lv]; }
    const std::vector<int>& id_tasks_lbm_wo_bc(int lv) const { return  id_tasks_lbm_wo_bc_[lv]; }

    const int* id_task_array_lbm_w_bc (int lv) const { return  id_task_array_lbm_w_bc_ [lv]; }
    const int* id_task_array_lbm_wo_bc(int lv) const { return  id_task_array_lbm_wo_bc_[lv]; }

    int  num_task_array_lbm_w_bc (int lv) const { return  num_task_array_lbm_w_bc_ [lv]; }
    int  num_task_array_lbm_wo_bc(int lv) const { return  num_task_array_lbm_wo_bc_[lv]; }

    // with halo //
    const std::vector<int>& id_tasks_with_halo_w_bc (int lv) const { return  id_tasks_with_halo_w_bc_ [lv]; }
    const std::vector<int>& id_tasks_with_halo_wo_bc(int lv) const { return  id_tasks_with_halo_wo_bc_[lv]; }

    const int* id_task_array_with_halo_w_bc (int lv) const { return  id_task_array_with_halo_w_bc_ [lv]; }
    const int* id_task_array_with_halo_wo_bc(int lv) const { return  id_task_array_with_halo_wo_bc_[lv]; }

    int  num_task_array_with_halo_w_bc (int lv) const { return  num_task_array_with_halo_w_bc_ [lv]; }
    int  num_task_array_with_halo_wo_bc(int lv) const { return  num_task_array_with_halo_wo_bc_[lv]; }

    // mlups //
    int  num_task_lbm_per_step() const;

public:
    void  set_taskID(const Tree&  tree);
    void  set_taskID_opt(const Tree&  tree, const MeshValue* meshValues);

    void  writeTaskID(const std::string filename) const;
    void  readTaskID(const std::string filename);

    void  copyTaskID(const TaskID& taskID);

    std::string filename(const int rank, const int step) const;

private:
    void  clear_taskID();
    void  clear_taskID_opt();

    void  init_taskID_vector(const Tree& tree);
    void  init_taskID_array ();

    void  init_taskID_vector_opt(const Tree&  tree, const MeshValue*  meshValue);
    bool leaf_has_object(
        const int    idl,
        const int*   mesh_offsets3x3x3,
        const real*  lv_obs
        );

    void  init_taskID_array_opt();


    void  write_vector(std::ofstream& fout, const std::vector<int>& vec) const;
    void  read_vector(std::ifstream& fin, std::vector<int>& vec);

};


#endif
