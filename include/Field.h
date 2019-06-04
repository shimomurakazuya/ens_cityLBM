#pragma once
#ifndef FIELD_H_
#define FIELD_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <fstream>
#include <mpi.h>

#include "defineAMR.h"
#include "Vector3d.h"
#include "Parser.h"
#include "Grid.h"
#include "Tree.h"
#include "Parameters.h"
#include "MeshValue.h"
#include "TaskID.h"
#include "ElapsedTimeInfo.h"
#include "TimerSimple.h"
#include "PostprocessMonitor.h"


class  Field {
private:
    const OptionParser* optionParser_;
    int  rank_;

    // grid information //
    Grid        grids_[DefAMR::LV_MAX];

    // connenctions //
    Tree        tree0_;
    Tree        tree1_;
    Tree*       tree_ptr_[2];

    // parameters //
    Parameters  parameters_;

    // values //
    MeshValue   meshValues0_[DefAMR::LV_MAX];
    MeshValue   meshValues1_[DefAMR::LV_MAX];
    MeshValue*  meshValues_ptr_[2][DefAMR::LV_MAX];

    // other //
    TimerSimple     timerSimple_;
    ElapsedTimeInfo allTimeInfo_;
    ElapsedTimeInfo funcTimeInfo_;
    ElapsedTimeInfo mpiTimeInfo_;
    TaskID          taskID_;

    // io //
    Parameters  parameters_io_;
    Grid        grids_io_[DefAMR::LV_MAX];
    Tree        tree_io_;
    TaskID      taskID_io_;
    MeshValue   meshValues_io_[DefAMR::LV_MAX];

    PostprocessMonitor postprocessmonitor_;

private: // swap counter //
    int   itp_meshval_ [DefAMR::LV_MAX];
    int   itpn_meshval_[DefAMR::LV_MAX];

    int   itp_mesh_val(int lv) const { return itp_meshval_ [lv]; }
    int   itpn_meshval(int lv) const { return itpn_meshval_[lv]; }

private:
    Field () : optionParser_(nullptr) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        init_ptrs();

        PostprocessMonitor postprocessmonitor_();
    }

public:
    Field (const OptionParser*  optionParser) :
        optionParser_(optionParser)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        init_ptrs();

        PostprocessMonitor postprocessmonitor_();
    }

    Field (const OptionParser*  optionParser, const MPI_Comm mpi_comm) :
        optionParser_(optionParser)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        init_ptrs();

        PostprocessMonitor postprocessmonitor_(mpi_comm);
    }

    ~Field () {}

public:
    const OptionParser* optionParser() const { return optionParser_; }

    // grid information //
    const Grid&  grid_i(const int i)  const { return grids_[i]; }

          Grid*  grids()        { return  grids_; }
    const Grid*  grids()  const { return  grids_; }

    // connections //
          Tree& tree()        { return  *tree_ptr_[0]; }
    const Tree& tree()  const { return  *tree_ptr_[0]; }

    // parameters //
    const Parameters& parameters() const { return parameters_; }

    // values //
          MeshValue* meshValues0()        { return  meshValues0_; }
    const MeshValue* meshValues0()  const { return  meshValues0_; }
          MeshValue* meshValues1()        { return  meshValues1_; }
    const MeshValue* meshValues1()  const { return  meshValues1_; }
          MeshValue* meshValues_io()        { return  meshValues_io_; }
    const MeshValue* meshValues_io()  const { return  meshValues_io_; }

          MeshValue& meshValue(const int lv)        { return  *meshValues_ptr_[itp_meshval_[lv]][lv]; }
    const MeshValue& meshValue(const int lv)  const { return  *meshValues_ptr_[itp_meshval_[lv]][lv]; }

          MeshValue& meshValue_new(const int lv)        { return  *meshValues_ptr_[itpn_meshval_[lv]][lv]; }
    const MeshValue& meshValue_new(const int lv)  const { return  *meshValues_ptr_[itpn_meshval_[lv]][lv]; }

    // other //
          TaskID& taskID()       { return taskID_; }
    const TaskID& taskID() const { return taskID_; }

    TimerSimple& timerSimple()  { return timerSimple_; }

    ElapsedTimeInfo& allTimeInfo()  { return allTimeInfo_; }
    ElapsedTimeInfo& funcTimeInfo() { return funcTimeInfo_; }
    ElapsedTimeInfo& mpiTimeInfo()  { return mpiTimeInfo_; }

public:
    void preset_field();
    void init_field();

    void read_field (const int step);
    void write_field(const int step) const;
    void write_monitor(const int step) const;
    void copy_io_field();

    void swap_meshValue(const int lv);

    void parameter_time(const int t);

    // copy //
    void copy_tree(Tree&  tree_new, const Tree&  tree, const bool  is_reset = false);

    void copy_meshValues(
              MeshValue*    meshValue_new,
        const MeshValue*    meshValue,
        const Grid*         grids,
        const Tree&         tree,
        const bool          is_reset = false
        );

    void update_meshValues(
        const Grid*         grids,
        const Tree&         tree
        );

    void init_taskID_opt     (TaskID& taskID, const Tree& tree, const MeshValue* meshValues);

private:
    void  init_ptrs();

    // preset //
    void preset_tree         (Tree&   tree);

    // init //
    void init_grid           (Grid*   grid);
    void init_tree           (Tree&   tree);
    void init_parameters     (Parameters& parameters);
    void init_taskID         (TaskID& taskID, const Tree& tree);

    void init_elapsedTimeInfo(const Parameters& parameter);
    void init_meshValues     (const Grid* grids, const Tree& tree, MeshValue*  meshValues);

    // read //
    int  read_latest_step    (const int step);
    void read_parameters     (const int step, Parameters&  parameters);
    void read_grid           (const int step, Grid*   grid, const int rank);
    void read_tree           (const int step, Tree&   tree, const int rank);
    void read_elapsedTimeInfo(const Parameters& parameter, const int rank);
    void read_taskID         (const int step, TaskID& taskID, const int rank);

    void read_meshValues(
        const int           step,
        const Grid*         grids,
        const Tree&         tree,
              MeshValue*    meshValues,
        const int           rank
        )
    const;

    // write //
    void write_latest_step(const int step, const int rank)  const;
    void write_parameters (const int step, const Parameters&  parameters, const int rank)  const;
    void write_grid       (const int step, const Grid*   grid, const int rank)  const;
    void write_tree       (const int step, const Tree&   tree, const int rank)  const;
    void write_taskID     (const int step, const TaskID& taskID, const int rank)  const;

    void write_meshValues(
        const int           step,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
        const MeshValue*    meshValues,
        const int           rank
        )
    const;

public:
    PostprocessMonitor& monitor(){ return postprocessmonitor_; }

};


#endif
