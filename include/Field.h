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
#include "ValueBuff.h"
#include "ValueStat.h"
#include "ValueTimeAverage.h"
#include "ValuePBVR.h"
#include "TaskID.h"
#include "ElapsedTimeInfo.h"
#include "TimerSimple.h"
#include "FuncMapData.h"
#include "MPICommEnsemble.h"
#include "mpi_wrapper.hpp"

#include "FastPostprocess.h"

class  Field {
private:
    const OptionParser& optionParser_;
    const MPICommEnsemble comm_;

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

    ValueBuff   valueBuff_;

    ValueStat valueStat_[DefAMR::LV_MAX];

    ValuePBVR valuePBVR_;

    ValueTimeAverage valueTimeAverage1min_ [DefAMR::LV_MAX];
//    ValueTimeAverage valueTimeAverage15min_[DefAMR::LV_MAX];

    // other //
    TimerSimple     timerSimple_;
    ElapsedTimeInfo allTimeInfo_;
    ElapsedTimeInfo funcTimeInfo_;
    ElapsedTimeInfo mpiTimeInfo_;
    TaskID          taskID_;

    // io //
    #ifndef NO_FIELD_WRITEDAT
    Parameters  parameters_io_;
    Grid        grids_io_[DefAMR::LV_MAX];
    Tree        tree_io_;
    TaskID      taskID_io_;
    MeshValue   meshValues_io_[DefAMR::LV_MAX];
    ValueStat   valueStat_io_[DefAMR::LV_MAX];
    ValueTimeAverage valueTimeAverage1min_io_ [DefAMR::LV_MAX];
    #endif

    FastPostprocess postprocessmonitor_;

private: // swap counter //
    int   itp_meshval_ [DefAMR::LV_MAX];
    int   itpn_meshval_[DefAMR::LV_MAX];

    int   itp_mesh_val(int lv) const { return itp_meshval_ [lv]; }
    int   itpn_meshval(int lv) const { return itpn_meshval_[lv]; }

public:

    Field() = delete;

    Field (const OptionParser&  optionParser, const MPICommEnsemble comm) :
        optionParser_(optionParser),
        comm_(comm),
        // 
        tree0_(comm), tree1_(comm), 
        parameters_(comm), 
        allTimeInfo_(comm), funcTimeInfo_(comm), mpiTimeInfo_(comm),
        postprocessmonitor_(comm),
        timerSimple_(comm)
        #ifndef NO_FIELD_WRITEDAT
        , tree_io_(comm), parameters_io_(comm)
        #endif
    {
        init_ptrs();
        init_parameters(parameters_);
    }

    ~Field () {}

public:
    const OptionParser& optionParser() const { return optionParser_; }

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
    #ifndef NO_FIELD_WRITEDAT
          MeshValue* meshValues_io()        { return  meshValues_io_; }
    const MeshValue* meshValues_io()  const { return  meshValues_io_; }
    #endif

          MeshValue& meshValue(const int lv)        { return  *meshValues_ptr_[itp_meshval_[lv]][lv]; }
    const MeshValue& meshValue(const int lv)  const { return  *meshValues_ptr_[itp_meshval_[lv]][lv]; }

          ValueBuff& valueBuff()        { return valueBuff_; }
    const ValueBuff& valueBuff()  const { return valueBuff_; }

          MeshValue& meshValue_new(const int lv)        { return  *meshValues_ptr_[itpn_meshval_[lv]][lv]; }
    const MeshValue& meshValue_new(const int lv)  const { return  *meshValues_ptr_[itpn_meshval_[lv]][lv]; }

          ValueStat& valueStat   (const int lv)        { return  valueStat_   [lv]; }
    const ValueStat& valueStat   (const int lv)  const { return  valueStat_   [lv]; }
    #ifndef NO_FIELD_WRITEDAT
          ValueStat& valueStat_io(const int lv)        { return  valueStat_io_[lv]; }
    const ValueStat& valueStat_io(const int lv)  const { return  valueStat_io_[lv]; }
    #endif

          ValuePBVR& valuePBVR()        { return valuePBVR_; }
    const ValuePBVR& valuePBVR()  const { return valuePBVR_; }

          ValueTimeAverage* valueTimeAverage1min   ()        { return  valueTimeAverage1min_; }
    const ValueTimeAverage* valueTimeAverage1min   ()  const { return  valueTimeAverage1min_; }
          ValueTimeAverage& valueTimeAverage1min   (const int lv)        { return  valueTimeAverage1min_[lv]; }
    const ValueTimeAverage& valueTimeAverage1min   (const int lv)  const { return  valueTimeAverage1min_[lv]; }

    // other //
          TaskID& taskID()       { return taskID_; }
    const TaskID& taskID() const { return taskID_; }

    TimerSimple& timerSimple()  { return timerSimple_; }

    ElapsedTimeInfo& allTimeInfo()  { return allTimeInfo_; }
    ElapsedTimeInfo& funcTimeInfo() { return funcTimeInfo_; }
    ElapsedTimeInfo& mpiTimeInfo()  { return mpiTimeInfo_; }

public:
    void read_field_full (const int step);
    void read_field_params (const int step);
    void read_field_meshval (const int step);
    void write_field(const int step);
    void write_monitor(const int step);
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

    void copy_valueStat_io();
    void copy_valueTimeAverage_io();

    void update_meshValues(
        const Grid*         grids,
        const Tree&         tree
        ) = delete;

    void update_meshValues_io(
        const Grid*         grids,
        const Tree&         tree,
        const bool          is_reset
        );

    void init_taskID_opt     (TaskID& taskID, const Tree& tree, const MeshValue* meshValues);

private:
    void  init_ptrs();

    // preset //
    void preset_field() = delete;
    void init_field() = delete;
    void preset_tree(Tree&) = delete;
    public: void init_field_wo_map();

    // init with Oklahoma or other map //
    public: void init_field_with_map(const MapData& map);
private:
    void init_tree_with_map  (Tree&   tree, const MapData& map);

    // init (maybe enable to reuse) //
    void init_others(); // after init_tree_*()
    void init_grid           (Grid*   grid); // init grid for level 0 tree_roots
    void init_tree           (Tree&   tree) = delete;
    void init_parameters     (Parameters& parameters);
    void init_taskID         (TaskID& taskID, const Tree& tree);

    void init_elapsedTimeInfo(const Parameters& parameter);
    void init_meshValues     (const Grid* grids, const Tree& tree, MeshValue*  meshValues);
    void init_ValueStat      (const Tree& tree);
    void init_valueBuff     (const Grid* grids, const Tree& tree, ValueBuff& valueBuff);
    void init_valuePBVR     (ValuePBVR& valuePBVR);
    void init_ValueTimeAverage(const Tree& tree, const MeshValue* meshValues);

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
              ValueTimeAverage* valueTimeAverage,
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
    auto& monitor(){ return postprocessmonitor_; }

    void write_ValueStatCsv(int step);
};


#endif
