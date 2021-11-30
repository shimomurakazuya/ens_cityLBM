#pragma once
#ifndef POSTPROCESSNATURALCONVECTION3D_H_
#define POSTPROCESSNATURALCONVECTION3D_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <mpi.h>
#include "defineAMR.h"
#include "Field.h"


class  PostprocessNaturalConvection3d {
private:
    struct value_cf {
        bool   is_exist;
        bool   is_existg;
        int    num;
        int    lv;
        real   z;
        real   u,  v,  w;
        real   uu, vv, ww;
        real   T;
    };

    std::vector<value_cf>  values_  [DefAMR::LV_MAX];
    std::vector<value_cf>  values_g_[DefAMR::LV_MAX];

    int   numz_[DefAMR::LV_MAX];
    real  dx_[DefAMR::LV_MAX];

    int    rank_;

    const real  x_min_nc3d_;
    const real  y_min_nc3d_;
    const real  z_min_nc3d_;

    const real  x_max_nc3d_;
    const real  y_max_nc3d_;
    const real  z_max_nc3d_;

public:
    PostprocessNaturalConvection3d () :
        x_min_nc3d_(-0.0025),
        x_max_nc3d_( 0.0025),
        y_min_nc3d_(-0.0025),
        y_max_nc3d_( 0.0025),
        z_min_nc3d_( 0.000),
        z_max_nc3d_( 0.800)
    {
        //MPI_Comm_rank_(MPI_COMM_WORLD, &rank_);
        rank_ = -1;
    }
    ~PostprocessNaturalConvection3d () {}

public:
    void  OutputStatistics(int t, Field& field, const MeshValue*    meshValues);

private:
    void  exist_flag(
        const int        lv,
        const Field&     field,
        const MeshValue* meshValues,
        const TaskID::vector_type& id_tasks
        );

    void  mean_values(
        const int        lv,
        const Field&     field,
        const MeshValue* meshValues,
        const TaskID::vector_type& id_tasks
        );

    void  rms_values(
        const int        lv,
        const Field&     field,
        const MeshValue* meshValues,
        const TaskID::vector_type& id_tasks
        );

    void  init_value_cf();
    void  add_exist_flag (const int lv, const int idz, const real z_obj);
    void  add_mean_values(const int lv, const int idz, const real z_obj, const real u, const real v, const real w, const real T);
    void  add_rms_values (const int lv, const int idz, const real z_obj, const real u, const real v, const real w);

    void  reduce_exist_flag(const int lv);
    void  reduce_z         (const int lv);
    void  average_num        (const int lv);
    void  average_mean_values(const int lv);
    void  average_rms_values (const int lv);

    void  output_values(int step);

    void  get_idz_zobj_nc3d(int& idz, real& z_obj, const real z, const real dx, const real lv_obj);

    std::string filename(const int t){
        return  Foldernames::output_folder + "/" + Filenames::NaturalConvection3d_statics + "_step" + std::to_string(t) + ".csv";
    }

    bool is_center(const real x, const real y);
};


#endif
