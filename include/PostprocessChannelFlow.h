#pragma once
#ifndef POSTPROCESSCHANNELFLOW_H_
#define POSTPROCESSCHANNELFLOW_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <mpi.h>
#include "defineAMR.h"
#include "Field.h"


class  PostprocessChannelFlow {
private:
    struct value_cf {
        bool    is_exist;
        int     num;
        int     lv;
        real   z;
        real   u,  v,  w;
        real   uu, vv, ww;
        real   uv, vw, wu;


        real   T;
        real   TT;
        real   uT, vT, wT;
    };

    std::vector<value_cf>  values_  [DefAMR::LV_MAX];
    std::vector<value_cf>  values_g_[DefAMR::LV_MAX];

    int    numz_[DefAMR::LV_MAX];
    real  dx_[DefAMR::LV_MAX];

    int    rank_;

    const real  z_min_channel_;
    const real  z_max_channel_;

public:
    PostprocessChannelFlow () :
        z_min_channel_(-1.0),
        z_max_channel_( 1.0)
    {
        //MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        rank_ = -1;
    }
    ~PostprocessChannelFlow () {}

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
    void  add_rms_values (const int lv, const int idz, const real z_obj, const real u, const real v, const real w, const real T);
    void  average_num        (const int lv);
    void  average_mean_values(const int lv);
    void  average_rms_values (const int lv);

    void  output_values(int step);

    void  z_obj_channel_flow(int& idz, real& z_obj, const real z, const real dx, const real lv_obj);

    std::string filename(const int t){
        return  Foldernames::output_folder + "/" + Filenames::ChannelFlow_statics + "_step" + std::to_string(t) + ".csv";
    }
};


#endif
