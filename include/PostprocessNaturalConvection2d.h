#pragma once
#ifndef POSTPROCESSNATURALCONVECTION2D_H_
#define POSTPROCESSNATURALCONVECTION2D_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <mpi.h>
#include "defineAMR.h"
#include "Field.h"


class  PostprocessNaturalConvection2d {
private:
    real  dx_[DefAMR::LV_MAX];
    double uT_;
    double Tx_;
    double Nu_;

    int    rank_;

public:
    PostprocessNaturalConvection2d () { 
        //MPI_Comm_rank_(MPI_COMM_WORLD, &rank_); 
        rank_ = -1;
    }
    ~PostprocessNaturalConvection2d () {}

public:
    void  OutputStatistics(int t, Field& field, const MeshValue*    meshValues);

private:
    void  Nusselt_number_lv(
              double&    uT,
              double&    Tx,
              double&    vol,
        const int        lv,
        const Field&     field,
        const MeshValue* meshValues,
        const TaskID::vector_type& id_tasks
        );

    void  output_values(int step);

    std::string filename(const int t){
        return  Foldernames::output_folder + "/" + Filenames::NaturalConvection2d_statics + "_step" + std::to_string(t) + ".csv";
    }
};


#endif
