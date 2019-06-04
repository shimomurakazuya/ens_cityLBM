#pragma once
#ifndef POSTPROCESSTAYLORGREEN_H_
#define POSTPROCESSTAYLORGREEN_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <mpi.h>
#include "defineAMR.h"
#include "Field.h"


class  PostprocessTaylorGreen {
private:
    int  rank_;

public:
    PostprocessTaylorGreen () { MPI_Comm_rank(MPI_COMM_WORLD, &rank_); }
    ~PostprocessTaylorGreen () {}

public:
    void  OutputStatistics(int t, float time, Field& field, const MeshValue*    meshValues);

private:

    std::string filename(const int t){
//        return  Foldernames::output_folder + "/" + Filenames::TaylorGreen_statics + "_step" + std::to_string(t) + ".csv";
        return  Foldernames::output_folder + "/" + Filenames::TaylorGreen_statics + ".csv";
    }
};


#endif

