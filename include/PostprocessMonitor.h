#pragma once
#ifndef POSTPROCESSMONITOR_H_
#define POSTPROCESSMONITOR_H_

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <mpi.h>

#include "defineAMR.h"
#include "Vector3d.h"
#include "Parser.h"
#include "Grid.h"
#include "Tree.h"
#include "Parameters.h"
#include "MeshValue.h"

#include "MonitorData.h"


class PostprocessMonitor {
private:
    MPI_Comm mpi_communicator_{MPI_COMM_WORLD};
    std::string in_filename_;
    std::string out_filename_;
    std::vector<MonitorData> monitordata_;

private:
    void Sendmonitordata(const int);

public:
    PostprocessMonitor(){
//        MPI_Comm_dup(MPI_COMM_WORLD, &mpi_communicator_);
        mpi_communicator_ = MPI_COMM_WORLD;

        in_filename_ = "input-";
        out_filename_= "output-timestep";
        MPI_Barrier(mpi_communicator_);
    }
    PostprocessMonitor(MPI_Comm mpi_comm){
        mpi_communicator_ = mpi_comm;
        in_filename_ = "input-";
        out_filename_= "output-timestep";
    }

    MonitorData& monitor(int index){
        return monitordata_.at(index);
    }

    void
    monitorset(
        int index,
        float u, float v, float w,
        float levelset_obj,
        float scalar,
        float T,
        int cflg
        )
    {
        monitordata_.at(index).setvalue(
            u,v,w,
            levelset_obj,
            scalar,
            T,
            cflg
            );
    }

    int monitorcount() const { return monitordata_.size();}

    void setupdata(int,const Tree&, const MeshValue*);
    void readdata();

    void OutputMonitorData(int step, const  Tree& tree, const Parameters& parameters, const  MeshValue* meshValues, int rank) const;

private:
	std::vector<std::string> split(std::string& input, char delimiter);
};


#endif
