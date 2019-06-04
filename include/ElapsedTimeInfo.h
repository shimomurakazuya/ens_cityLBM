#pragma once
#ifndef ELAPSEDTIMEINFO_H_
#define ELAPSEDTIMEINFO_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <functional>
#include <fstream>
#include <map>
#include <mpi.h>

#include "defineAMR.h"
#include "defineFilenames.h"
#include "stElapsedTimeInfo.h"
#include "Parameters.h"
#include "Timer.h"


class  ElapsedTimeInfo {
private:
    const bool  flag_measure_;

    int   output_period_ = 1;
    int   rank_{0};
    bool  is_ios_base_app_{false};

    std::string  fname_;

    std::vector<Timer>  timer_;
    std::map<std::string, Timer>  map_timer_;

    std::vector<stElapsedTimeInfo>  st_elapsedTimeInfo_;

public:
    ElapsedTimeInfo () :
        flag_measure_(true),
        output_period_(1)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    ~ElapsedTimeInfo () {}

public:
    void InitElapsedTimeInfo(const std::string& fname, const Parameters& parameters);
    void ReadElapsedTimeInfo(const std::string& fname, const Parameters& parameters);
    void OutputElapsedTimeInfo(int t);

    void StartTimer(const std::string& func_name);
    void StopTimer(const std::string& func_name);

    void SubmitElapsedTimeInfo(const std::string& func_name, int step);
    void SubmitElapsedTimeInfo(const std::string& func_name, int lv, int count, int step);

private:
    void setOutputPeriod(const Parameters& parameters);

private:
    void check_namelist(const std::string& func_name);
    std::string filename(const int rank) const;

};


#endif
