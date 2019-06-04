#pragma once
#ifndef OUTPUTFUNC_H_
#define OUTPUTFUNC_H_


#include <iostream>
#include <cstdlib>
#include <functional>
#include <sys/time.h>
#include "Field.h"
#include "WorkerThread.h"


class  OutputFunc {
private:
    int  rank_;

    struct timeval  t_begin_;
    struct timeval  t_mid_;
    struct timeval  t_end_;


public:
    OutputFunc(){
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

        gettimeofday(&t_begin_, NULL);
        gettimeofday(&t_mid_,   NULL);
        gettimeofday(&t_end_,   NULL);
    }
    ~OutputFunc(){}

public:
    void OutputFluidData(int t, Field& field, WorkerThread& iothread);

private:
    void cout_step(int t, const Parameters& parameters) const;
    void cout_monitor(int t, const Field& field) const;
    void cout_timer(int t, Field& field);

    void output_data(int t, Field& field, WorkerThread& iothread) const;

    void output_monitor_data(int t, Field& field) const;
    void output_channel_flow_data(int t, Field& field) const;
    void output_taylor_green_data(int t, Field& field) const;
    void output_natural_convection2d_data(int t, Field& field) const;
    void output_natural_convection3d_data(int t, Field& field) const;

    double  cal_elapsed_time_msec (
        const struct timeval&   begin,
        const struct timeval&   end
        )
        const;

};


#endif
