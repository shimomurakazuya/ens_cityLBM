#pragma once
#ifndef TIMERSIMPLE_H_
#define TIMERSIMPLE_H_


#include <iostream>
#include <cstdlib>
#include <functional>
#include <sys/time.h>
#include "MPICommEnsemble.h"

class  TimerSimple {
private:
    int  rank_;

    struct timeval  t_cal_begin_;
    struct timeval  t_cal_end_;
    double  t_cal_total_;

    struct timeval  t_comm_begin_;
    struct timeval  t_comm_end_;
    double  t_comm_total_;

public:
    TimerSimple() = delete;

    TimerSimple(const MPICommEnsemble comm): rank_(comm.world().rank())
    {
        gettimeofday(&t_cal_begin_, NULL);
        gettimeofday(&t_cal_end_,   NULL);
        t_cal_total_ = 0.0;

        gettimeofday(&t_comm_begin_, NULL);
        gettimeofday(&t_comm_end_,   NULL);
        t_comm_total_ = 0.0;
    }

    ~TimerSimple(){}

public:
    double t_cal_total () const { return t_cal_total_; }
    double t_comm_total() const { return t_comm_total_; }

public:
    void  start_cal_timer() { gettimeofday(&t_cal_begin_, NULL); }
    void  stop_cal_timer()  { gettimeofday(&t_cal_end_,   NULL); }
    void  add_cal_timer()   { t_cal_total_ += cal_elapsed_time_sec(t_cal_begin_, t_cal_end_); }
    void  reset_cal_timer() { t_cal_total_  = 0.0; }

    void  start_comm_timer() { gettimeofday(&t_comm_begin_, NULL); }
    void  stop_comm_timer()  { gettimeofday(&t_comm_end_,   NULL); }
    void  add_comm_timer()   { t_comm_total_ += cal_elapsed_time_sec(t_comm_begin_, t_comm_end_); }
    void  reset_comm_timer() { t_comm_total_  = 0.0; }

private:
    double cal_elapsed_time_sec (
        const struct timeval&   begin,
        const struct timeval&   end
        ) const
    {
        // sec //
        return (   (end.tv_sec  - begin.tv_sec) * 1000
                 + (end.tv_usec - begin.tv_usec) / 1000.0 ) / 1000.0;
    }

};


#endif
