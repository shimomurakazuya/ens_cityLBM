#pragma once
#ifndef TIMEEVOLUTIONLOOP_H_
#define TIMEEVOLUTIONLOOP_H_


#include <iostream>
#include <functional>
#include <vector>
#include <mpi.h>
#include "definePrecision.h"
#include "Field.h"
#include "WorkerThread.h"
#include "LBMCalculation.h"
#include "OutputFunc.h"


class  TimeEvolutionLoop {
private:
    int  rank_;

    const int  start_step_;
    const int  finish_step_;
    int  step_;

    const real start_time_;
    const real finish_time_;
    const real dt_;
    real  time_;

    // functions //
    OutputFunc      outputFunc_;
    LBMCalculation  lbmCalculation_;

public:
    TimeEvolutionLoop(
            int  start_step,
            int  finish_step,
            real start_time,
            real finish_time,
            real dt
        ) :
        start_step_(start_step),
        finish_step_(finish_step),
        start_time_(start_time),
        finish_time_(finish_time),
        dt_(dt)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

        // check dt //
//        const double _dt = (double)(finish_time - start_time) / (finish_step - start_step);
//        const double  ep = 1.0e-16;
//        if ( ! ((dt > _dt - ep) && (dt < _dt + ep)) ) {
//            std::cout << "error dt : " << __PRETTY_FUNCTION__ << std::endl;
//            std::cout << "dt = " << dt_ << ", " << _dt << std::endl;
//            exit(-1);
//        }

//        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    int  step()        const { return  step_; }
    int  start_step()  const { return  start_step_; }
    int  finish_step() const { return  finish_step_; }

    real time()        const { return  time_; }
    real start_time()  const { return  start_time_; }
    real finish_time() const { return  finish_time_; }

public:
    void time_evolution(
        Field&          field,
        WorkerThread&   workerThread
        );

private:

};


#endif
