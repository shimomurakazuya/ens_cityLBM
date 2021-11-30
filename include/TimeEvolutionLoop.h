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
#include "MPICommEnsemble.h"

class  TimeEvolutionLoop {
private:
    const MPICommEnsemble comm_;

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
        const MPICommEnsemble comm,
            int  start_step,
            int  finish_step,
            real start_time,
            real finish_time,
            real dt
        ) :
        comm_(comm),
        start_step_(start_step),
        finish_step_(finish_step),
        start_time_(start_time),
        finish_time_(finish_time),
        dt_(dt),
        //
        outputFunc_(comm),
        lbmCalculation_(comm)
    {
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
