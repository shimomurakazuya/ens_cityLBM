#pragma once
#ifndef VALUELBM_H_
#define VALUELBM_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"
#include "defineLBM.h"


class  ValueLBM {
private:
    MemType memType_{MemType::NullPtr};

    int     nQ_;
	real*	f_lbm_ = nullptr;

public:
    ValueLBM () : nQ_(LBM_velocity_model::nQ) {}

    ~ValueLBM () {
        release();
    }

public:
    MemType  memType() const { return memType_; }

          real*  f_lbm()        { return  f_lbm_; }
    const real*  f_lbm()  const { return  f_lbm_; }

public:
    void  init(const int  nn_max, const enum MemType  memType);
    void  copy(const int  nn_max, const ValueLBM&  other);

    void  reallocate(const MemType  memType, const int  nn_max)
    {
        release();
        memType_ = memType;
        allocate(nn_max);
    }

private:
    void  allocate(const int nn_max);
    void  release();

    void  fill(const int nn_max);

};


#endif
