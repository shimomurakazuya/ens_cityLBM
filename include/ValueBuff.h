#pragma once
#ifndef VALUEBUFF_H_
#define VALUEBUFF_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"
#include "defineLBM.h"


class  ValueBuff {
private:
    MemType memType_{MemType::NullPtr};

    int     nQ_;
	real*	sbuff_lbm_ = nullptr;
	real*	rbuff_lbm_ = nullptr;

	real*	sbuff_lbm_host_ = nullptr;
	real*	rbuff_lbm_host_ = nullptr;

    const real  ndiv_buff_;

public:
     ValueBuff () : nQ_(LBM_velocity_model::nQ), ndiv_buff_(0.90)
     {
        memType_ = MemType::Device;
     }
    ~ValueBuff () { release(); }

public:
    MemType  memType() const { return memType_; }

          real*  sbuff_lbm()        { return  sbuff_lbm_; }
          real*  rbuff_lbm()        { return  rbuff_lbm_; }
    const real*  sbuff_lbm()  const { return  sbuff_lbm_; }
    const real*  rbuff_lbm()  const { return  rbuff_lbm_; }

          real*  sbuff_lbm_host()        { return  sbuff_lbm_host_; }
          real*  rbuff_lbm_host()        { return  rbuff_lbm_host_; }
    const real*  sbuff_lbm_host()  const { return  sbuff_lbm_host_; }
    const real*  rbuff_lbm_host()  const { return  rbuff_lbm_host_; }

public:
    void  init(const int  nn_max, const enum MemType  memType);

    void  reallocate(const MemType  memType, const int  nn_max)
    {
        release();
        memType_ = memType; // MemType::Managed //
#ifdef GPU_CALCULATION__
//        memType_ = MemType::Host;
//        memType_ = MemType::Managed;
        memType_ = MemType::Device;
#endif
        allocate(nn_max);
    }

private:
    void  allocate(const int nn_max);
    void  release();

};


#endif
