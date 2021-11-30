#pragma once
#ifndef VALUETIMEAVERAGE_H_
#define VALUETIMEAVERAGE_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "FuncAllocate.h"
#include "ValueNS.h"


class  ValueTimeAverage {
private:
    MemType memType_{MemType::NullPtr};

    int   nn_max_;
    real  time_min_av_ = 1.0;

    real* u_mean_ = nullptr;
    real* v_mean_ = nullptr;
    real* w_mean_ = nullptr;
    real* T_mean_ = nullptr;

    real* vel2_fluc_ = nullptr;
//    real* uu_fluc_ = nullptr;
//    real* vv_fluc_ = nullptr;
//    real* ww_fluc_ = nullptr;
//    real* TT_fluc_ = nullptr;

public:
    ValueTimeAverage () {
        memType_ = MemType::Host;
#ifdef GPU_CALCULATION__
        memType_ = MemType::Managed;
#endif
    }

    ~ValueTimeAverage () { release(); }

public:
    MemType  memType() const { return memType_; }

    real   time_min_av() const { return  time_min_av_; }

          real*  u_mean()        { return  u_mean_; }
          real*  v_mean()        { return  v_mean_; }
          real*  w_mean()        { return  w_mean_; }
          real*  T_mean()        { return  T_mean_; }

          real*  vel2_fluc()        { return  vel2_fluc_; }
//          real*  uu_fluc()        { return  uu_fluc_; }
//          real*  vv_fluc()        { return  vv_fluc_; }
//          real*  ww_fluc()        { return  ww_fluc_; }
//          real*  TT_fluc()        { return  TT_fluc_; }


    const real*  u_mean()  const { return  u_mean_; }
    const real*  v_mean()  const { return  v_mean_; }
    const real*  w_mean()  const { return  w_mean_; }
    const real*  T_mean()  const { return  T_mean_; }

    const real*  vel2_fluc()  const { return  vel2_fluc_; }
//    const real*  uu_fluc()  const { return  uu_fluc_; }
//    const real*  vv_fluc()  const { return  vv_fluc_; }
//    const real*  ww_fluc()  const { return  ww_fluc_; }
//    const real*  TT_fluc()  const { return  TT_fluc_; }

    int nn_max() const { return nn_max_; }

public:
    void  init(real time_min_av, const int  nn_max)
    {
        nn_max_ = nn_max;
        time_min_av_ = time_min_av;
    
        allocate(nn_max);
        fill(nn_max);
    }

    void  copy(const int  nn_max, const ValueTimeAverage&  other)
    {
        nn_max_ = nn_max;
        time_min_av_ = other.time_min_av();

        FuncAllocate::copy_values( u_mean(), other.u_mean(), nn_max, memType_ );
        FuncAllocate::copy_values( v_mean(), other.v_mean(), nn_max, memType_ );
        FuncAllocate::copy_values( w_mean(), other.w_mean(), nn_max, memType_ );
        FuncAllocate::copy_values( T_mean(), other.T_mean(), nn_max, memType_ );

        FuncAllocate::copy_values( vel2_fluc(), other.vel2_fluc(), nn_max, memType_ );
//        FuncAllocate::copy_values( uu_fluc(), other.uu_fluc(), nn_max, memType_ );
//        FuncAllocate::copy_values( vv_fluc(), other.vv_fluc(), nn_max, memType_ );
//        FuncAllocate::copy_values( ww_fluc(), other.ww_fluc(), nn_max, memType_ );
//        FuncAllocate::copy_values( TT_fluc(), other.TT_fluc(), nn_max, memType_ );
    }

    void  reallocate(const int  nn_max)
    {
        release();
        nn_max_  = nn_max;
        allocate(nn_max);
    }

    void copy_from_ValueNS(const int nn_max, const ValueNS&  valueNS)
    {
        nn_max_ = nn_max;
        FuncAllocate::copy_values( u_mean(), valueNS.u(), nn_max, memType_ );
        FuncAllocate::copy_values( v_mean(), valueNS.v(), nn_max, memType_ );
        FuncAllocate::copy_values( w_mean(), valueNS.w(), nn_max, memType_ );
        FuncAllocate::copy_values( T_mean(), valueNS.T(), nn_max, memType_ );
    }

private:
    void  allocate(const int nn_max)
    {
        FuncAllocate::allocate_value<real>(&u_mean_, nn_max, memType_);
        FuncAllocate::allocate_value<real>(&v_mean_, nn_max, memType_);
        FuncAllocate::allocate_value<real>(&w_mean_, nn_max, memType_);
        FuncAllocate::allocate_value<real>(&T_mean_, nn_max, memType_);

        FuncAllocate::allocate_value<real>(&vel2_fluc_, nn_max, memType_);
//        FuncAllocate::allocate_value<real>(&uu_fluc_, nn_max, memType_);
//        FuncAllocate::allocate_value<real>(&vv_fluc_, nn_max, memType_);
//        FuncAllocate::allocate_value<real>(&ww_fluc_, nn_max, memType_);
//        FuncAllocate::allocate_value<real>(&TT_fluc_, nn_max, memType_);
    }

    void  release()
    {
        FuncAllocate::release_value<real>(u_mean_, memType_);
        FuncAllocate::release_value<real>(v_mean_, memType_);
        FuncAllocate::release_value<real>(w_mean_, memType_);
        FuncAllocate::release_value<real>(T_mean_, memType_);

        FuncAllocate::release_value<real>(vel2_fluc_, memType_);
//        FuncAllocate::release_value<real>(uu_fluc_, memType_);
//        FuncAllocate::release_value<real>(vv_fluc_, memType_);
//        FuncAllocate::release_value<real>(ww_fluc_, memType_);
//        FuncAllocate::release_value<real>(TT_fluc_, memType_);
    }

    void  fill(const int nn_max)
    {
        FuncAllocate::fill_value<real>(u_mean_, 0.0, nn_max, memType_);
        FuncAllocate::fill_value<real>(v_mean_, 0.0, nn_max, memType_);
        FuncAllocate::fill_value<real>(w_mean_, 0.0, nn_max, memType_);
        FuncAllocate::fill_value<real>(T_mean_, 0.0, nn_max, memType_);

        FuncAllocate::fill_value<real>(vel2_fluc_, 0.0, nn_max, memType_);
//        FuncAllocate::fill_value<real>(uu_fluc_, 0.0, nn_max, memType_);
//        FuncAllocate::fill_value<real>(vv_fluc_, 0.0, nn_max, memType_);
//        FuncAllocate::fill_value<real>(ww_fluc_, 0.0, nn_max, memType_);
//        FuncAllocate::fill_value<real>(TT_fluc_, 0.0, nn_max, memType_);
    }

};


#endif
