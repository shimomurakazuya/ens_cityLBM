#pragma once
#ifndef VALUENS_H_
#define VALUENS_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"


class  ValueNS {
private:
    MemType memType_{MemType::NullPtr};

    // veloctiy //
    real*   u_ = nullptr;
    real*   v_ = nullptr;
    real*   w_ = nullptr;

    // density //
    real*   rho_ = nullptr;

    // scalar values //
    real*   scalar_ = nullptr;

    // temperature //
    real*   T_ = nullptr;

    // sgs //
    real*   sgs_vis_ = nullptr;

public:
    ValueNS () {}

    ~ValueNS () {
        release();
    }

public:
    MemType  memType() const { return memType_; }

          real*  u()        { return  u_; }
          real*  v()        { return  v_; }
          real*  w()        { return  w_; }
    const real*  u()  const { return  u_; }
    const real*  v()  const { return  v_; }
    const real*  w()  const { return  w_; }

          real*  rho()        { return  rho_; }
    const real*  rho()  const { return  rho_; }

          real*  scalar()        { return  scalar_; }
    const real*  scalar()  const { return  scalar_; }

          real*  T()        { return  T_; }
    const real*  T()  const { return  T_; }

          real*  sgs_vis()        { return  sgs_vis_; }
    const real*  sgs_vis()  const { return  sgs_vis_; }

    real  u(const int i)  const { return  u()[i]; }
    real  v(const int i)  const { return  v()[i]; }
    real  w(const int i)  const { return  w()[i]; }

    real  rho(const int i)  const { return  rho()[i]; }
    real  T  (const int i)  const { return  T  ()[i]; }

public:
    void  init(const int  nn_max, const enum MemType  memType);
    void  copy(const int  nn_max, const ValueNS&  other);

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
