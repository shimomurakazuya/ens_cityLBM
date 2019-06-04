#pragma once
#ifndef VALUEOBJLS_H_
#define VALUEOBJLS_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "defineObjBC.h"
#include "Grid.h"


class  ValueObjLS {
private:
    MemType memType_{MemType::NullPtr};

    // solid wall //
    real*   lv_obj_  = nullptr;

    real*   rho_obj_      = nullptr;
    real*   u_obj_        = nullptr;
    real*   v_obj_        = nullptr;
    real*   w_obj_        = nullptr;
    real*   T_obj_        = nullptr;
    real*   scalar_obj_   = nullptr;

    real*   rhon_obj_      = nullptr;
    real*   un_obj_        = nullptr;
    real*   vn_obj_        = nullptr;
    real*   wn_obj_        = nullptr;
    real*   Tn_obj_        = nullptr;
    real*   scalarn_obj_   = nullptr;

    real    time_obj_ { 0.0 };
    real    timen_obj_{ 1.0 };

    // inflow & outflow bc //
//    int*  bcTypes_scaler_;
//    int*  bcTypes_vector_;
    int*  bcTypes_f_;
    int*  bcTypes_T_;

//    real*     bc_rho_ = nullptr;
//    real*     bc_u_   = nullptr;
//    real*     bc_v_   = nullptr;
//    real*     bc_w_   = nullptr;
//    real*     bc_weight_ = nullptr; // weight factor of bc values : 1 is the same as the bc value //

    real*     dirichlet_weight_ = nullptr; // weight factor of bc values : 1 is the same as the bc value //
    real*     viscosity_weight_ = nullptr; // weight factor of bc values : 1 is the same as the bc value //

    // source & diffuse //
    real*     sc_scalar_ = nullptr;

public:
    ValueObjLS () {
    }

    ~ValueObjLS () {
        release();
    }

public:
    MemType  memType() const { return memType_; }

    // solid wall //
          real*  lv_obj()        { return  lv_obj_; }
    const real*  lv_obj()  const { return  lv_obj_; }

          real*  rho_obj()        { return  rho_obj_; }
    const real*  rho_obj()  const { return  rho_obj_; }

          real*  scalar_obj()        { return  scalar_obj_; }
    const real*  scalar_obj()  const { return  scalar_obj_; }

          real*  u_obj()        { return  u_obj_; }
          real*  v_obj()        { return  v_obj_; }
          real*  w_obj()        { return  w_obj_; }
    const real*  u_obj()  const { return  u_obj_; }
    const real*  v_obj()  const { return  v_obj_; }
    const real*  w_obj()  const { return  w_obj_; }
          real*  T_obj()        { return  T_obj_; }
    const real*  T_obj()  const { return  T_obj_; }

    // new //
          real*  rhon_obj()        { return  rhon_obj_; }
    const real*  rhon_obj()  const { return  rhon_obj_; }

          real*  scalarn_obj()        { return  scalarn_obj_; }
    const real*  scalarn_obj()  const { return  scalarn_obj_; }

          real*  un_obj()        { return  un_obj_; }
          real*  vn_obj()        { return  vn_obj_; }
          real*  wn_obj()        { return  wn_obj_; }
    const real*  un_obj()  const { return  un_obj_; }
    const real*  vn_obj()  const { return  vn_obj_; }
    const real*  wn_obj()  const { return  wn_obj_; }
          real*  Tn_obj()        { return  Tn_obj_; }
    const real*  Tn_obj()  const { return  Tn_obj_; }

    // inflow & outflow bc //
//          int*  bcTypes_scaler()        { return  bcTypes_scaler_; }
//          int*  bcTypes_vector()        { return  bcTypes_vector_; }
//    const int*  bcTypes_scaler()  const { return  bcTypes_scaler_; }
//    const int*  bcTypes_vector()  const { return  bcTypes_vector_; }

          int*  bcTypes_f()        { return  bcTypes_f_; }
    const int*  bcTypes_f()  const { return  bcTypes_f_; }

          int*  bcTypes_T()        { return  bcTypes_T_; }
    const int*  bcTypes_T()  const { return  bcTypes_T_; }

//          real*  bc_rho()        { return  bc_rho_; }
//    const real*  bc_rho()  const { return  bc_rho_; }
//
//          real*  bc_u()        { return  bc_u_; }
//          real*  bc_v()        { return  bc_v_; }
//          real*  bc_w()        { return  bc_w_; }
//    const real*  bc_u()  const { return  bc_u_; }
//    const real*  bc_v()  const { return  bc_v_; }
//    const real*  bc_w()  const { return  bc_w_; }
//
//          real*  bc_weight()        { return  bc_weight_; }
//    const real*  bc_weight()  const { return  bc_weight_; }

          real*  dirichlet_weight()        { return  dirichlet_weight_; }
    const real*  dirichlet_weight()  const { return  dirichlet_weight_; }

          real*  viscosity_weight()        { return  viscosity_weight_; }
    const real*  viscosity_weight()  const { return  viscosity_weight_; }

    // source //
          real*  sc_scalar()        { return  sc_scalar_; }
    const real*  sc_scalar()  const { return  sc_scalar_; }


public:
    void  init(const int  nn_max, const enum MemType  memType);
    void  copy(const int  nn_max, const ValueObjLS&  other);

    void  reallocate(const MemType  memType, const int  nn_max)
    {
        release();
        memType_ = memType;
        allocate(nn_max);
    }


    // time //
    void set_obj_time (const real time) { time_obj_  = time; }
    void set_objn_time(const real time) { timen_obj_ = time; }

    void update_obj_time(const real time) {
        time_obj_  = timen_obj_;
        timen_obj_ = time;
    }

    void swap_objs();

    bool  check_obj_time(const real t, const real time0, const real time1)
    {
        return (t >= time_obj_ && t < timen_obj_);
    }


private:
    void  allocate(const int nn_max);
    void  release ();

    void  fill(const int nn_max);

};


#endif
