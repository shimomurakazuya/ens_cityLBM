#pragma once
#ifndef FUNCOBJ_H_
#define FUNCOBJ_H_


#include "defineCUDA.h"
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <functional>
#include "definePrecision.h"
#include "FuncMath.h"


namespace  FuncObj {


inline
__HOST__ __DEVICE__
real  lv_fluid()
{
    return (real)(-1.0);
}


inline
__HOST__ __DEVICE__
real  lv_obj  ()
{
    return (real)1.0;
}


template<typename T>
__HOST__ __DEVICE__
bool  is_fluid(const T  lv_obj)
{
    return (lv_obj < (T)0.0) ? true : false;
}


template<typename T>
__HOST__ __DEVICE__
bool  is_included_in_fluid(const T  lv_obj, const T dx, const T factor)
{
    return (lv_obj + dx*factor < (T)0.0) ? true : false;
}


template<typename T>
__HOST__ __DEVICE__
bool  is_obj  (const T  lv_obj)
{
    return !is_fluid(lv_obj);
}


template<typename T>
__HOST__ __DEVICE__
T  overlap_objects(const T lv_obj0, const T lv_obj1)
{
    // solid or fluid //
    // s & s  ->  min //
    // s & f  ->  max //
    // f & f  ->  max //
    return ( is_obj(lv_obj0)  &&  is_obj(lv_obj1) ) ? FuncMath::min(lv_obj0, lv_obj1)
                                                    : FuncMath::max(lv_obj0, lv_obj1);
}


inline
__HOST__ __DEVICE__
real
length_from_plane(
    const real  x[],    /* x,y,z */
    const real  plane,
    const int   dim
    );


inline
__HOST__ __DEVICE__
double
length_from_box(
    const real  x[],    /* x,y,z */
    const real  box_min[],
    const real  box_max[]
    );


inline
__HOST__ __DEVICE__
double
deco_boco(
    const double  x,
    const double  y,
    const double  width,
    const double  offset_x,
    const double  offset_y
    );


inline
__HOST__ __DEVICE__
bool
is_even_position(int val)
{
    if (val >= 0) { return (std::abs(val)%2 == 0); }
    else          { return (std::abs(val)%2 == 1); }
}


};


#include "FuncObj.hpp"


#endif
