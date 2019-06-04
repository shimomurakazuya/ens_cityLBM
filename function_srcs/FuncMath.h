#pragma once
#ifndef FUNCMATH_H_
#define FUNCMATH_H_


#include <iostream>
#include <cstdlib>
#include "defineCal.h"
#include "FuncLoop.h"


namespace  FuncMath {


template <typename T0>
void copy_array(T0* val0, const T0* val1, const int n)
{
    std::copy(val1, val1+n, val0);
}


template <typename T0>
void copy_i(T0* val0, const T0* val1, const int n)
{
    for (int i=0; i<n; i++) {  val0[i] = val1[i];  }
}


template <typename T0>
void copy_i_omp(T0* val0, const T0* val1, const int n)
{
#pragma omp parallel for
    for (int i=0; i<n; i++) {  val0[i] = val1[i];  }
}


template <typename T0, typename T1>
void copy_i(T0* val0, const T1* val1, const int n)
{
    for (int i=0; i<n; i++) {  val0[i] = val1[i];  }
}


template <typename T0, typename T1>
void copy_i_omp(T0* val0, const T1* val1, const int n)
{
#pragma omp parallel for
    for (int i=0; i<n; i++) {  val0[i] = val1[i];  }
}


template <typename T>
void
swap (T**   f, T**  fn)
{
    T*  tmp;

    tmp = *f;
    *f  = *fn;
    *fn = tmp;
}


template <typename T>
__HOST__ __DEVICE__
T min(const T x, const T y)
{
    return x < y ?  x : y;
}


template <typename T>
__HOST__ __DEVICE__
T max(const T x, const T y)
{
    return x > y ?  x : y;
}


template <typename T>
__HOST__ __DEVICE__
T max3(const T x, const T y, const T z)
{
    const T tmp = (x > y) ?  x : y;
    return z > tmp ?  z : tmp;
}


template <typename T>
__HOST__ __DEVICE__
T abs(const T x)
{
    return x < 0 ? -x : x;
}


template <int x>
__HOST__ __DEVICE__
int abs()
{
    return x < 0 ? -x : x;
}


template <typename T>
__HOST__ __DEVICE__
T abs(const T x, const T y, const T z)
{
    return  abs(x) + abs(y) + abs(z);
}


template <int x, int y, int z>
__HOST__ __DEVICE__
int abs()
{
    return  abs<x>() + abs<y>() + abs<z>();
}


//template <typename First, typename... Rest>
//constexpr First abs(const First& first, const Rest&... rest)
//{
//    return  abs(first) + abs(rest...);
//}


template <typename T>
constexpr T norm(const T x)
{
    return  x*x;
}


template <typename First, typename... Rest>
constexpr First norm(const First& first, const Rest&... rest)
{
    return  norm(first) + norm(rest...);
}


};


#endif
