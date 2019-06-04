#pragma once
#ifndef DEFINE_PRECISION_H_
#define DEFINE_PRECISION_H_


typedef unsigned char BYTE;


#define SINGLE_PRECISION_CALCULATION
//#define   DOUBLE_PRECISION_CALCULATION


// single precision
#ifdef SINGLE_PRECISION_CALCULATION
    typedef float   real;
    #define MFLOAT  MPI_FLOAT
#endif


// double precision
#ifdef DOUBLE_PRECISION_CALCULATION
    typedef double  real;
    #define MFLOAT  MPI_DOUBLE
#endif


#endif
