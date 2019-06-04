#pragma once
#ifndef FUNCLBMSGS_H_
#define FUNCLBMSGS_H_


#include "defineCUDA.h"
#include <cmath>
#include "definePrecision.h"


namespace  FuncLBMSGS {


inline
__HOST__ __DEVICE__
real
sgs_viscosity_D3Q27(
    const real  Fcs,
    const real  Csgs,
    const real  rho,
    const real  vis,
    const real  fs_deq[]
    );


inline
__HOST__ __DEVICE__
real
SS_D3Q27(
    const real  rho,
    const real  vis,
    const real  fs_deq[]
    );


};


#include "FuncLBMSGS.hpp"


#endif
