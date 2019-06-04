#pragma once
#ifndef FUNCLBMAMR_H_
#define FUNCLBMAMR_H_


#include "defineCUDA.h"
#include <iostream>
#include "definePrecision.h"
#include "defineLBM.h"


namespace  FuncLBMAMR {


inline
__HOST__ __DEVICE__
real  feqAMR_L2X(
    const real  fs,
    const int   idv,
    const real  rhos,
    const real  us,
    const real  vs,
    const real  ws,
    const real  coefAMR
    );


inline
__HOST__ __DEVICE__
void  feqAMR_L2X(
          real  fsAMR[],
    const real  fs[],
    const real  coefAMR
    );


};

#include "FuncLBMAMR.hpp"


#endif
