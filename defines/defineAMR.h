#pragma once
#ifndef DEFINELEAF_H_
#define DEFINELEAF_H_


#include "defineCUDA.h"


namespace  DefAMR {

    // LV_MAX: by DEFINEAMR_LV_MAX, or set Oklahoma config as default
    #ifdef DEFINEAMR_LV_MAX
    constexpr int  LV_MAX  = DEFINEAMR_LV_MAX;
    #else
    constexpr int  LV_MAX  = 3;
    #endif

    const int  NX_LEAF = 4;
//    const int  NX_LEAF = 8;
    const int  NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    inline bool is_level_max(const int lv){ return  (lv == LV_MAX-1); }
};


#endif
