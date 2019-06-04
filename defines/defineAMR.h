#pragma once
#ifndef DEFINELEAF_H_
#define DEFINELEAF_H_


#include "defineCUDA.h"


namespace  DefAMR {
//    const int  LV_MAX  = 1;
//    const int  LV_MAX  = 2;
    const int  LV_MAX  = 3;

    const int  NX_LEAF = 4;
//    const int  NX_LEAF = 8;
    const int  NN_LEAF = NX_LEAF*NX_LEAF*NX_LEAF;

    inline bool is_level_max(const int lv){ return  (lv == LV_MAX-1); }
};


#endif
