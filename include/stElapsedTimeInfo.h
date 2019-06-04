#pragma once
#ifndef STELAPSEDTIMEINFO_H_
#define STELAPSEDTIMEINFO_H_


#include <string>


struct  stElapsedTimeInfo {
    std::string  func_name;

    int   rank;
    int   step;
    int   nn_leaf;
    int   lv;
    int   count;

    float elapsed_time_msec;
};


#endif
