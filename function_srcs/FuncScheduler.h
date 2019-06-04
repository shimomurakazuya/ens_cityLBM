#pragma once
#ifndef FUNCSCHEDULER_H_
#define FUNCSCHEDULER_H_


#include <iostream>
#include <string>
#include <vector>
#include <functional>


namespace  FuncScheduler {


void
launcher(
    const int  lv_max,
    std::function<void(int)>&  func_cal,
    std::function<void(int)>&  func_L2F,
    std::function<void(int)>&  func_L2FA,
    std::function<void(int)>&  func_F2L
    );


void
worker(
    const int  lv,
    const int  lv_max,
    std::function<void(int)>&  func_cal,
    std::function<void(int)>&  func_L2F,
    std::function<void(int)>&  func_L2FA,
    std::function<void(int)>&  func_F2L
    );


};


#endif
