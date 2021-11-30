#include "FuncScheduler.h"


namespace  FuncScheduler {


void
launcher(
    const int  lv_max,
    std::function<void(int)>&  func_cal,
    std::function<void(int)>&  func_L2F,
    std::function<void(int)>&  func_L2FA,
    std::function<void(int)>&  func_F2L
    )
{
    worker(0, lv_max, func_cal, func_L2F, func_L2FA, func_F2L);
}


void
worker(
    const int  lv,
    const int  lv_max,
    std::function<void(int)>&  func_cal,
    std::function<void(int)>&  func_L2F,
    std::function<void(int)>&  func_L2FA,
    std::function<void(int)>&  func_F2L
    )
{
    // node @ lv //
    // i :: smaller is coarser //
    if (lv == lv_max-1) { // top level //
        func_cal(lv);
    }
    else {
        // fine level //
        worker(lv+1, lv_max, func_cal, func_L2F, func_L2FA, func_F2L);

        // this level //
        func_cal(lv);
        func_L2FA(lv);

        // fine level //
        worker(lv+1, lv_max, func_cal, func_L2F, func_L2FA, func_F2L);
        func_F2L(lv);  func_L2F(lv);
    }
}


};
