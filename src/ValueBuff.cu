#include "ValueBuff.h"
#include "FuncMath.h"
#include "FuncAllocate.h"


// public //
void  ValueBuff::
init(const int  nn_max, const enum MemType  memType)
{
    memType_ = memType;
#ifdef GPU_CALCULATION__
    memType_ = MemType::Device;
#endif

    nQ_ = LBM_velocity_model::nQ;

    allocate(nn_max);
}


// private //
void  ValueBuff::
allocate(const int nn_max)
{
    const int  nn_max_lbm_buff = nn_max*nQ_ * ndiv_buff_;

    FuncAllocate::allocate_value<real>(&sbuff_lbm_, nn_max_lbm_buff, memType_);
    FuncAllocate::allocate_value<real>(&rbuff_lbm_, nn_max_lbm_buff, memType_);

    FuncAllocate::allocate_value<real>(&sbuff_lbm_host_, nn_max_lbm_buff, MemType::Host);
    FuncAllocate::allocate_value<real>(&rbuff_lbm_host_, nn_max_lbm_buff, MemType::Host);
}


void  ValueBuff::
release()
{
    FuncAllocate::release_value<real>(sbuff_lbm_, memType_);
    FuncAllocate::release_value<real>(rbuff_lbm_, memType_);

    FuncAllocate::release_value<real>(sbuff_lbm_host_, MemType::Host);
    FuncAllocate::release_value<real>(rbuff_lbm_host_, MemType::Host);
}
