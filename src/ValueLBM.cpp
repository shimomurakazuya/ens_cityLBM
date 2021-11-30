#include "ValueLBM.h"
#include "FuncMath.h"
#include "FuncAllocate.h"


// public //
void  ValueLBM::
init(const int  nn_max, const enum MemType  memType)
{
    memType_ = memType;

    nQ_ = LBM_velocity_model::nQ;

    allocate(nn_max);
    fill(nn_max);
}


void  ValueLBM::
copy(const int  nn_max, const ValueLBM&  other)
{
    const int  nn_max_lbm = nn_max*nQ_;

    FuncAllocate::copy_values( f_lbm(), other.f_lbm(), nn_max_lbm, memType_ );
}


// private //
void  ValueLBM::
allocate(const int nn_max)
{
    const int  nn_max_lbm = nn_max*nQ_;

    FuncAllocate::allocate_value<real>(&f_lbm_, nn_max_lbm, memType_);
}


void  ValueLBM::
release()
{
    FuncAllocate::release_value<real>(f_lbm_, memType_);
}


void  ValueLBM::
fill(const int nn_max)
{
    const int  nn_max_lbm = nn_max*nQ_;
    const real lbm_const  = 1.0e-4;

    FuncAllocate::fill_value<real>(f_lbm_, lbm_const, nn_max_lbm, memType_);
}
