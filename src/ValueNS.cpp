#include "ValueNS.h"
#include "FuncMath.h"
#include "FuncAllocate.h"


// public //
void  ValueNS::
init(const int  nn_max, const int n_scalars, const enum MemType  memType)
{
    memType_ = memType;
    nn_max_ = nn_max;

    allocate(nn_max, n_scalars);
    fill(nn_max, n_scalars);
}


void  ValueNS::
copy(const int  nn_max, const int n_scalars, const ValueNS&  other)
{
    nn_max_ = nn_max;
    FuncAllocate::copy_values( u(), other.u(), nn_max, memType_ );
    FuncAllocate::copy_values( v(), other.v(), nn_max, memType_ );
    FuncAllocate::copy_values( w(), other.w(), nn_max, memType_ );

    FuncAllocate::copy_values( rho(), other.rho(), nn_max, memType_ );

    FuncAllocate::copy_values( scalar (), other.scalar (), nn_max * n_scalars, memType_ );

    FuncAllocate::copy_values( T(),   other.T(),   nn_max, memType_ );

    FuncAllocate::copy_values( sgs_vis(),   other.sgs_vis(),   nn_max, memType_ );
}


// private //
void  ValueNS::
allocate(const int nn_max, const int n_scalars)
{
    FuncAllocate::allocate_value<real>(&u_, nn_max, memType_);
    FuncAllocate::allocate_value<real>(&v_, nn_max, memType_);
    FuncAllocate::allocate_value<real>(&w_, nn_max, memType_);

    FuncAllocate::allocate_value<real>(&rho_, nn_max, memType_);

    FuncAllocate::allocate_value<real>(&scalar_ , nn_max * n_scalars, memType_);

    FuncAllocate::allocate_value<real>(&T_,   nn_max, memType_);

    FuncAllocate::allocate_value<real>(&sgs_vis_,   nn_max, memType_);
}


void  ValueNS::
release()
{
    FuncAllocate::release_value<real>(u_, memType_);
    FuncAllocate::release_value<real>(v_, memType_);
    FuncAllocate::release_value<real>(w_, memType_);

    FuncAllocate::release_value<real>(rho_, memType_);

    FuncAllocate::release_value<real>(scalar_ , memType_);

    FuncAllocate::release_value<real>(T_,   memType_);

    FuncAllocate::release_value<real>(sgs_vis_,   memType_);
}


void  ValueNS::
fill(const int nn_max, const int n_scalars)
{
    const real  tmp_vel = 0.0;
    const real  tmp_rho = 1.0;
    const real  tmp_scalar = 0.0;
    const real  tmp_T      = 0.0;
    const real  tmp_sgs_vis = 0.0;

    FuncAllocate::fill_value<real>(u_, tmp_vel, nn_max, memType_);
    FuncAllocate::fill_value<real>(v_, tmp_vel, nn_max, memType_);
    FuncAllocate::fill_value<real>(w_, tmp_vel, nn_max, memType_);

    FuncAllocate::fill_value<real>(rho_, tmp_rho, nn_max, memType_);

    FuncAllocate::fill_value<real>(scalar_ , tmp_scalar, nn_max * n_scalars, memType_);

    FuncAllocate::fill_value<real>(T_,   tmp_T, nn_max, memType_);

    FuncAllocate::fill_value<real>(sgs_vis_,   tmp_sgs_vis, nn_max, memType_);
}
