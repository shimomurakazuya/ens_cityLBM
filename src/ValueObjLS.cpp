#include "ValueObjLS.h"
#include "FuncMath.h"
#include "FuncAllocate.h"
#include "FuncObj.h"
#include "defineFluidProperty.h"


// public //
void  ValueObjLS::
init(const int  nn_max, const int n_scalars, const enum MemType  memType)
{
    memType_ = memType;

    allocate(nn_max, n_scalars);
    fill(nn_max, n_scalars);
}


void  ValueObjLS::
copy (const int  nn_max, const int n_scalars, const ValueObjLS&  other)
{
    // solid wall //
    FuncAllocate::copy_values( lv_obj(),  other.lv_obj(),  nn_max, memType_ );
    FuncAllocate::copy_values( pad_obj(), other.pad_obj(),  nn_max, memType_ );

    FuncAllocate::copy_values( rho_obj(),    other.rho_obj(),    nn_max, memType_ );
    FuncAllocate::copy_values( u_obj(),      other.u_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( v_obj(),      other.v_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( w_obj(),      other.w_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( T_obj(),      other.T_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( scalar_obj(), other.scalar_obj(), nn_max * n_scalars, memType_ );

    FuncAllocate::copy_values( rhon_obj(),    other.rhon_obj(),    nn_max, memType_ );
    FuncAllocate::copy_values( un_obj(),      other.un_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( vn_obj(),      other.vn_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( wn_obj(),      other.wn_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( Tn_obj(),      other.Tn_obj(),      nn_max, memType_ );
    FuncAllocate::copy_values( scalarn_obj(), other.scalarn_obj(), nn_max * n_scalars, memType_ );

    FuncAllocate::copy_values( hflux_z_obj(),    other.hflux_z_obj(),    nn_max, memType_ );
    FuncAllocate::copy_values( hfluxn_z_obj(),   other.hfluxn_z_obj(),   nn_max, memType_ );

    // inflow & outflow bc //
//  FuncAllocate::copy_values( bcTypes_scaler(),  other.bcTypes_scaler(),  nn_max, memType_ );
//  FuncAllocate::copy_values( bcTypes_vector(),  other.bcTypes_vector(),  nn_max, memType_ );
    FuncAllocate::copy_values( bcTypes_f(),  other.bcTypes_f(),  nn_max, memType_ );

//  FuncAllocate::copy_values( bc_rho(),  other.bc_rho(),  nn_max, memType_ );
//  FuncAllocate::copy_values( bc_u(),    other.bc_u(),    nn_max, memType_ );
//  FuncAllocate::copy_values( bc_v(),    other.bc_v(),    nn_max, memType_ );
//  FuncAllocate::copy_values( bc_w(),    other.bc_w(),    nn_max, memType_ );
//
//  FuncAllocate::copy_values( bc_weight(),  other.bc_weight(),  nn_max, memType_ );

    FuncAllocate::copy_values( dirichlet_weight(),  other.dirichlet_weight(),  nn_max, memType_ );

    // source //
    FuncAllocate::copy_values( sc_scalar (),  other.sc_scalar (),  nn_max * n_scalars, memType_ );
}


void  ValueObjLS::
swap_objs ()
{
    FuncMath::swap(&rhon_obj_,    &rho_obj_);
    FuncMath::swap(&un_obj_,      &u_obj_  );
    FuncMath::swap(&vn_obj_,      &v_obj_  );
    FuncMath::swap(&wn_obj_,      &w_obj_  );
    FuncMath::swap(&Tn_obj_,      &T_obj_  );
    FuncMath::swap(&scalarn_obj_, &scalar_obj_);

    FuncMath::swap(&hfluxn_z_obj_, &hflux_z_obj_);
}


// private //
void  ValueObjLS::
allocate(const int nn_max, const int n_scalars)
{
    // solid wall //
    FuncAllocate::allocate_value<real>(&lv_obj_,  nn_max, memType_);
    FuncAllocate::allocate_value<real>(&pad_obj_, nn_max, memType_);

    FuncAllocate::allocate_value<real>(&rho_obj_,    nn_max, memType_);
    FuncAllocate::allocate_value<real>(&u_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&v_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&w_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&T_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&scalar_obj_, nn_max * n_scalars, memType_);

    FuncAllocate::allocate_value<real>(&rhon_obj_,    nn_max, memType_);
    FuncAllocate::allocate_value<real>(&un_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&vn_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&wn_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&Tn_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&scalarn_obj_, nn_max * n_scalars, memType_);

    FuncAllocate::allocate_value<real>(&hflux_z_obj_,    nn_max, memType_);
    FuncAllocate::allocate_value<real>(&hfluxn_z_obj_,   nn_max, memType_);

    // inflow & outflow bc //
//    FuncAllocate::allocate_value<int>( &bcTypes_scaler_,  nn_max, memType_);
//    FuncAllocate::allocate_value<int>( &bcTypes_vector_,  nn_max, memType_);
    FuncAllocate::allocate_value<int>( &bcTypes_f_,       nn_max, memType_);

//    FuncAllocate::allocate_value<real>( &bc_rho_,  nn_max, memType_);
//    FuncAllocate::allocate_value<real>( &bc_u_,    nn_max, memType_);
//    FuncAllocate::allocate_value<real>( &bc_v_,    nn_max, memType_);
//    FuncAllocate::allocate_value<real>( &bc_w_,    nn_max, memType_);
//
//    FuncAllocate::allocate_value<real>( &bc_weight_,  nn_max, memType_);

    FuncAllocate::allocate_value<real>( &dirichlet_weight_,  nn_max, memType_);

    // source //
    FuncAllocate::allocate_value<real>( &sc_scalar_ ,  nn_max * n_scalars, memType_);
}


void  ValueObjLS::
release ()
{
    // solid wall //
    FuncAllocate::release_value<real>(lv_obj_,  memType_);
    FuncAllocate::release_value<real>(pad_obj_, memType_);

    FuncAllocate::release_value<real>(rho_obj_,    memType_);
    FuncAllocate::release_value<real>(u_obj_,      memType_);
    FuncAllocate::release_value<real>(v_obj_,      memType_);
    FuncAllocate::release_value<real>(w_obj_,      memType_);
    FuncAllocate::release_value<real>(T_obj_,      memType_);
    FuncAllocate::release_value<real>(scalar_obj_, memType_);

    FuncAllocate::release_value<real>(rhon_obj_,    memType_);
    FuncAllocate::release_value<real>(un_obj_,      memType_);
    FuncAllocate::release_value<real>(vn_obj_,      memType_);
    FuncAllocate::release_value<real>(wn_obj_,      memType_);
    FuncAllocate::release_value<real>(Tn_obj_,      memType_);
    FuncAllocate::release_value<real>(scalarn_obj_, memType_);

    FuncAllocate::release_value<real>(hflux_z_obj_,    memType_);
    FuncAllocate::release_value<real>(hfluxn_z_obj_,   memType_);

    // inflow & outflow bc //
//    FuncAllocate::release_value<int>( bcTypes_scaler_, memType_);
//    FuncAllocate::release_value<int>( bcTypes_vector_, memType_);
    FuncAllocate::release_value<int>( bcTypes_f_,      memType_);

//    FuncAllocate::release_value<real>( bc_rho_,  memType_);
//    FuncAllocate::release_value<real>( bc_u_,    memType_);
//    FuncAllocate::release_value<real>( bc_v_,    memType_);
//    FuncAllocate::release_value<real>( bc_w_,    memType_);
//
//    FuncAllocate::release_value<real>( bc_weight_, memType_);

    FuncAllocate::release_value<real>( dirichlet_weight_, memType_);

    // source //
    FuncAllocate::release_value<real>( sc_scalar_ , memType_);
}


void  ValueObjLS::
fill(const int nn_max, const int n_scalars)
{
    constexpr real T0 = fluid_property::Temperature_bouyancy;

    // solid wall //
    FuncAllocate::fill_value<real>(lv_obj_,  FuncObj::lv_fluid()*10.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(pad_obj_, 0.0, nn_max, memType_);

    FuncAllocate::fill_value<real>(rho_obj_,    1.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(u_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(v_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(w_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(T_obj_,      T0 , nn_max, memType_);
    FuncAllocate::fill_value<real>(scalar_obj_, 0.0, nn_max * n_scalars, memType_);

    FuncAllocate::fill_value<real>(rhon_obj_,    1.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(un_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(vn_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(wn_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(Tn_obj_,      T0 , nn_max, memType_);
    FuncAllocate::fill_value<real>(scalarn_obj_, 0.0, nn_max * n_scalars, memType_);

    FuncAllocate::fill_value<real>(hflux_z_obj_,    0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(hfluxn_z_obj_,   0.0, nn_max, memType_);

    // inflow & outflow bc //
//    FuncAllocate::fill_value<int>( bcTypes_scaler_, BCTypes::CalRegion, nn_max, memType_);
//    FuncAllocate::fill_value<int>( bcTypes_vector_, BCTypes::CalRegion, nn_max, memType_);
    FuncAllocate::fill_value<int>( bcTypes_f_,      BCTypes::CalRegion, nn_max, memType_);

//    FuncAllocate::fill_value<real>( bc_rho_,  1.0, nn_max, memType_);
//    FuncAllocate::fill_value<real>( bc_u_,    0.0, nn_max, memType_);
//    FuncAllocate::fill_value<real>( bc_v_,    0.0, nn_max, memType_);
//    FuncAllocate::fill_value<real>( bc_w_,    0.0, nn_max, memType_);
//
//    FuncAllocate::fill_value<real>( bc_weight_, 0.0, nn_max, memType_);

    FuncAllocate::fill_value<real>( dirichlet_weight_, 0.0, nn_max, memType_);

    // source //
    FuncAllocate::fill_value<real>( sc_scalar_ , 0.0, nn_max * n_scalars, memType_);
}
