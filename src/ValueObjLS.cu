#include "ValueObjLS.h"
#include "FuncMath.h"
#include "FuncAllocate.h"
#include "FuncObj.h"
#include "defineFluidProperty.h"


// public //
void  ValueObjLS::
init(const int  nn_max, const enum MemType  memType)
{
    memType_ = memType;

    allocate(nn_max);
    fill(nn_max);
}


void  ValueObjLS::
copy (const int  nn_max, const ValueObjLS&  other)
{
    // solid wall //
    FuncMath::copy_array( lv_obj(),  other.lv_obj(),  nn_max );

    FuncMath::copy_array( rho_obj(),    other.rho_obj(),    nn_max );
    FuncMath::copy_array( u_obj(),      other.u_obj(),      nn_max );
    FuncMath::copy_array( v_obj(),      other.v_obj(),      nn_max );
    FuncMath::copy_array( w_obj(),      other.w_obj(),      nn_max );
    FuncMath::copy_array( T_obj(),      other.T_obj(),      nn_max );
    FuncMath::copy_array( scalar_obj(), other.scalar_obj(), nn_max );

    FuncMath::copy_array( rhon_obj(),    other.rhon_obj(),    nn_max );
    FuncMath::copy_array( un_obj(),      other.un_obj(),      nn_max );
    FuncMath::copy_array( vn_obj(),      other.vn_obj(),      nn_max );
    FuncMath::copy_array( wn_obj(),      other.wn_obj(),      nn_max );
    FuncMath::copy_array( Tn_obj(),      other.Tn_obj(),      nn_max );
    FuncMath::copy_array( scalarn_obj(), other.scalarn_obj(), nn_max );

    // inflow & outflow bc //
//    FuncMath::copy_array( bcTypes_scaler(),  other.bcTypes_scaler(),  nn_max );
//    FuncMath::copy_array( bcTypes_vector(),  other.bcTypes_vector(),  nn_max );
    FuncMath::copy_array( bcTypes_f(),  other.bcTypes_f(),  nn_max );
    FuncMath::copy_array( bcTypes_T(),  other.bcTypes_T(),  nn_max );

//    FuncMath::copy_array( bc_rho(),  other.bc_rho(),  nn_max );
//    FuncMath::copy_array( bc_u(),    other.bc_u(),    nn_max );
//    FuncMath::copy_array( bc_v(),    other.bc_v(),    nn_max );
//    FuncMath::copy_array( bc_w(),    other.bc_w(),    nn_max );
//
//    FuncMath::copy_array( bc_weight(),  other.bc_weight(),  nn_max );

    FuncMath::copy_array( dirichlet_weight(),  other.dirichlet_weight(),  nn_max );
    FuncMath::copy_array( viscosity_weight(),  other.viscosity_weight(),  nn_max );

    // source //
    FuncMath::copy_array( sc_scalar(),  other.sc_scalar(),  nn_max );
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
}


// private //
void  ValueObjLS::
allocate(const int nn_max)
{
    // solid wall //
    FuncAllocate::allocate_value<real>(&lv_obj_,  nn_max, memType_);

    FuncAllocate::allocate_value<real>(&rho_obj_,    nn_max, memType_);
    FuncAllocate::allocate_value<real>(&u_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&v_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&w_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&T_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&scalar_obj_, nn_max, memType_);

    FuncAllocate::allocate_value<real>(&rhon_obj_,    nn_max, memType_);
    FuncAllocate::allocate_value<real>(&un_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&vn_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&wn_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&Tn_obj_,      nn_max, memType_);
    FuncAllocate::allocate_value<real>(&scalarn_obj_, nn_max, memType_);

    // inflow & outflow bc //
//    FuncAllocate::allocate_value<int>( &bcTypes_scaler_,  nn_max, memType_);
//    FuncAllocate::allocate_value<int>( &bcTypes_vector_,  nn_max, memType_);
    FuncAllocate::allocate_value<int>( &bcTypes_f_,       nn_max, memType_);
    FuncAllocate::allocate_value<int>( &bcTypes_T_,       nn_max, memType_);

//    FuncAllocate::allocate_value<real>( &bc_rho_,  nn_max, memType_);
//    FuncAllocate::allocate_value<real>( &bc_u_,    nn_max, memType_);
//    FuncAllocate::allocate_value<real>( &bc_v_,    nn_max, memType_);
//    FuncAllocate::allocate_value<real>( &bc_w_,    nn_max, memType_);
//
//    FuncAllocate::allocate_value<real>( &bc_weight_,  nn_max, memType_);

    FuncAllocate::allocate_value<real>( &dirichlet_weight_,  nn_max, memType_);
    FuncAllocate::allocate_value<real>( &viscosity_weight_,  nn_max, memType_);

    // source //
    FuncAllocate::allocate_value<real>( &sc_scalar_,  nn_max, memType_);
}


void  ValueObjLS::
release ()
{
    // solid wall //
    FuncAllocate::release_value<real>(lv_obj_,  memType_);

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

    // inflow & outflow bc //
//    FuncAllocate::release_value<int>( bcTypes_scaler_, memType_);
//    FuncAllocate::release_value<int>( bcTypes_vector_, memType_);
    FuncAllocate::release_value<int>( bcTypes_f_,      memType_);
    FuncAllocate::release_value<int>( bcTypes_T_,      memType_);

//    FuncAllocate::release_value<real>( bc_rho_,  memType_);
//    FuncAllocate::release_value<real>( bc_u_,    memType_);
//    FuncAllocate::release_value<real>( bc_v_,    memType_);
//    FuncAllocate::release_value<real>( bc_w_,    memType_);
//
//    FuncAllocate::release_value<real>( bc_weight_, memType_);

    FuncAllocate::release_value<real>( dirichlet_weight_, memType_);
    FuncAllocate::release_value<real>( viscosity_weight_, memType_);

    // source //
    FuncAllocate::release_value<real>( sc_scalar_, memType_);
}


void  ValueObjLS::
fill(const int nn_max)
{
    constexpr real T0 = fluid_property::Temperature_bouyancy;

    // solid wall //
    FuncAllocate::fill_value<real>(lv_obj_, FuncObj::lv_fluid()*10.0, nn_max, memType_);

    FuncAllocate::fill_value<real>(rho_obj_,    1.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(u_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(v_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(w_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(T_obj_,      T0 , nn_max, memType_);
    FuncAllocate::fill_value<real>(scalar_obj_, 0.0, nn_max, memType_);

    FuncAllocate::fill_value<real>(rhon_obj_,    1.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(un_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(vn_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(wn_obj_,      0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>(Tn_obj_,      T0 , nn_max, memType_);
    FuncAllocate::fill_value<real>(scalarn_obj_, 0.0, nn_max, memType_);

    // inflow & outflow bc //
//    FuncAllocate::fill_value<int>( bcTypes_scaler_, BCTypes::CalRegion, nn_max, memType_);
//    FuncAllocate::fill_value<int>( bcTypes_vector_, BCTypes::CalRegion, nn_max, memType_);
    FuncAllocate::fill_value<int>( bcTypes_f_,      BCTypes::CalRegion, nn_max, memType_);
    FuncAllocate::fill_value<int>( bcTypes_T_,      BCTypes::CalRegion, nn_max, memType_);

//    FuncAllocate::fill_value<real>( bc_rho_,  1.0, nn_max, memType_);
//    FuncAllocate::fill_value<real>( bc_u_,    0.0, nn_max, memType_);
//    FuncAllocate::fill_value<real>( bc_v_,    0.0, nn_max, memType_);
//    FuncAllocate::fill_value<real>( bc_w_,    0.0, nn_max, memType_);
//
//    FuncAllocate::fill_value<real>( bc_weight_, 0.0, nn_max, memType_);

    FuncAllocate::fill_value<real>( dirichlet_weight_, 0.0, nn_max, memType_);
    FuncAllocate::fill_value<real>( viscosity_weight_, 0.0, nn_max, memType_);

    // source //
    FuncAllocate::fill_value<real>( sc_scalar_, 0.0, nn_max, memType_);
}
