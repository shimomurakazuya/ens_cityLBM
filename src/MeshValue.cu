#include "MeshValue.h"


// public //
void  MeshValue::
init(const int  nn_max)
{
    nn_max_  = nn_max;

    // coordinates xyz //
    init_coordinates(nn_max);

    // values //
    init_mesh_values(nn_max);
}


void  MeshValue::
set_coordinate(const int  lv, const Grid& grid, const Tree& tree)
{
    // coordinates xyz //
    coordinates_.set_uniform(lv, grid, tree);
}


void  MeshValue::
copy_MeshValue(
    const MeshValue&    other,
    const Grid&         grid,
    const Tree&         tree,
    const bool          is_reset
    )
{
    nn_max_ = other.nn_max();
//    std::cout << "nn_max_ = " << nn_max_ << std::endl;

    if (is_reset) {
        // coordinates xyz //
        coordinates().reallocate( other.coordinates().memType(), other.nn_max() );

        // values //
        valueObjLS(). reallocate( other.valueObjLS(). memType(), other.nn_max() );
        valueNS().    reallocate( other.valueNS().    memType(), other.nn_max() );
        valueLBM().   reallocate( other.valueLBM().   memType(), other.nn_max() );
        valueBuff().  reallocate( other.valueBuff().  memType(), other.nn_max() );
    }

    // coordinates xyz //
    coordinates().copy( other.nn_max(), other.coordinates() );

    // values //
    valueObjLS(). copy( other.nn_max(), other.valueObjLS() );
    valueNS().    copy( other.nn_max(), other.valueNS() );
    valueLBM().   copy( other.nn_max(), other.valueLBM() );
//    valueBuff().  copy( other.nn_max(), other.valueBuff() );
}


// private //
void  MeshValue::
init_coordinates(const int  nn_max)
{
    coordinates_.init(nn_max, MemType::Managed);
}


void  MeshValue::
init_mesh_values(const int  nn_max)
{
    valueObjLS_.init(nn_max, memType_);
    valueNS_   .init(nn_max, memType_);
    valueLBM_  .init(nn_max, memType_);

    valueBuff_ .init(nn_max, memType_);
}
