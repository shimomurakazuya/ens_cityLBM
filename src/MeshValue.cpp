#include "MeshValue.h"
#include "FuncObj.h"


// public //
void  MeshValue::
init(const int  nn_max, const int n_scalars)
{
    nn_max_    = nn_max;
    n_scalars_ = n_scalars;

    // coordinates xyz //
    init_coordinates(nn_max);

    // values //
    init_mesh_values(nn_max, n_scalars);
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
    n_scalars_ = other.n_scalars();
//    std::cout << "nn_max_ = " << nn_max_ << std::endl;

    if (is_reset) {
        // coordinates xyz //
        coordinates().reallocate( other.coordinates().memType(), other.nn_max() );

        // values //
        valueObjLS(). reallocate( other.valueObjLS(). memType(), other.nn_max(), other.n_scalars() );
        valueNS().    reallocate( other.valueNS().    memType(), other.nn_max(), other.n_scalars() );
        valueLBM().   reallocate( other.valueLBM().   memType(), other.nn_max() );
    }

    // coordinates xyz //
    coordinates().copy( other.nn_max(), other.coordinates() );

    // values //
    valueObjLS(). copy( other.nn_max(), other.n_scalars(), other.valueObjLS() );
    valueNS().    copy( other.nn_max(), other.n_scalars(), other.valueNS() );
    valueLBM().   copy( other.nn_max(), other.valueLBM() );
}

void MeshValue::
mask_meshValue_wo_halo(const Tree& tree, const int lv) {
    for(const auto& node: tree.nodes()) {
        if( node.level() == lv && (!node.nodeCalFlags().Cal()) ) { // if halo
            const auto offset = node.neighbor_mesh_offsets().offset0();
            std::vector<real*> real_ptrs_all {
                valueObjLS_.lv_obj(), valueObjLS_.pad_obj(),
                valueObjLS_.rho_obj(), valueObjLS_.rhon_obj(),
                valueObjLS_.u_obj(), valueObjLS_.un_obj(), 
                valueObjLS_.v_obj(), valueObjLS_.vn_obj(),
                valueObjLS_.w_obj(), valueObjLS_.wn_obj(),
                valueObjLS_.T_obj(), valueObjLS_.Tn_obj(),
                valueObjLS_.scalar_obj(), valueObjLS_.scalarn_obj(),
                valueObjLS_.hflux_z_obj(), 
                valueObjLS_.dirichlet_weight(),
                valueNS_.rho(),
                valueNS_.u(), valueNS_.v(), valueNS_.w(), 
                valueNS_.T(), valueNS_.scalar(), valueNS_.sgs_vis()
            };
            std::vector<real*> real_lbm_ptrs_all {
                valueLBM_.f_lbm()
            };
            std::vector<int*> int_ptrs_all {
                //coordinates_.i(), coordinates_.j(), coordinates_.k(), // ignore coordinates
                //valueObjLS_.bcTypes_f(), valueObjLS_.bcTypes_T() // ignore bcTypes (∵ init only, no comm)
            };
            for(int ijk=0; ijk<DefAMR::NN_LEAF; ijk++) {
                const auto id = offset + ijk;
                const auto id_lbm = offset*LBM_velocity_model::nQ + ijk;
                for(real* p: real_ptrs_all) { p[id] = 999999.9; }
                for(real* p: real_lbm_ptrs_all) {
                    for(int q=0; q<LBM_velocity_model::nQ; q++) {
                        p[id_lbm + q*DefAMR::NN_LEAF] = 999999.9;
                    }
                }
                for(int* p: int_ptrs_all) { p[id] = 999999; }
            }
        } // if halo
    }
}

void MeshValue::
mask_meshValue_wo_obj(const Tree& tree, const int lv) {
    for(const auto& node: tree.nodes()) {
        if( node.level() == lv && node.nodeCalFlags().Cal() ) { // if obj
            const auto offset = node.neighbor_mesh_offsets().offset0();
            const auto lv_obj = valueObjLS_.lv_obj();
            std::vector<real*> real_ptrs_all {
                valueObjLS_.lv_obj(), valueObjLS_.pad_obj(),
                valueObjLS_.rho_obj(), valueObjLS_.rhon_obj(),
                valueObjLS_.u_obj(), valueObjLS_.un_obj(), 
                valueObjLS_.v_obj(), valueObjLS_.vn_obj(),
                valueObjLS_.w_obj(), valueObjLS_.wn_obj(),
                valueObjLS_.T_obj(), valueObjLS_.Tn_obj(),
                valueObjLS_.scalar_obj(), valueObjLS_.scalarn_obj(),
                valueObjLS_.hflux_z_obj(), 
                valueObjLS_.dirichlet_weight(),
                valueNS_.rho(),
                valueNS_.u(), valueNS_.v(), valueNS_.w(), 
                valueNS_.T(), valueNS_.scalar(), valueNS_.sgs_vis()
            };
            std::vector<real*> real_lbm_ptrs_all {
                valueLBM_.f_lbm()
            };
            std::vector<int*> int_ptrs_all {
                //coordinates_.i(), coordinates_.j(), coordinates_.k(), // ignore coordinates
                //valueObjLS_.bcTypes_f(), valueObjLS_.bcTypes_T() // ignore bcTypes (∵ init only, no comm)
            };
            for(int ijk=0; ijk<DefAMR::NN_LEAF; ijk++) {
                const auto id = offset + ijk;
                const auto id_lbm = offset*LBM_velocity_model::nQ + ijk;
                if(FuncObj::is_obj(lv_obj[id])) {
                    for(real* p: real_ptrs_all) { p[id] = 999999.9; }
                    for(real* p: real_lbm_ptrs_all) {
                        for(int q=0; q<LBM_velocity_model::nQ; q++) {
                            p[id_lbm + q*DefAMR::NN_LEAF] = 999999.9;
                        }
                    }
                    for(int* p: int_ptrs_all) { p[id] = 999999; }
                }
            }
        } // if obj
    }
}

// private //
void  MeshValue::
init_coordinates(const int  nn_max)
{
    coordinates_.init(nn_max, memType_);
    //coordinates_.init(nn_max, MemType::Managed);
}


void  MeshValue::
init_mesh_values(const int  nn_max, const int n_scalars)
{
    valueObjLS_.init(nn_max, n_scalars, memType_);
    valueNS_   .init(nn_max, n_scalars, memType_);
    valueLBM_  .init(nn_max, memType_);
}
