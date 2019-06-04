#include "MonitorData.h"


void MonitorData::
CoordinatesToCell(
    const Tree& tree,
    const MeshValue* meshValues,
    const int rank
    )
{
//    if  (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf = tree.number_of_nodes();
    const int  nx = DefAMR::NX_LEAF;
    const int  ny = DefAMR::NX_LEAF;
    const int  nz = DefAMR::NX_LEAF;

    int searchflg = 0;
    for (int l=0; l<n_leaf; l++) {
        if ( tree.nodes(l)->nodeCalFlags().Cal() == false ) { continue; }

        const int  lv = tree.nodes(l)->level();

        const MeshValue&    meshValue = meshValues[lv];
        const Array3D<int>  offsets   = tree.nodes(l)->neighbor_mesh_offsets();

        const int  offset     = offsets.offset0();

        // coordinates //
        const real  dx0 = (meshValue.coordinates().x(offset + 1    ) - meshValue.coordinates().x(offset));
        const real  dy0 = (meshValue.coordinates().y(offset + nx   ) - meshValue.coordinates().y(offset));
        const real  dz0 = (meshValue.coordinates().z(offset + nx*ny) - meshValue.coordinates().z(offset));

        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real x_leaf[2] = { x[0], x[0] + dx0*DefAMR::NX_LEAF };
        const real y_leaf[2] = { y[0], y[0] + dy0*DefAMR::NX_LEAF };
        const real z_leaf[3] = { z[0], z[0] + dz0*DefAMR::NX_LEAF };

        if ( x_ < x_leaf[0] || x_ > x_leaf[1]
          || y_ < y_leaf[0] || y_ > y_leaf[1]
          || z_ < z_leaf[0] || z_ > z_leaf[1] ) { continue; }

        // xyz_vtu //
        for (int k=0; k<nz; k++) { // leftside index //
            for (int j=0; j<ny; j++) {
                for (int i=0; i<nx; i++) {
//        for (int k=0; k<nz+1; k++) {
//            for (int j=0; j<ny+1; j++) {
//                for (int i=0; i<nx+1; i++) {
                    const float cx = meshValue.coordinates().x(offset) + dx0*i;
                    const float cy = meshValue.coordinates().y(offset) + dy0*j;
                    const float cz = meshValue.coordinates().z(offset) + dz0*k;

                    if ( ( cx <= x_ && x_ < cx+dx0) && (cy <= y_ && y_ < cy +dy0 ) && ( cz <= z_ && z_ < cz+dz0 )){
                        l_ = l;
                        i_ = i;
                        j_ = j;
                        k_ = k;

                        searchflg = 1;
                        break;
                    }

                }
                if (searchflg == 1) break;
            }
            if (searchflg == 1) break;
        }

        if (searchflg == 1) {
            const float x0 = meshValue.coordinates().x(offset)+dx0*i_;
            const float y0 = meshValue.coordinates().y(offset)+dy0*j_;
            const float z0 = meshValue.coordinates().z(offset)+dz0*k_;

            dist_coef_[0] = (x_-x0)/dx0;
            dist_coef_[1] = (y_-y0)/dy0;
            dist_coef_[2] = (z_-z0)/dz0;

            iscalc_ = 1;
            init_rank_ = rank;

            break;
        }
    }

}


void MonitorData::
Approximatecalculation(
    const Tree&          tree,
    const Parameters&    parameters,
    const MeshValue*     meshValues
    )
{
    const int  lv = tree.nodes(l_)->level();
    const MeshValue& meshValue = meshValues[lv];

    const Array3D<int>  offsets = tree.nodes(l_)->neighbor_mesh_offsets();
    const int  offset  = offsets.offset0();

    const int  nx_leaf = DefAMR::NX_LEAF;
    const real c_ref = parameters.c_ref_lbm();

    const int calcflg = tree.nodes(l_)->nodeCalFlags().Cal() ? 1: 0;

    u_ = 0.0f;
    v_ = 0.0f;
    w_ = 0.0f;
    levelset_obj_ = 0.0f;
    scalar_       = 0.0f;
    T_            = 0.0f;

    if (calcflg == 0) {
        iscalc_ = calcflg;
        return;
    }
    for (int k=0; k<2; k++) {
        for (int j=0; j<2; j++) {
            for (int i=0; i<2; i++) {
                float coef0 = dist_coef_[0]*i+((i+1)%2)*(1.0f-dist_coef_[0]);
                float coef1 = dist_coef_[1]*j+((j+1)%2)*(1.0f-dist_coef_[1]);
                float coef2 = dist_coef_[2]*k+((k+1)%2)*(1.0f-dist_coef_[2]);


                u_ += (coef0*coef1*coef2) * FuncAMRMesh::RawData( meshValue.valueNS().u(), i_+i,j_+j,k_+k, nx_leaf, offsets )*c_ref;
                v_ += (coef0*coef1*coef2) * FuncAMRMesh::RawData( meshValue.valueNS().v(), i_+i,j_+j,k_+k, nx_leaf, offsets )*c_ref;
                w_ += (coef0*coef1*coef2) * FuncAMRMesh::RawData( meshValue.valueNS().w(), i_+i,j_+j,k_+k, nx_leaf, offsets )*c_ref;

                levelset_obj_ += (coef0*coef1*coef2) * FuncAMRMesh::RawData( meshValue.valueObjLS().lv_obj(), i_+i,j_+j,k_+k, nx_leaf, offsets );

                scalar_ += (coef0*coef1*coef2) * FuncAMRMesh::RawData( meshValue.valueNS().scalar(), i_+i,j_+j,k_+k, nx_leaf, offsets );
                T_      += (coef0*coef1*coef2) * FuncAMRMesh::RawData( meshValue.valueNS().T(),      i_+i,j_+j,k_+k, nx_leaf, offsets );
            }
        }
    }

}
