#include "BoundaryConditions.h"

#include <sstream>
#include <iomanip>
#include "FuncObj.h"
#include "FuncThermalConvection.h"


void  BoundaryConditions::InitializeByWRFData(
    const int           time_num,
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << std::endl << __PRETTY_FUNCTION__ << " : time = " << time_num << " minutes" << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

    const int  nx_da = 16;
    const int  ny_da = 16;
    const int  nz_da = 50;

    const real pitch_da[3] = { 500.0, 500.0, 50.0 };

    float* u_da = new float[nx_da*ny_da*nz_da];
    float* v_da = new float[nx_da*ny_da*nz_da];
    float* w_da = new float[nx_da*ny_da*nz_da];
    float* T_da = new float[nx_da*ny_da*nz_da];
    float* E_da = new float[nx_da*ny_da*nz_da];

    ReadOklahoma_uvwTE(u_da, v_da, w_da, T_da, E_da,
                       nx_da, ny_da, nz_da,
                       time_num);

//    auto fx_da = [nx_da, ny_da, nz_da](const float* f, int id_da){ return (f[id_da+1          ] - f[id_da-1          ])*(float)0.5; };
//    auto fy_da = [nx_da, ny_da, nz_da](const float* f, int id_da){ return (f[id_da+nx_da      ] - f[id_da-nx_da      ])*(float)0.5; };
//    auto fz_da = [nx_da, ny_da, nz_da](const float* f, int id_da){ return (f[id_da+nx_da*ny_da] - f[id_da-nx_da*ny_da])*(float)0.5; };

//    auto fx_up_da = [nx_da, ny_da, nz_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+1          ] - f[id_da] :  f[id_da] - f[id_da-1          ]; };
//    auto fy_up_da = [nx_da, ny_da, nz_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+nx_da      ] - f[id_da] :  f[id_da] - f[id_da-nx_da      ]; };
//    auto fz_up_da = [nx_da, ny_da, nz_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+nx_da*ny_da] - f[id_da] :  f[id_da] - f[id_da-nx_da*ny_da]; };

//#pragma omp parallel for schedule(dynamic)
#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real   dx = meshValues[lv].coordinates().dx();
//        const real   dt  = dx/c_ref;


        // check position from map offsets //
        const real dx_map    =     4.0;
        const real x_min_map = -4000.0;
        const real y_min_map = -4000.0;
        const real z_min_map =     0.0;
        const real x_leaf[2] = { x[0] - x_min_map, x[0] + dx*DefAMR::NX_LEAF - x_min_map };
        const real y_leaf[2] = { y[0] - y_min_map, y[0] + dx*DefAMR::NX_LEAF - y_min_map };
        const real z_leaf[3] = { z[0] - z_min_map, z[0] + dx*DefAMR::NX_LEAF*(real)0.5 - z_min_map, z[0] + dx*DefAMR::NX_LEAF - z_min_map };


        std::vector<int>  ids_da;
        for (int kk=0; kk<3; kk++) {
        for (int jj=0; jj<2; jj++) {
        for (int ii=0; ii<2; ii++) {
            const int  i_da = x_leaf[ii]/pitch_da[0];
            const int  j_da = y_leaf[jj]/pitch_da[1];
            const int  k_da = z_leaf[kk]/pitch_da[2];

            const bool is_cal_region = ( i_da>=0 && i_da <=nx_da-1 && j_da>=0 && j_da <=ny_da-1 && k_da>=0 && k_da <=nz_da-1 ) ? true : false;

            if (is_cal_region) { ids_da.push_back( i_da + j_da*nx_da + k_da*nx_da*ny_da ); }
        }
        }
        }
        std::sort(ids_da.begin(), ids_da.end());
        ids_da.erase(std::unique(ids_da.begin(), ids_da.end()), ids_da.end());


        for (const int id_da : ids_da) {
            const int   i_da = (id_da%(nx_da*ny_da))%nx_da;
            const int   j_da = (id_da%(nx_da*ny_da))/nx_da;
            const int   k_da =  id_da/(nx_da*ny_da);


            // box //
            const real width_da[3] = { pitch_da[0] + (real)2.0*dx_map, pitch_da[1] + (real)2.0*dx_map, pitch_da[2] + (real)2.0*dx_map };

            const real center_da[3] = { (i_da + (real)0.5)*pitch_da[0] + x_min_map,
                                        (j_da + (real)0.5)*pitch_da[1] + y_min_map,
                                        (k_da + (real)0.5)*pitch_da[2] + z_min_map };

            const real box_min[3] = { center_da[0] - (real)0.5*width_da[0],
                                      center_da[1] - (real)0.5*width_da[1],
                                      center_da[2] - (real)0.5*width_da[2] };

            const real box_max[3] = { center_da[0] + (real)0.5*width_da[0],
                                      center_da[1] + (real)0.5*width_da[1],
                                      center_da[2] + (real)0.5*width_da[2] };


            const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];

            real*  rho   = &meshValues[lv].valueNS().rho()[offset];
            real*  u     = &meshValues[lv].valueNS().u()  [offset];
            real*  v     = &meshValues[lv].valueNS().v()  [offset];
            real*  w     = &meshValues[lv].valueNS().w()  [offset];
            real*  T     = &meshValues[lv].valueNS().T()  [offset];

            for (int id=0; id<nn_cell; id++) {
                // lv //
                const real x_tmp = x[id] + (real)0.5*dx;
                const real y_tmp = y[id] + (real)0.5*dx;
                const real z_tmp = z[id] + (real)0.5*dx;

                const real xyz[3] = { x_tmp, y_tmp, z_tmp };

                const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);

                // in the box area //
                if ( FuncObj::is_fluid(lv_obs[id]) && FuncObj::is_obj(box_lv_obs) ) {
                    const real ox = (x_tmp - center_da[0])/pitch_da[0];
                    const real oy = (y_tmp - center_da[1])/pitch_da[1];
                    const real oz = (z_tmp - center_da[2])/pitch_da[2];

                    // interpolation //
                    const real u_da_int = f_interpolate3d(u_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);
                    const real v_da_int = f_interpolate3d(v_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);
                    const real w_da_int = f_interpolate3d(w_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);
                    const real T_da_int = f_interpolate3d(T_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);

//                    if (T_da_int < 273.0 || T_da_int > 313.0) { std::cout << "\n----------\nT_da_int = " << T_da_int << "\n----------\n" << std::endl; }


//                    rho[id] = (real)1.0;
                    rho[id] = FuncThermalConvection::rho_buoyancy(T_da_int);
                    T  [id] = T_da_int;

                    u[id] = u_da_int / c_ref;
                    v[id] = v_da_int / c_ref;
                    w[id] = w_da_int / c_ref;
                }
            }


        }

    }

    delete [] u_da;
    delete [] v_da;
    delete [] w_da;
    delete [] T_da;
    delete [] E_da;
}


void  BoundaryConditions::InitializeByGroundData(
    const int           time_num,
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << std::endl << __PRETTY_FUNCTION__ << " : time = " << time_num << " minutes" << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;

    const int  nx_da = 16;
    const int  ny_da = 16;

    const real pitch_da[2] = { 500.0, 500.0 };

    float* T_da = new float[nx_da*ny_da];
    ReadOklahoma_GroundTb(T_da, nx_da, ny_da, time_num);

//    auto fx_up_da = [nx_da, ny_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+1          ] - f[id_da] :  f[id_da] - f[id_da-1          ]; };
//    auto fy_up_da = [nx_da, ny_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+nx_da      ] - f[id_da] :  f[id_da] - f[id_da-nx_da      ]; };

//#pragma omp parallel for schedule(dynamic)
#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real   dx = meshValues[lv].coordinates().dx();


        // check position from map offsets //
//        const real dx_map    =  4.0;
        const real x_min_map = -4000.0;
        const real y_min_map = -4000.0;
        const real x_leaf[2] = { x[0] - x_min_map, x[0] + dx*DefAMR::NX_LEAF - x_min_map };
        const real y_leaf[2] = { y[0] - y_min_map, y[0] + dx*DefAMR::NX_LEAF - y_min_map };


        std::vector<int>  ids_da;
        for (int jj=0; jj<2; jj++) {
        for (int ii=0; ii<2; ii++) {
            const int  i_da = x_leaf[ii]/pitch_da[0];
            const int  j_da = y_leaf[jj]/pitch_da[1];

            const bool is_cal_region = ( i_da>=0 && i_da <=nx_da-1 && j_da>=0 && j_da <=ny_da-1 ) ? true : false;

            if (is_cal_region) { ids_da.push_back( i_da + j_da*nx_da ); }
        }
        }
        std::sort(ids_da.begin(), ids_da.end());
        ids_da.erase(std::unique(ids_da.begin(), ids_da.end()), ids_da.end());


        for (const int id_da : ids_da) {
            const int   i_da = id_da%nx_da;
            const int   j_da = id_da/nx_da;

            // box //
            const real width_da[3] = { pitch_da[0],
                                       pitch_da[1],
                                       fabs(parameters.coefDomain().z_global_domain_max - parameters.coefDomain().z_global_domain_min) };

            const real center_da[3] = { (i_da + (real)0.5)*pitch_da[0] + x_min_map,
                                        (j_da + (real)0.5)*pitch_da[1] + y_min_map,
                                        width_da[2]*(real)0.5 };

            const real box_min[3] = { center_da[0] - (real)0.5*width_da[0],
                                      center_da[1] - (real)0.5*width_da[1],
                                      center_da[2] - (real)1.0*width_da[2] };

            const real box_max[3] = { center_da[0] + (real)0.5*width_da[0],
                                      center_da[1] + (real)0.5*width_da[1],
                                      center_da[2] };
//                                      center_da[2] + (real)0.5*width_da[2] + dx };


            const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];

            real*  T     = &meshValues[lv].valueNS().T()  [offset];

            for (int id=0; id<nn_cell; id++) {
                // lv //
                const real x_tmp = x[id] + (real)0.5*dx;
                const real y_tmp = y[id] + (real)0.5*dx;
                const real z_tmp = z[id] + (real)0.5*dx;

                const real xyz[3] = { x_tmp, y_tmp, z_tmp };

                const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);

                // in the box area //
                if ( FuncObj::is_obj(lv_obs[id]) && FuncObj::is_obj(box_lv_obs) ) {
                    const real ox = (x_tmp - center_da[0])/pitch_da[0];
                    const real oy = (y_tmp - center_da[1])/pitch_da[1];

                    // interpolation //
                    const real T_da_int = f_interpolate2d(T_da, ox,oy, i_da,j_da, nx_da, ny_da);

                    T  [id] = T_da_int;
                }
            }


        }

    }

    delete [] T_da;
}


void  BoundaryConditions::ReadWRFData(
    const int           time_num,
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << std::endl << __PRETTY_FUNCTION__ << " : time = " << time_num << " minutes" << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

    const int  nx_da = 16;
    const int  ny_da = 16;
    const int  nz_da = 50;

    const real pitch_da[3] = { 500.0, 500.0, 50.0 };

    float* u_da = new float[nx_da*ny_da*nz_da];
    float* v_da = new float[nx_da*ny_da*nz_da];
    float* w_da = new float[nx_da*ny_da*nz_da];
    float* T_da = new float[nx_da*ny_da*nz_da];
    float* E_da = new float[nx_da*ny_da*nz_da];

    ReadOklahoma_uvwTE(u_da, v_da, w_da, T_da, E_da,
                       nx_da, ny_da, nz_da,
                       time_num);

//    auto fx_da = [nx_da, ny_da, nz_da](const float* f, int id_da){ return (f[id_da+1          ] - f[id_da-1          ])*(float)0.5; };
//    auto fy_da = [nx_da, ny_da, nz_da](const float* f, int id_da){ return (f[id_da+nx_da      ] - f[id_da-nx_da      ])*(float)0.5; };
//    auto fz_da = [nx_da, ny_da, nz_da](const float* f, int id_da){ return (f[id_da+nx_da*ny_da] - f[id_da-nx_da*ny_da])*(float)0.5; };

//    auto fx_up_da = [nx_da, ny_da, nz_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+1          ] - f[id_da] :  f[id_da] - f[id_da-1          ]; };
//    auto fy_up_da = [nx_da, ny_da, nz_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+nx_da      ] - f[id_da] :  f[id_da] - f[id_da-nx_da      ]; };
//    auto fz_up_da = [nx_da, ny_da, nz_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+nx_da*ny_da] - f[id_da] :  f[id_da] - f[id_da-nx_da*ny_da]; };

#pragma omp parallel for schedule(dynamic)
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real   dx = meshValues[lv].coordinates().dx();
        const real   dt = dx/c_ref;


        // check position from map offsets //
        const real dx_map    =  4.0;
        const real x_min_map = -4000.0;
        const real y_min_map = -4000.0;
        const real z_min_map =     0.0;
        const real x_leaf[2] = { x[0] - x_min_map, x[0] + dx*DefAMR::NX_LEAF - x_min_map };
        const real y_leaf[2] = { y[0] - y_min_map, y[0] + dx*DefAMR::NX_LEAF - y_min_map };
        const real z_leaf[3] = { z[0] - z_min_map, z[0] + dx*DefAMR::NX_LEAF*(real)0.5 - z_min_map, z[0] + dx*DefAMR::NX_LEAF - z_min_map };


        std::vector<int>  ids_da;
        for (int kk=0; kk<3; kk++) {
        for (int jj=0; jj<2; jj++) {
        for (int ii=0; ii<2; ii++) {
            const int  i_da = x_leaf[ii]/pitch_da[0];
            const int  j_da = y_leaf[jj]/pitch_da[1];
            const int  k_da = z_leaf[kk]/pitch_da[2];

            const bool is_cal_region = ( i_da>=0 && i_da <=nx_da-1 && j_da>=0 && j_da <=ny_da-1 && k_da>=0 && k_da <=nz_da-1 ) ? true : false;

            if (is_cal_region) { ids_da.push_back( i_da + j_da*nx_da + k_da*nx_da*ny_da ); }
        }
        }
        }
        std::sort(ids_da.begin(), ids_da.end());
        ids_da.erase(std::unique(ids_da.begin(), ids_da.end()), ids_da.end());


        for (const int id_da : ids_da) {
            const int   i_da = (id_da%(nx_da*ny_da))%nx_da;
            const int   j_da = (id_da%(nx_da*ny_da))/nx_da;
            const int   k_da =  id_da/(nx_da*ny_da);

            const bool is_outer_boundary = (
                                       i_da==0 || i_da==nx_da-1
                                    || j_da==0 || j_da==ny_da-1
                                    || k_da==nz_da-1
                                    || i_da==int( (parameters.coefDomain().x_global_domain_min - x_min_map)/pitch_da[0]   )
                                    || j_da==int( (parameters.coefDomain().y_global_domain_min - y_min_map)/pitch_da[1]   )
                                    || i_da==int( (parameters.coefDomain().x_global_domain_max - x_min_map)/pitch_da[0]   )
                                    || j_da==int( (parameters.coefDomain().y_global_domain_max - y_min_map)/pitch_da[1]   )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]   )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-1 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-2 )
                                    ) ? true : false;

            const bool is_inner_boundary = (
                                       k_da==nz_da-2
                                    || k_da==nz_da-3
                                    || k_da==nz_da-4
                                    || k_da==nz_da-5
                                    || i_da==int( (parameters.coefDomain().x_global_domain_min - x_min_map)/pitch_da[0]+1 )
                                    || j_da==int( (parameters.coefDomain().y_global_domain_min - y_min_map)/pitch_da[1]+1 )
                                    || i_da==int( (parameters.coefDomain().x_global_domain_max - x_min_map)/pitch_da[0]-1 )
                                    || j_da==int( (parameters.coefDomain().y_global_domain_max - y_min_map)/pitch_da[1]-1 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-3 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-4 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-5 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-6 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-7 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-8 )
                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-9 )
                                    ) ? true : false;

            const bool is_boundary = ( is_outer_boundary || is_inner_boundary ) ? true : false;

            constexpr real minutes   = 60.0;
            const real _coef_dirichlet_weight_1min = (real)1.0/(minutes / dt);
            const real _coef_dirichlet_weight = ( is_outer_boundary ) ? _coef_dirichlet_weight_1min * (real)1.0   // coefficient1
                                                                      : _coef_dirichlet_weight_1min / (real)10.0;  // coefficient2


//            const bool is_boundary = ( i_da==0 || i_da==nx_da-1
//                                    || j_da==0 || j_da==ny_da-1
//                                    || k_da==nz_da-1
//                                    || k_da==nz_da-2
//                                    || k_da==nz_da-3
//                                    || k_da==nz_da-4
//                                    || k_da==nz_da-5
//                                    || i_da==int( (parameters.coefDomain().x_global_domain_min - x_min_map)/pitch_da[0]   ) || i_da==int( (parameters.coefDomain().x_global_domain_min - x_min_map)/pitch_da[0]+1 )
//                                    || j_da==int( (parameters.coefDomain().y_global_domain_min - y_min_map)/pitch_da[1]   ) || j_da==int( (parameters.coefDomain().y_global_domain_min - y_min_map)/pitch_da[1]+1 )
//                                    || i_da==int( (parameters.coefDomain().x_global_domain_max - x_min_map)/pitch_da[0]   ) || i_da==int( (parameters.coefDomain().x_global_domain_max - x_min_map)/pitch_da[0]-1 )
//                                    || j_da==int( (parameters.coefDomain().y_global_domain_max - y_min_map)/pitch_da[1]   ) || j_da==int( (parameters.coefDomain().y_global_domain_max - y_min_map)/pitch_da[1]-1 )
//                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]   )
//                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-1 )
//                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-2 )
//                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-3 )
//                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-4 )
//                                    || k_da==int( (parameters.coefDomain().z_global_domain_max - z_min_map)/pitch_da[2]-5 )
//                                    ) ? true : false;

            if (!is_boundary) { continue; }

//            const real coef_wid = is_boundary ? _coef_wid[0] : _coef_wid[1];
            const real coef_dirichlet_weight = _coef_dirichlet_weight > 0.5 ? 0.5 : _coef_dirichlet_weight;
            const real coef_viscosity_weight = 0.0;
//            const real coef_viscosity_weight = coef_dirichlet_weight*0.1;


            // box //
//            const real width_da[3] = { std::fmin(coef_wid*dx, pitch_da[0]) + (real)dx, std::fmin(coef_wid*dx, pitch_da[1]) + (real)dx, std::fmin(coef_wid*dx, pitch_da[2]) + (real)dx };
            const real width_da[3] = { pitch_da[0] + (real)2.5*dx_map, pitch_da[1] + (real)2.5*dx_map, pitch_da[2] + (real)2.5*dx_map };
//            const real width_da_max = std::max(width_da[0], std::max(width_da[1], width_da[2]));

            const real center_da[3] = { (i_da + (real)0.5)*pitch_da[0] + x_min_map,
                                        (j_da + (real)0.5)*pitch_da[1] + y_min_map,
                                        (k_da + (real)0.5)*pitch_da[2] + z_min_map };

            const real box_min[3] = { center_da[0] - (real)0.5*width_da[0],
                                      center_da[1] - (real)0.5*width_da[1],
                                      center_da[2] - (real)0.5*width_da[2] };

            const real box_max[3] = { center_da[0] + (real)0.5*width_da[0],
                                      center_da[1] + (real)0.5*width_da[1],
                                      center_da[2] + (real)0.5*width_da[2] };


            const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];

            int*   bcTypes_f = &meshValues[lv].valueObjLS().bcTypes_f()[offset];

            real*  rho_obs   = &meshValues[lv].valueObjLS().rho_obj()[offset];
            real*  u_obs     = &meshValues[lv].valueObjLS().u_obj()  [offset];
            real*  v_obs     = &meshValues[lv].valueObjLS().v_obj()  [offset];
            real*  w_obs     = &meshValues[lv].valueObjLS().w_obj()  [offset];
            real*  T_obs     = &meshValues[lv].valueObjLS().T_obj()  [offset];

            real*  dirichlet_weight = &meshValues[lv].valueObjLS().dirichlet_weight()[offset];
            real*  viscosity_weight = &meshValues[lv].valueObjLS().viscosity_weight()[offset];

            for (int id=0; id<nn_cell; id++) {
                if ( !FuncObj::is_fluid(lv_obs[id]) ) { continue; }

                // lv //
                const real x_tmp = x[id] + (real)0.5*dx;
                const real y_tmp = y[id] + (real)0.5*dx;
                const real z_tmp = z[id] + (real)0.5*dx;

                const real xyz[3] = { x_tmp, y_tmp, z_tmp };

                const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);

                // in the box area //
                if ( FuncObj::is_included_in_fluid(lv_obs[id], dx, (real)1.0) && FuncObj::is_obj(box_lv_obs) ) {
//                if ( FuncObj::is_fluid(lv_obs[id]) && FuncObj::is_obj(box_lv_obs) ) { // bug : effect on the bounce-back boundary condition //
                    const real ox = (x_tmp - center_da[0])/pitch_da[0];
                    const real oy = (y_tmp - center_da[1])/pitch_da[1];
                    const real oz = (z_tmp - center_da[2])/pitch_da[2];

                    const real u_da_int = f_interpolate3d(u_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);
                    const real v_da_int = f_interpolate3d(v_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);
                    const real w_da_int = f_interpolate3d(w_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);
                    const real T_da_int = f_interpolate3d(T_da, ox,oy,oz, i_da,j_da,k_da, nx_da,ny_da,nz_da);

//                    if (T_da_int < 273.0 || T_da_int > 320.0) { std::cout << "\n----------\nT_da_int = " << T_da_int << "\n----------\n" << std::endl; }

                    bcTypes_f[id] = BCTypes::BCRegion;

//                    rho_obs[id] = (real)1.0;
                    rho_obs[id] = FuncThermalConvection::rho_buoyancy(T_da_int);
                    T_obs  [id] = T_da_int;

                    u_obs[id] = u_da_int / c_ref;
                    v_obs[id] = v_da_int / c_ref;
                    w_obs[id] = w_da_int / c_ref;


                    dirichlet_weight[id] = coef_dirichlet_weight;
                    viscosity_weight[id] = coef_viscosity_weight;
                }
            }


        }

    }

    delete [] u_da;
    delete [] v_da;
    delete [] w_da;
    delete [] T_da;
    delete [] E_da;
}


void  BoundaryConditions::ReadGroundData(
    const int           time_num,
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << std::endl << __PRETTY_FUNCTION__ << " : time = " << time_num << " minutes" << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
//    const real c_ref   = parameters.c_ref_lbm();

    const int  nx_da = 16;
    const int  ny_da = 16;

//    const real pitch_da[2] = { 500.0, 500.0 };
    const real pitch_da[3] = { 500.0, 500.0, 50.0 };

    float* T_da = new float[nx_da*ny_da];
    ReadOklahoma_GroundTb(T_da, nx_da, ny_da, time_num);

    auto fx_up_da = [nx_da, ny_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+1          ] - f[id_da] :  f[id_da] - f[id_da-1          ]; };
    auto fy_up_da = [nx_da, ny_da](const float* f, const float xup, int id_da){ return (xup >= 0.0f) ?  f[id_da+nx_da      ] - f[id_da] :  f[id_da] - f[id_da-nx_da      ]; };

#pragma omp parallel for schedule(dynamic)
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real   dx = meshValues[lv].coordinates().dx();
//        const real   dt = dx/c_ref;


        // check position from map offsets //
//        const real dx_map    =  4.0;
        const real x_min_map = -4000.0;
        const real y_min_map = -4000.0;
        const real z_min_map =     0.0;
        const real x_leaf[2] = { x[0] - x_min_map, x[0] + dx*DefAMR::NX_LEAF - x_min_map };
        const real y_leaf[2] = { y[0] - y_min_map, y[0] + dx*DefAMR::NX_LEAF - y_min_map };
        const real z_leaf[3] = { z[0] - z_min_map, z[0] + dx*DefAMR::NX_LEAF*(real)0.5 - z_min_map, z[0] + dx*DefAMR::NX_LEAF - z_min_map };

        // add //
//        if (z_leaf[0]/pitch_da[2] <= 0 || z_leaf[1]/pitch_da[2] <= 0 || z_leaf[2]/pitch_da[2] <= 0) { /* do nothing */; }
        if (z_leaf[0]/pitch_da[2] <= 1 || z_leaf[1]/pitch_da[2] <= 1 || z_leaf[2]/pitch_da[2] <= 1) { /* do nothing */; }
        else    { continue; }
        // add //

        std::vector<int>  ids_da;
        for (int jj=0; jj<2; jj++) {
        for (int ii=0; ii<2; ii++) {
            const int  i_da = x_leaf[ii]/pitch_da[0];
            const int  j_da = y_leaf[jj]/pitch_da[1];

            const bool is_cal_region = ( i_da>=0 && i_da <=nx_da-1 && j_da>=0 && j_da <=ny_da-1 ) ? true : false;

            if (is_cal_region) { ids_da.push_back( i_da + j_da*nx_da ); }
        }
        }
        std::sort(ids_da.begin(), ids_da.end());
        ids_da.erase(std::unique(ids_da.begin(), ids_da.end()), ids_da.end());


        for (const int id_da : ids_da) {
            const int   i_da = id_da%nx_da;
            const int   j_da = id_da/nx_da;

            // box //
            const real width_da[3] = { pitch_da[0] + (real)2.5*dx,
                                       pitch_da[1] + (real)2.5*dx,
                                       fabs(parameters.coefDomain().z_global_domain_max - parameters.coefDomain().z_global_domain_min) };

            const real center_da[3] = { (i_da + (real)0.5)*pitch_da[0] + x_min_map,
                                        (j_da + (real)0.5)*pitch_da[1] + y_min_map,
                                        width_da[2]*(real)0.5 };

            const real box_min[3] = { center_da[0] - (real)0.5*width_da[0],
                                      center_da[1] - (real)0.5*width_da[1],
                                      center_da[2] - (real)1.0*width_da[2] - (real)4.0*dx };

            const real box_max[3] = { center_da[0] + (real)0.5*width_da[0],
                                      center_da[1] + (real)0.5*width_da[1],
                                      center_da[2] };
//                                      center_da[2] + (real)1.0*width_da[2] + dx };


            const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];

            int*   bcTypes_f = &meshValues[lv].valueObjLS().bcTypes_f()[offset];

            real*  rho_obs   = &meshValues[lv].valueObjLS().rho_obj()[offset];
            real*  u_obs     = &meshValues[lv].valueObjLS().u_obj()  [offset];
            real*  v_obs     = &meshValues[lv].valueObjLS().v_obj()  [offset];
            real*  w_obs     = &meshValues[lv].valueObjLS().w_obj()  [offset];
            real*  T_obs     = &meshValues[lv].valueObjLS().T_obj()  [offset];

            real*  dirichlet_weight = &meshValues[lv].valueObjLS().dirichlet_weight()[offset];
            real*  viscosity_weight = &meshValues[lv].valueObjLS().viscosity_weight()[offset];

            for (int id=0; id<nn_cell; id++) {
                if ( !FuncObj::is_obj(lv_obs[id]) ) { continue; }

                // lv //
                const real x_tmp = x[id] + (real)0.5*dx;
                const real y_tmp = y[id] + (real)0.5*dx;
                const real z_tmp = z[id] + (real)0.5*dx;

                if (z_tmp > box_max[2]) { continue; }

                const real xyz[3] = { x_tmp, y_tmp, z_tmp };

                const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);

                // in the box area //
                if ( FuncObj::is_obj(lv_obs[id]) && FuncObj::is_obj(box_lv_obs) ) {
                    const real ox = (x_tmp - center_da[0])/pitch_da[0];
                    const real oy = (y_tmp - center_da[1])/pitch_da[1];

                    // interpolation //
                    const real T_da_int = f_interpolate2d(T_da, ox,oy, i_da,j_da, nx_da, ny_da);

                    bcTypes_f[id] = BCTypes::BCRegion;

                    rho_obs[id] = (real)1.0;
                    T_obs  [id] = T_da_int;

                    u_obs[id] = (real)0.0;
                    v_obs[id] = (real)0.0;
                    w_obs[id] = (real)0.0;

                    dirichlet_weight[id] = (real)1.0;
                    viscosity_weight[id] = (real)0.0;
                }
            }


        }

    }

    delete [] T_da;
}


void  BoundaryConditions::ReadOklahoma_uvwTE(float* u, float* v, float* w, float* T, float* E, int nx, int ny, int nz, int time_num)
{
    if (nx != 16 || ny != 16 || nz != 50) { std::cout << "error : number of grid points in wrf data\n"; exit(-1); }

    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(4) << time_num;

    const std::string fname = Foldernames::input_folder + "/oklahoma/boundary_uvwTE/wrf_data_16_16_50_" + sout.str() + ".ssv";
    std::ifstream  fin;
    fin.open(fname.c_str(), std::ios::in);
    if (!fin) { std::cout << __PRETTY_FUNCTION__ << " : error fin" << std::endl; exit(0); }

    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        float _u, _v, _w, _T, _E;

        fin >> _u >> _v >> _w >> _T >> _E;

//        if (rank_ == 0) { std::cout << _u << "," << _v << "," << _w << "," << _T << "," << _E << ","; }

        const int id = i + nx*j + nx*ny*k;
        u[id] = _u;
        v[id] = _v;
        w[id] = _w;
        T[id] = _T;
        E[id] = _E;

//        if (rank_ == 0) { std::cout << kk << ","; }
    }
    }
    }
}


void  BoundaryConditions::ReadOklahoma_GroundTb(float* Tb, int nx, int ny, int time_num)
{
    if (nx != 16 || ny != 16) { std::cout << "error : number of grid points in wrf data\n"; exit(-1); }

    std::ostringstream sout;
    sout << std::setfill('0') << std::setw(4) << time_num;

    const std::string fname = Foldernames::input_folder + "/oklahoma/boundary_Tb/wrf_data_Tb_16_16_" + sout.str() + ".ssv";
    std::ifstream  fin;
    fin.open(fname.c_str(), std::ios::in);
    if (!fin) { std::cout << __PRETTY_FUNCTION__ << " : error fin" << std::endl; exit(0); }

    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        float _T;

        fin >> _T;

//        if (rank_ == 0) { std::cout << _u << "," << _v << "," << _w << "," << _T << "," << _E << ","; }

        const int id = i + nx*j;
        Tb[id] = _T;
//        Tb[id] = _T * 10.0;
//        Tb[id] = i + nx*j;

//        if (rank_ == 0) { std::cout << kk << ","; }
    }
    }
}


real  BoundaryConditions::f_interpolate2d(const float* f, float xup, float yup, int i, int j, int nx, int ny)
{
    // xup,yup : 0 ~ 1 //
    auto index = [nx,ny](int i, int j){
        const int _i = (i+nx)%nx;
        const int _j = (j+ny)%ny;

        return  _i + nx*_j;
    };

    const int iup = (xup >= 0.0f) ? i+1 : i-1;
    const int jup = (yup >= 0.0f) ? j+1 : j-1;

    float  f2d[4] = {
        f[ index(i  ,j  ) ],
        f[ index(iup,j  ) ],
        f[ index(i  ,jup) ],
        f[ index(iup,jup) ],
    };

    const float xratio = std::fmax( std::fmin(fabs(xup), (float)1.0), (float)0.0 );
    const float yratio = std::fmax( std::fmin(fabs(yup), (float)1.0), (float)0.0 );

    const float ftmpx[2] = { f2d[0]*( (float)1.0 - xratio ) + f2d[1]*xratio,
                             f2d[2]*( (float)1.0 - xratio ) + f2d[3]*xratio };

    const float ftmpy = ftmpx[0]*( (float)1.0 - yratio ) + ftmpx[1]*yratio;

    return ftmpy;
}


real  BoundaryConditions::f_interpolate3d(const float* f, float xup, float yup, float zup, int i, int j, int k, int nx, int ny, int nz)
{
    // xup,yup : 0 ~ 1 //
    auto index = [nx,ny,nz](int i, int j, int k){
        const int _i = (i+nx)%nx;
        const int _j = (j+ny)%ny;
//        const int _k = (k+nz)%nz;
        const int _k = (k < 0   ) ? 0   :
                       (k > nz-1) ? nz-1:
                                    k;

        return  _i + nx*_j + nx*ny*_k;
    };

    const int iup = (xup >= 0.0f) ? i+1 : i-1;
    const int jup = (yup >= 0.0f) ? j+1 : j-1;
    const int kup = (zup >= 0.0f) ? k+1 : k-1;

    const float  f3d[8] = {
        f[ index(i  ,j  ,k  ) ],
        f[ index(iup,j  ,k  ) ],
        f[ index(i  ,jup,k  ) ],
        f[ index(iup,jup,k  ) ],
        f[ index(i  ,j  ,kup) ],
        f[ index(iup,j  ,kup) ],
        f[ index(i  ,jup,kup) ],
        f[ index(iup,jup,kup) ],
    };

    const float xratio = std::fmax( std::fmin(fabs(xup), (float)1.0), (float)0.0 );
    const float yratio = std::fmax( std::fmin(fabs(yup), (float)1.0), (float)0.0 );
    const float zratio = std::fmax( std::fmin(fabs(zup), (float)1.0), (float)0.0 );

    const float ftmpx[4] = { f3d[0]*( (float)1.0 - xratio ) + f3d[1]*xratio,
                             f3d[2]*( (float)1.0 - xratio ) + f3d[3]*xratio,
                             f3d[4]*( (float)1.0 - xratio ) + f3d[5]*xratio,
                             f3d[6]*( (float)1.0 - xratio ) + f3d[7]*xratio };

    const float ftmpy[2] = { ftmpx[0]*( (float)1.0 - yratio ) + ftmpx[1]*yratio,
                             ftmpx[2]*( (float)1.0 - yratio ) + ftmpx[3]*yratio };

    const float ftmpz = ftmpy[0]*( (float)1.0 - zratio ) + ftmpy[1]*zratio;

    return ftmpz;
}
