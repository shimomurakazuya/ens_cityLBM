#include "TestChannelFlow.h"

#include <random>
#include "defineAMR.h"
#include "defineFluidProperty.h"
#include "Parser.h"
#include "WorkerThread.h"
#include "FuncLBM.h"
#include "FuncObj.h"
#include "FuncMapData.h"
#include "OutputFunc.h"

#include "TimeEvolutionLoop.h"
#include "LBMCalculation.h"

#include "IndexLBM.h"


void  TestChannelFlow::Flow(int argc, char* argv[])
{
    Parser  parser(argc, argv);

    Field  field(parser.optionParser());
    MPI_Barrier(MPI_COMM_WORLD);

    if ( !field.optionParser()->restart_flags_and_step(0) ) {
//        field.init_field();
        if (rank_ == 0) { std::cout << "preset_field\n"; }
        field.preset_field();

        // initialize //
        if (rank_ == 0) { std::cout << "Init\n"; }
        Init( field.grids(),
              field.tree(),
              field.parameters(),
              field.meshValues0() );

        if (rank_ == 0) { std::cout << "init_taskID_opt\n"; }
        field.init_taskID_opt(field.taskID(), field.tree(), field.meshValues0() );
    }
    else {
        // restart //
        std::cout << __PRETTY_FUNCTION__ << " : " << "restart" << std::endl;

        field.read_field( field.optionParser()->restart_flags_and_step(1) );

        InitSource ( field.grids(),
                     field.tree(),
                     field.parameters(),
                     field.meshValues0() );

//        InitBC( field.grids(),
//                          field.tree(),
//                          field.parameters(),
//                          field.meshValues0() );
    }
    MPI_Barrier(MPI_COMM_WORLD);

    field.
    copy_meshValues(
        field.meshValues1(),
        field.meshValues0(),
        field.grids(),
        field.tree()
        );

    MPI_Barrier(MPI_COMM_WORLD);

    /// submit functions ///
//    if (rank_ == 0) { std::cout << "std::vector< std::function<void(int)> >  cal_funcs" << std::endl; }
    std::vector< std::function<void(int)> >  cal_funcs;
    WorkerThread        iothread;

    MPI_Barrier(MPI_COMM_WORLD);

    TimeEvolutionLoop  timeEvolutionLoop(
        field.parameters().step_init(),  field.parameters().step_end(),
        field.parameters().time_init(),  field.parameters().time_end(),
        field.parameters().dt0()
        );

    MPI_Barrier(MPI_COMM_WORLD);
    timeEvolutionLoop.time_evolution( field, iothread );
}


void  TestChannelFlow::Init(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    InitObjectFloarRoof(grids, tree, parameters, meshValues);

    InitValueChannelFlow (grids, tree, parameters, meshValues);
    InitSource (grids, tree, parameters, meshValues);

//    InitBC(grids, tree, parameters, meshValues);
}


void  TestChannelFlow::InitValueChannelFlow(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;


    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        // NS //
        real*  u   = &meshValues[lv].valueNS().u()  [offset];
        real*  v   = &meshValues[lv].valueNS().v()  [offset];
        real*  w   = &meshValues[lv].valueNS().w()  [offset];
        real*  rho = &meshValues[lv].valueNS().rho()[offset];

        const real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        const real  c_ref  = parameters.c_ref_lbm();


        auto vortex_xy = [](real& u_vortex, real& v_vortex, real x, real y, real xc, real yc, real length_lim) {
            const real  length     = sqrt( pow(x-xc, 2) + pow(y-yc,2) + 1.0e-5 );
            const real  vel_vortex = (length < length_lim) ? length/length_lim : 0.0;
        
            u_vortex   =  vel_vortex * y;
            v_vortex   = -vel_vortex * x;
        };

        auto vortex_yz = [](real& v_vortex, real& w_vortex, real y, real z, real yc, real zc, real length_lim) {
            const real  length     = sqrt( pow(y-yc, 2) + pow(z-zc,2) + 1.0e-5 );
            const real  vel_vortex = (length < length_lim) ? length/length_lim : 0.0;
        
            v_vortex   =  vel_vortex * z;
            w_vortex   = -vel_vortex * y;
        };

        auto vortex_zx = [](real& w_vortex, real& u_vortex, real z, real x, real zc, real xc, real length_lim) {
            const real  length     = sqrt( pow(z-zc, 2) + pow(x-xc,2) + 1.0e-5 );
            const real  vel_vortex = (length < length_lim) ? length/length_lim : 0.0;
        
            w_vortex   =  vel_vortex * x;
            u_vortex   = -vel_vortex * z;
        };


        for (int j=0; j<nn_cell; j++) {
            std::random_device rnd;
            std::mt19937 mt(rnd());
            constexpr int rand_half_max = 5;
            std::uniform_int_distribution<int> rand11(-rand_half_max, rand_half_max);

            const real  Re_tau = air_property::Re_tau;
            auto vel_channel = [Re_tau](const real zz) {
                const real _zz = fabs( (real)1.0 - fabs(zz) );

                return  (real)2.4 * log(_zz*Re_tau) + 5.5;
            };


            const real  amp_channel = 1.00;
            const real  c_vortex    = 0.6;
//            const real  c_vortex    = 1.5;
//            const real  c_vortex    = 2.0;
            const real  length_lim  = 0.5;
//            const real  length_lim  = 1.5;


            real u_vortex = 0.0, _u_vortex;
            real v_vortex = 0.0, _v_vortex;
            real w_vortex = 0.0, _w_vortex;
            for (int jj=0; jj<2; jj++) {
            for (int ii=0; ii<2; ii++) {
                const real xc = -3.0 + 2.0*ii;
                const real yc = -1.0 + 2.0*jj;
                const real zc = -0.5 + 1.0*jj;

                vortex_xy(_u_vortex, _v_vortex, x[j], y[j], xc,     yc,     length_lim);
                u_vortex += _u_vortex; v_vortex += _v_vortex;

                vortex_xy(_u_vortex, _v_vortex, x[j], y[j], xc+1.5, yc+0.5, length_lim);
                u_vortex -= _u_vortex; v_vortex -= _v_vortex;

                vortex_yz(_v_vortex, _w_vortex, y[j], z[j], yc,     zc,     length_lim);
                v_vortex += _v_vortex; w_vortex += _w_vortex;

                vortex_yz(_v_vortex, _w_vortex, y[j], z[j], yc+0.5, zc+0.1, length_lim);
                v_vortex -= _v_vortex; w_vortex -= _w_vortex;

                vortex_zx(_w_vortex, _u_vortex, z[j], x[j], zc,     xc,     length_lim);
                w_vortex += _w_vortex; u_vortex += _u_vortex;

                vortex_zx(_w_vortex, _u_vortex, z[j], x[j], zc+0.1, xc+1.1, length_lim);
                w_vortex -= _w_vortex; u_vortex -= _u_vortex;
            }
            }

            const real  rho_tmp = 1.0;
            const real  _u_tmp   = c_vortex*u_vortex + vel_channel(z[j]) * amp_channel;
            const real  _v_tmp   = c_vortex*v_vortex;
            const real  _w_tmp   = c_vortex*w_vortex;
//            const real  u_tmp   = 0.5*rand11(mt)/(2*rand_half_max) + u_vortex + vel_channel(z[j]) * amp_channel;
//            const real  v_tmp   = 0.5*rand11(mt)/(2*rand_half_max) + v_vortex;
//            const real  w_tmp   = 0.5*rand11(mt)/(2*rand_half_max);

            const real  u_tmp   = _u_tmp;
            const real  v_tmp   = _v_tmp;
            const real  w_tmp   = _w_tmp;
//            const real  u_tmp   = FuncObj::is_fluid(lv_obj[j] + 4.0*dx) ? _u_tmp : 0.0;
//            const real  v_tmp   = FuncObj::is_fluid(lv_obj[j] + 4.0*dx) ? _v_tmp : 0.0;
//            const real  w_tmp   = FuncObj::is_fluid(lv_obj[j] + 4.0*dx) ? _w_tmp : 0.0;

            rho[j] = rho_tmp;
            u  [j] = u_tmp / c_ref;
            v  [j] = v_tmp / c_ref;
            w  [j] = w_tmp / c_ref;
        }

        // LBM //
        const int  nQ  = LBM_velocity_model::nQ;
        const int  offset_lbm = offset * nQ;

        real*  f_lbm = &meshValues[lv].valueLBM().f_lbm()[offset_lbm];
        for (int ii=0; ii<nQ; ii++) {
            int  iv, jv, kv;
            IndexLBM::lbm_index(iv,jv,kv, ii);

            for (int j=0; j<nn_cell; j++) {
                f_lbm[j + nn_cell*ii] = FuncLBM::feq_D3Q27(rho[j], u[j],v[j],w[j], iv,jv,kv);
            }
        }
    }
}


void  TestChannelFlow::InitObjectFloarRoof(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues

    )
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real  c_ref  = parameters.c_ref_lbm();

    const real dxc = meshValues[0].coordinates().dx();

    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obj  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obj    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obj    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obj    = &meshValues[lv].valueObjLS().w_obj()  [offset];

        for (int j=0; j<nn_cell; j++) {
            const double zcell = z[j] + 0.5*dx;

            const bool  is_bottom_side = (zcell <= 0.0) ? true : false;

            const double z_floar = -1.0;
            const double z_roof  =  1.0;

            const double tmp_nlv_obj = (is_bottom_side) ? zcell - z_floar:
                                                          z_roof - zcell;

//            const int    ii = floor(tmp_nlv_obj/dx);
//            const double nlv_obj = (ii + 0.5)*dx;
            const double nlv_obj = floor( tmp_nlv_obj/dx )*dx + 0.5*dx;

            lv_obj[j] = - nlv_obj;

            // vel //
            rho_obj[j] = 1.0;
            u_obj  [j] = 0.0;
            v_obj  [j] = 0.0;
            w_obj  [j] = 0.0;
        }
    }
}


void  TestChannelFlow::InitSource(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;

    const real sc_dx     = 6.0*meshValues[DefAMR::LV_MAX-1].coordinates().dx();
    const real sc_xyz[3] =  { ( parameters.coefDomain().x_global_domain_min + parameters.coefDomain().x_global_domain_max )*(real)0.5,
                              ( parameters.coefDomain().y_global_domain_min + parameters.coefDomain().y_global_domain_max )*(real)0.5,
                              ( parameters.coefDomain().z_global_domain_min + parameters.coefDomain().z_global_domain_max )*(real)0.5 };

#pragma omp parallel for
    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        const real   dx = meshValues[lv].coordinates().dx();

        const real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        // initialize //
        real*  sc_scaler = &meshValues[lv].valueObjLS().sc_scalar()[offset];

        for (int j=0; j<nn_cell; j++) {
            const real  xyz[] = { x[j], y[j], z[j] };

            sc_scaler[j] = 0.0;

            if ( !FuncObj::is_fluid(lv_obj[j]) ) { continue; }

            if (    xyz[0] - 0.5*sc_dx <= sc_xyz[0] && sc_xyz[0] < xyz[0] + 0.5*sc_dx
                 && xyz[1] - 0.5*sc_dx <= sc_xyz[1] && sc_xyz[1] < xyz[1] + 0.5*sc_dx
                 && xyz[2] - 0.5*sc_dx <= sc_xyz[2] && sc_xyz[2] < xyz[2] + 0.5*sc_dx ) {
                sc_scaler[j] = sc_dscalar_ / (dx*dx*dx) * 1000.0; // ppm = mg/L : kg / m3 =  1000 * 1000 mg / 1000 L //

//                std::cout << "**************************************************" << std::endl;
//                std::cout << "sc_scaler[j] = " << sc_scaler[j] << std::endl;
//                std::cout << "**************************************************" << std::endl;

            }

        }
    }
}
