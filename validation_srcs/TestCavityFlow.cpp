#include "TestCavityFlow.h"

#include <random>
#include "defineAMR.h"
#include "Parser.h"
#include "WorkerThread.h"
#include "FuncLBM.h"
#include "FuncObj.h"
#include "OutputFunc.h"

#include "TimeEvolutionLoop.h"
#include "LBMCalculation.h"

#include "IndexLBM.h"


void  TestCavityFlow::Flow()
{

    Field  field(opt_, comm_);

    if (rank_ == 0) { std::cout << "preset_field\n"; }
    field.init_field_wo_map();

    // initialize //
    if (rank_ == 0) { std::cout << "Init\n"; }
    Init( field.grids(),
            field.tree(),
            field.parameters(),
            field.meshValues0() );

    if (rank_ == 0) { std::cout << "init_taskID_opt\n"; }
    field.init_taskID_opt(field.taskID(), field.tree(), field.meshValues0() );

    field.
    copy_meshValues(
            field.meshValues1(),
            field.meshValues0(),
            field.grids(),
            field.tree()
            );

    if (comm_.is_rank0()) { std::cout << "init_taskID_opt\n"; }
    field.init_taskID_opt(field.taskID(), field.tree(), field.meshValues0() );



#if 1
    /// submit functions ///
    std::vector< std::function<void(int)> >  cal_funcs;
    WorkerThread        iothread;

    // calculation //
//    LBMCalculation      lbmCalculation;
//    cal_funcs.push_back( [&field, &lbmCalculation](int t){ lbmCalculation.LBM_incompressible_flow(t, field); } );


    TimeEvolutionLoop  timeEvolutionLoop(
        comm_,
        field.parameters().step_init(),  field.parameters().step_end(),
        field.parameters().time_init(),  field.parameters().time_end(),
        field.parameters().dt0()
        );

    timeEvolutionLoop.time_evolution( field, iothread );
#endif
}


void  TestCavityFlow::Init(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    constexpr real  box_min[] = { -0.5, -0.5, -1.0 };
    constexpr real  box_max[] = {  0.5,  0.5,  1.0 };

    Init_Value_Cavity2D (grids, tree, parameters, meshValues);
    Init_Object_Cavity2D(grids, tree, parameters, meshValues, box_min, box_max);
}


void  TestCavityFlow::Init_Value_Cavity2D(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;


    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
//        const real*  x = &meshValues[lv].coordinates().x()[offset];
//        const real*  y = &meshValues[lv].coordinates().y()[offset];
//        const real*  z = &meshValues[lv].coordinates().z()[offset];
//        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        // NS //
        real*  u   = &meshValues[lv].valueNS().u()  [offset];
        real*  v   = &meshValues[lv].valueNS().v()  [offset];
        real*  w   = &meshValues[lv].valueNS().w()  [offset];
        real*  rho = &meshValues[lv].valueNS().rho()[offset];

        const real  rho_tmp = 1.0;
        const real  u_tmp   = 0.0;
        const real  v_tmp   = 0.0;
        const real  w_tmp   = 0.0;
        const real  c_ref  = parameters.c_ref_lbm();
        for (int j=0; j<nn_cell; j++) {
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


void  TestCavityFlow::Init_Object_Cavity2D(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues,
    const real          box_min[],
    const real          box_max[]
    )
{
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real  c_ref  = parameters.c_ref_lbm();

    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const auto& coordinates = meshValues[lv].coordinates();

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        real*  lv_obj = &meshValues[lv].valueObjLS().lv_obj()[offset];
        real*  u_obj  = &meshValues[lv].valueObjLS().u_obj() [offset];
        real*  v_obj  = &meshValues[lv].valueObjLS().v_obj() [offset];
        real*  w_obj  = &meshValues[lv].valueObjLS().w_obj() [offset];

        for (int j=0; j<nn_cell; j++) {
            const real xc = coordinates.xcell(offset + j);
            const real yc = coordinates.ycell(offset + j);
            const real zc = coordinates.zcell(offset + j);
            const real  xyz[] = { xc, yc, zc };

            // lv //
            const double lv_tmp  = FuncObj::length_from_box(xyz, box_min, box_max);
            const double nlv_obj = floor( lv_tmp/dx )*dx + 0.5*dx;

            lv_obj[j] = -nlv_obj;

            // vel //
            const int  xdim = 0;
            const int  ydim = 1;
            if (   fabs(lv_obj[j]) < 2.0*dx && FuncObj::length_from_plane(xyz, box_max[ydim], ydim) < 2.0*dx
                && xyz[xdim] > box_min[xdim] + dx
                && xyz[xdim] < box_max[xdim] - dx ) {
                u_obj[j] = 1.0/c_ref;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;
            }
            else {
                u_obj[j] = 0.0;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;
            }

        }
    }
}


