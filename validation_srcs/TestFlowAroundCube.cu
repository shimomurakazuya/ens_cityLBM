#include "TestFlowAroundCube.h"
#if 0

#include "defineAMR.h"
#include "Parser.h"
#include "WorkerThread.h"
#include "FuncLBM.h"
#include "FuncObj.h"
#include "FuncMapData.h"
#include "OutputFunc.h"

#include "TimeEvolutionLoop.h"
#include "LBMCalculation.h"

#include "IndexLBM.h"


void  TestFlowAroundCube::FlowAroundCube(int argc, char* argv[])
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

        InitBCAroundCube( field.grids(),
                          field.tree(),
                          field.parameters(),
                          field.meshValues0() );
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


void  TestFlowAroundCube::Init(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    InitObjectFlowAroundCube(grids, tree, parameters, meshValues, box_min_, box_max_);

    InitValueFlowAroundCube (grids, tree, parameters, meshValues);
    InitSource (grids, tree, parameters, meshValues);

    InitBCAroundCube(grids, tree, parameters, meshValues);
}


void  TestFlowAroundCube::InitValueFlowAroundCube(
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
//        const real*  x = &meshValues[lv].coordinates().x()[offset];
//        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
//        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        // NS //
        real*  u   = &meshValues[lv].valueNS().u()  [offset];
        real*  v   = &meshValues[lv].valueNS().v()  [offset];
        real*  w   = &meshValues[lv].valueNS().w()  [offset];
        real*  rho = &meshValues[lv].valueNS().rho()[offset];

        const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        const real  c_ref  = parameters.c_ref_lbm();
        for (int j=0; j<nn_cell; j++) {
            const real  rho_tmp = 1.0;
            const real  u_tmp   = FuncObj::is_fluid(lv_obs[j]) ? FuncFlowAroundCube::inflow_velocity(z[j]) : 0.0; // @ header file //
//            const real  u_tmp   = 0.0;
            const real  v_tmp   = 0.0;
            const real  w_tmp   = 0.0;

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


void  TestFlowAroundCube::InitObjectFlowAroundCube(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues,
    const real          box_min[],
    const real          box_max[]

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
        real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obs  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obs    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obs    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obs    = &meshValues[lv].valueObjLS().w_obj()  [offset];

        for (int j=0; j<nn_cell; j++) {
            const real  xyz[] = { x[j], y[j], z[j] };

            // lv //
            const real z_floar = 0.0;
            const real z_roof  = xyz_max_[2];

            const real floar_lv_obs = z_floar - z[j];
            const real roof_lv_obs  = z[j] - z_roof;

            const real box_lv_obs   =   FuncObj::length_from_box(xyz, box_min, box_max);

//            const real top_min[3] = { parameters.coefDomain().x_global_domain_min, parameters.coefDomain().y_global_domain_min, parameters.coefDomain().z_global_domain_max - (real)dxc*2.0 };
//            const real top_max[3] = { parameters.coefDomain().x_global_domain_max, parameters.coefDomain().y_global_domain_max, parameters.coefDomain().z_global_domain_max - (real)dxc };
//            const real top_lv_obs = FuncObj::length_from_box(xyz, top_min, top_max);

            const real bottom_min[3] = { parameters.coefDomain().x_global_domain_min, parameters.coefDomain().y_global_domain_min, parameters.coefDomain().z_global_domain_min + (real)dxc };
            const real bottom_max[3] = { parameters.coefDomain().x_global_domain_max, parameters.coefDomain().y_global_domain_max, z_floar };
            const real bottom_lv_obs = FuncObj::length_from_box(xyz, bottom_min, bottom_max);

            real tmp_lv_obs = box_lv_obs;
            tmp_lv_obs = FuncObj::overlap_objects(tmp_lv_obs, floar_lv_obs);
            tmp_lv_obs = FuncObj::overlap_objects(tmp_lv_obs, roof_lv_obs);
            tmp_lv_obs = FuncObj::overlap_objects(tmp_lv_obs, bottom_lv_obs);
//            tmp_lv_obs = FuncObj::overlap_objects(tmp_lv_obs, top_lv_obs);

            // fin //
            tmp_lv_obs = levelset_from_Inflow_fin(tmp_lv_obs, xyz);

            lv_obs[j] = tmp_lv_obs;

            // vel //
            rho_obs[j] = 1.0;
            u_obs  [j] = 0.0;
            v_obs  [j] = 0.0;
            w_obs  [j] = 0.0;

            if ( fabs(roof_lv_obs) <= 2.0*dx ) {
                u_obs[j] = FuncFlowAroundCube::inflow_velocity(z_roof) / c_ref;
                v_obs[j] = 0.0;
                w_obs[j] = 0.0;
            }
        }
    }
}


void  TestFlowAroundCube::InitBCAroundCube(
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

    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        int*  bcTypes_scaler = &meshValues[lv].valueObjLS().bcTypes_scaler()[offset];
        int*  bcTypes_vector = &meshValues[lv].valueObjLS().bcTypes_vector()[offset];

        real*     bc_rho = &meshValues[lv].valueObjLS().bc_rho()[offset];
        real*     bc_u   = &meshValues[lv].valueObjLS().bc_u  ()[offset];
        real*     bc_v   = &meshValues[lv].valueObjLS().bc_v  ()[offset];
        real*     bc_w   = &meshValues[lv].valueObjLS().bc_w  ()[offset];

        real*     bc_weight = &meshValues[lv].valueObjLS().bc_weight()[offset];

        const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        for (int j=0; j<nn_cell; j++) {
            const real  xyz[] = { x[j], y[j], z[j] };


            bc_rho[j] = 1.0;
            bc_u  [j] = FuncObj::is_fluid(lv_obs[j]) ? FuncFlowAroundCube::inflow_velocity(xyz[2]) / c_ref : 0.0;
            bc_v  [j] = 0.0;
            bc_w  [j] = 0.0;

            // roof bc //
            if ( xyz[2] > xyz_max_[2] - 2.0*dx ) {
                bc_u[j] = 0.0;
            }

            // upper bc //
//            if ( xyz[2] > xyz_max_[2] - bc_layer_max_[2] ) {
//                bcTypes_scaler[j] = BCNeumann;
//                bcTypes_vector[j] = BCNeumann;
//
//                bc_weight[j] = 0.05;
//            }


//            // side bc //
//            if (xyz[1] < xyz_min_[1] + bc_layer_min_[1] || xyz[1] > xyz_max_[1] - bc_layer_max_[1]) { 
//                bcTypes_scaler[j] = BCNeumann;
//                bcTypes_vector[j] = BCNeumann;
//
//                bc_weight[j] = 0.05;
//            }


            // inflow & outflow //
            // neumann //
            if (xyz[0] < xyz_min_[0] + bc_layer_min_[0] || xyz[0] > xyz_max_[0] - bc_layer_max_[0]) {
                bcTypes_scaler[j] = BCNeumann;
                bcTypes_vector[j] = BCNeumann;

                bc_weight[j] = 0.05;
            }


            // dirichlet //
            if (xyz[0] < xyz_min_[0] + bc_dirichlet_min_[0]) { // inflow //
                bcTypes_scaler[j] = BCDirichlet;
                bcTypes_vector[j] = BCDirichlet;

                bc_rho[j] = 1.0;
                bc_u  [j] = FuncObj::is_fluid(lv_obs[j]) ? FuncFlowAroundCube::inflow_velocity(xyz[2]) / c_ref : 0.0;

                bc_weight[j] = 1.0;
//                bc_weight[j] = fmin( 1.0 - fabs( xyz_min_[0] - xyz[0] )/bc_dirichlet_min_[0], 1.0 ); // diverge? //
            }
            if ( xyz[0] > xyz_max_[0] - bc_dirichlet_max_[0]) { // outflow //
                bcTypes_scaler[j] = BCDirichlet;
                bcTypes_vector[j] = BCDirichlet;

                bc_rho[j] = 1.0;
                bc_u  [j] = FuncObj::is_fluid(lv_obs[j]) ? FuncFlowAroundCube::inflow_velocity(xyz[2]) / c_ref : 0.0;

                bc_weight[j] = 1.0;
//                bc_weight[j] = fmin( 1.0 - fabs( xyz_max_[0] - xyz[0] )/bc_dirichlet_max_[0], 1.0 ); // diverge //
            }


        }
    }
}


void  TestFlowAroundCube::InitSource(
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

    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        real*  sc_scaler = &meshValues[lv].valueObjLS().sc_scalar()[offset];

        for (int j=0; j<nn_cell; j++) {
            const real  xyz[] = { x[j], y[j], z[j] };

            sc_scaler[j] = 0.0;

            if (    xyz[0] - 0.5*dx <= sc_xyz_[0] && sc_xyz_[0] < xyz[0] + 0.5*dx
                 && xyz[1] - 0.5*dx <= sc_xyz_[1] && sc_xyz_[1] < xyz[1] + 0.5*dx
                 && xyz[2] - 0.5*dx <= sc_xyz_[2] && sc_xyz_[2] < xyz[2] + 0.5*dx ) {
                sc_scaler[j] = sc_dscalar_ / (dx*dx*dx) * 1000.0; // ppm = mg/L : kg / m3 =  1000 * 1000 mg / 1000 L //

                std::cout << "**************************************************" << std::endl;
                std::cout << "sc_scaler[j] = " << sc_scaler[j] << std::endl;
                std::cout << "**************************************************" << std::endl;

            }

        }
    }
}


real  TestFlowAroundCube::levelset_from_Inflow_fin(const real lv_obs, const real xyz[]) const
{
    const int  num_fins = num_fins_;
    const real length = (xyz_max_[1] - xyz_min_[1]) / (num_fins + 1);

    real  tmp_lv_obs = lv_obs;
    for (int i=0; i<num_fins; i++) {
        const real  box_min[3] = { fin_pos_             , xyz_min_[1] + length*(i+1) - fin_width_*(real)0.5, xyz_min_[2]              };
        const real  box_max[3] = { fin_pos_ + fin_width_, xyz_min_[1] + length*(i+1) + fin_width_*(real)0.5, xyz_min_[2] + fin_height_ };

        const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);

        tmp_lv_obs = FuncObj::overlap_objects(tmp_lv_obs, box_lv_obs);
    }

    return  tmp_lv_obs;
}
#endif
