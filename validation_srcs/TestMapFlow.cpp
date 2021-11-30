#include "TestMapFlow.h"

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


void  TestMapFlow::Flow(int argc, char* argv[])
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

        InitBC( field.grids(),
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


void  TestMapFlow::Init(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    InitObjectMap(grids, tree, parameters, meshValues);
    InitObjectFloarRoof(grids, tree, parameters, meshValues);
//    ReconstructLevelset(grids, tree, parameters, meshValues);

    InitValue  (grids, tree, parameters, meshValues);
    InitSource (grids, tree, parameters, meshValues);

    InitBC(grids, tree, parameters, meshValues);
}


void  TestMapFlow::InitValue(
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

        const real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        const real  c_ref  = parameters.c_ref_lbm();
        for (int j=0; j<nn_cell; j++) {
            const real  rho_tmp = 1.0;
            const real  u_tmp   = FuncObj::is_fluid(lv_obj[j]) ? FuncMapFlow::inflow_velocity(z[j]) : 0.0; // @ header file //
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


void  TestMapFlow::InitObjectFloarRoof(
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
            // lv //
            const real z_floar = 0.0;
            const real z_roof = parameters.coefDomain().z_global_domain_max - dxc*(real)5.0;

            const real floar_lv_obj = z_floar - z[j];
            const real roof_lv_obj  = z[j] - z_roof;


            real tmp_lv_obj = lv_obj[j];
            tmp_lv_obj = FuncObj::overlap_objects(tmp_lv_obj, floar_lv_obj);
            tmp_lv_obj = FuncObj::overlap_objects(tmp_lv_obj, roof_lv_obj);

            lv_obj[j] = tmp_lv_obj;

            // vel //
            rho_obj[j] = 1.0;
            u_obj  [j] = 0.0;
            v_obj  [j] = 0.0;
            w_obj  [j] = 0.0;

            if ( fabs(roof_lv_obj) <= 2.0*dx ) {
                u_obj[j] = FuncMapFlow::inflow_velocity(z_roof) / c_ref;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;
            }
        }
    }
}


void  TestMapFlow::InitObjectMap(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int  n_leaf  = tree.number_of_nodes();

    const real dx_map = 1.0;
//    const real dx_map = 2.0;
    const real dxc    = 4.0;

    // map //
    const std::string fname = Foldernames::input_folder + "/" + Filenames::Map_name0;


    // read //
    std::ifstream  fin;
    fin.open(fname.c_str(), std::ios::in);
    if (!fin) { std::cout << __PRETTY_FUNCTION__ << " : error fin" << std::endl; exit(0); }

    int  mx, my;
    fin >> mx;
    fin >> my;

    if (rank_ == 0) {
        std::cout << "read map data" << std::endl;
        std::cout << mx << ", " << my << std::endl;
    }

    int* height_map;
    height_map = new int[mx*my];

    for (int _j=0; _j<my; _j++) {
    for (int i=0; i<mx; i++) {
//        const int j  = _j;
        const int j  = my-1 - _j;
        const int id = i + mx*j;

        fin >> height_map[id];
//        std::cout << height_map[id] << ", ";
    }
//    std::cout << std::endl;
    }
    fin.close();
    // map //

//    const real  cal_min_x = parameters.coefDomain().x_global_domain_min + dxc*(real)15.0;
//    const real  cal_max_x = parameters.coefDomain().x_global_domain_max - dxc*(real)15.0;
//
//    const real  cal_min_y = parameters.coefDomain().y_global_domain_min + dxc*(real)15.0;
//    const real  cal_max_y = parameters.coefDomain().y_global_domain_max - dxc*(real)15.0;

    const real  cal_min_x = parameters.coefDomain().x_global_domain_min + dxc*(real)4.0 + bc_layer_min_[0];
    const real  cal_max_x = parameters.coefDomain().x_global_domain_max - dxc*(real)4.0 - bc_layer_max_[0];

    const real  cal_min_y = parameters.coefDomain().y_global_domain_min + dxc*(real)4.0 + bc_layer_min_[1];
    const real  cal_max_y = parameters.coefDomain().y_global_domain_max - dxc*(real)4.0 - bc_layer_max_[1];

    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real   dx = meshValues[lv].coordinates().dx();

        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        // initialize //
        real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obj  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obj    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obj    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obj    = &meshValues[lv].valueObjLS().w_obj()  [offset];

        // map //
        const real  map_offset_x = 0.0;
        const real  map_offset_y = 0.0;

        FuncMapData::levelset_map(lv_obj, x,y,z, DefAMR::NX_LEAF, lv, dx, 
                                  height_map, mx, my, dx_map, map_offset_x, map_offset_y, 
                                  cal_min_x, cal_min_y, cal_max_x, cal_max_y);
    }
    delete [] height_map;
}


void  TestMapFlow::InitBC(
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

//    const real dxc = meshValues[0].coordinates().dx();
    const real dxc = 4.0;

    const real xyz_bc_min[3] = { parameters.coefDomain().x_global_domain_min + dxc*(real)2.0,
                                 parameters.coefDomain().y_global_domain_min + dxc*(real)2.0,
                                 parameters.coefDomain().z_global_domain_min + dxc*(real)2.0 };

    const real xyz_bc_max[3] = { parameters.coefDomain().x_global_domain_max - dxc*(real)2.0,
                                 parameters.coefDomain().y_global_domain_max - dxc*(real)2.0,
                                 parameters.coefDomain().z_global_domain_max - dxc*(real)2.0 };

    int  num_cal = 0;
    int  num_obj = 0;
    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        int*  bcTypes_f = &meshValues[lv].valueObjLS().bcTypes_f()[offset];

        const real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        real*  rho_obj  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obj    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obj    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obj    = &meshValues[lv].valueObjLS().w_obj()  [offset];

        real*  dirichlet_weight = &meshValues[lv].valueObjLS().dirichlet_weight()[offset];

        for (int j=0; j<nn_cell; j++) {
            const real  xyz[] = { x[j], y[j], z[j] };

            if ( FuncObj::is_obj(lv_obj[j]) ) { num_obj++; continue; }


            /// viscosity ///
            // upper bc //
            if ( xyz[2] > xyz_bc_max[2] - bc_layer_max_[2] ) {
                bcTypes_f[j] = BCTypes::BCRegion;

                u_obj[j] = FuncObj::is_fluid(lv_obj[j]) ? FuncMapFlow::inflow_velocity(xyz[2]) / c_ref : 0.0;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;

                dirichlet_weight[j] = 0.001;
            }
            // layer region //
            if (  (xyz[0] < xyz_bc_min[0] + bc_layer_min_[0] || xyz[0] > xyz_bc_max[0] - bc_layer_max_[0]) 
               || (xyz[1] < xyz_bc_min[1] + bc_layer_min_[1] || xyz[1] > xyz_bc_max[1] - bc_layer_max_[1]) ) {
                bcTypes_f[j] = BCTypes::BCRegion;

                u_obj[j] = FuncObj::is_fluid(lv_obj[j]) ? FuncMapFlow::inflow_velocity(xyz[2]) / c_ref : 0.0;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;

                dirichlet_weight[j] = 0.001;
            }


            /// dirichlet ///
            // upper bc //
            if ( xyz[2] > xyz_bc_max[2] - bc_dirichlet_max_[2] ) {
                bcTypes_f[j] = BCTypes::BCRegion;

                u_obj[j] = FuncObj::is_fluid(lv_obj[j]) ? FuncMapFlow::inflow_velocity(xyz[2]) / c_ref : 0.0;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;

                dirichlet_weight[j] = 1.0;
            }

            // inflow outflow //
            if (  (xyz[0] < xyz_bc_min[0] + bc_dirichlet_min_[0]) ||  (xyz[0] > xyz_bc_max[0] - bc_dirichlet_max_[0]) 
              ||  (xyz[1] < xyz_bc_min[1] - bc_dirichlet_min_[1]) ||  (xyz[1] > xyz_bc_max[1] - bc_dirichlet_max_[1]) ) {
                bcTypes_f[j] = BCTypes::BCRegion;

                u_obj[j] = FuncObj::is_fluid(lv_obj[j]) ? FuncMapFlow::inflow_velocity(xyz[2]) / c_ref : 0.0;
                v_obj[j] = 0.0;
                w_obj[j] = 0.0;

                dirichlet_weight[j] = 1.0;
            }


            num_cal++;
        }
    }

    std::cout << rank_ << " obj/cal_total (%) = " << (float)100.0*(float)num_obj/(float)(num_cal + num_obj) << std::endl;
}


void  TestMapFlow::InitSource(
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
                              5.0 };

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


void  TestMapFlow::ReconstructLevelset(
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
        const real   dx = meshValues[lv].coordinates().dx();

        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        // initialize //
        real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        for (int j=0; j<nn_cell; j++) {
            real  tmp_lv_obj = lv_obj[j];

            if ( tmp_lv_obj >= 0.0 && tmp_lv_obj <  0.5*dx ) { tmp_lv_obj =  0.5*dx; }
            if ( tmp_lv_obj <  0.0 && tmp_lv_obj > -0.5*dx ) { tmp_lv_obj = -0.5*dx; }

            lv_obj[j] = tmp_lv_obj;
        }
    }
}
