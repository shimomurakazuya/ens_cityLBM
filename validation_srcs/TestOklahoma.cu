#include "TestOklahoma.h"
#include <sstream>
#include <iomanip>

#include "defineAMR.h"
#include "defineFluidProperty.h"
#include "Parser.h"
#include "WorkerThread.h"
#include "FuncLBM.h"
#include "FuncObj.h"
#include "OutputFunc.h"
#include "FuncMapData.h"

#include "TimeEvolutionLoop.h"
#include "LBMCalculation.h"

#include "IndexLBM.h"
#include "BoundaryConditions.h"


void  TestOklahoma::Flow(int argc, char* argv[])
{
    Parser  parser(argc, argv);

    Field  field(parser.optionParser());

    if ( !field.optionParser()->restart_flags_and_step(0) ) {
//        field.init_field();
        field.preset_field();

        // initialize //
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
    }

    field.
    copy_meshValues(
        field.meshValues1(),
        field.meshValues0(),
        field.grids(),
        field.tree()
        );

#if 1
    /// submit functions ///
//    std::vector< std::function<void(int)> >  cal_funcs;
    WorkerThread        iothread;

    // calculation //
//    LBMCalculation      lbmCalculation;
//    cal_funcs.push_back( [&field, &lbmCalculation](int t){ lbmCalculation.LBM_incompressible_flow(t, field); } );


    TimeEvolutionLoop  timeEvolutionLoop(
        field.parameters().step_init(),  field.parameters().step_end(),
        field.parameters().time_init(),  field.parameters().time_end(),
        field.parameters().dt0()
        );

    timeEvolutionLoop.time_evolution( field, iothread );
#endif
}


void  TestOklahoma::Init(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    // object //
//    InitObjectMapTokyo(grids, tree, parameters, meshValues);
    InitObjectMapOklahoma(grids, tree, parameters, meshValues);

    InitObjectFloarRoof(grids, tree, parameters, meshValues);
    InitDecoBoco       (grids, tree, parameters, meshValues);

    DigitalizeObject   (grids, tree, parameters, meshValues);
    // object //

    InitValue (grids, tree, parameters, meshValues);
    InitValueOklahoma (grids, tree, parameters, meshValues);


    InitObservationDataOklahoma (grids, tree, parameters, meshValues);
//    InitializeByWRFData (grids, tree, parameters, meshValues);

    InitSource (grids, tree, parameters, meshValues);
}


void  TestOklahoma::InitValue(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;


    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
//        const real*  z = &meshValues[lv].coordinates().z()[offset];
//        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        // NS //
        real*  u   = &meshValues[lv].valueNS().u()  [offset];
        real*  v   = &meshValues[lv].valueNS().v()  [offset];
        real*  w   = &meshValues[lv].valueNS().w()  [offset];
        real*  rho = &meshValues[lv].valueNS().rho()[offset];

        real*  T   = &meshValues[lv].valueNS().T()  [offset];

        const real  rho_tmp = 1.0;
        const real  c_ref  = parameters.c_ref_lbm();
        for (int id=0; id<nn_cell; id++) {
            const real  u_tmp   =  0.0;
            const real  v_tmp   =  0.0;
            const real  w_tmp   =  0.0;

//            const real  T_tmp   = fluid_property::Temperature0;
//            const real  T_tmp   = fluid_property::Temperature0 + (TemperatureBotoomWall - fluid_property::Temperature0)*(real)0.15;
            real T_tmp = TemperatureGround_;

            rho[id] = rho_tmp;
            u  [id] = u_tmp / c_ref;
            v  [id] = v_tmp / c_ref;
            w  [id] = w_tmp / c_ref;

            T [id] = T_tmp;
        }

        // LBM //
        const int  nQ  = LBM_velocity_model::nQ;
        const int  offset_lbm = offset * nQ;

        real*  f_lbm = &meshValues[lv].valueLBM().f_lbm()[offset_lbm];
        for (int ii=0; ii<nQ; ii++) {
            int  iv, jv, kv;
            IndexLBM::lbm_index(iv,jv,kv, ii);

            for (int id=0; id<nn_cell; id++) {
                f_lbm[id + nn_cell*ii] = FuncLBM::feq_D3Q27(rho[id], u[id],v[id],w[id], iv,jv,kv);
            }
        }
    }
}


void  TestOklahoma::InitValueOklahoma(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    BoundaryConditions  boundaryConditions;
    boundaryConditions.InitializeByWRFData   (start_minutes_, grids, tree, parameters, meshValues);
    boundaryConditions.InitializeByGroundData(start_minutes_, grids, tree, parameters, meshValues);
}


void  TestOklahoma::InitObjectFloarRoof(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

    const real dxc = meshValues[0].coordinates().dx();

#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

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

#pragma omp simd
        for (int id=0; id<nn_cell; id++) {
            // lv //
            const real z_floar = Zmin_;
            const real z_roof  = parameters.coefDomain().z_global_domain_max - dxc*(real)3.0;
//            const real z_roof  = 1000.0;

            const real z_tmp = z[id] + 0.5*dx;
            const real floar_lv_obs = z_floar - z_tmp;
            const real roof_lv_obs  = z_tmp - z_roof;


            real tmp_lv_obs = lv_obs[id];
            tmp_lv_obs      = FuncObj::overlap_objects(tmp_lv_obs, floar_lv_obs);
            tmp_lv_obs      = FuncObj::overlap_objects(tmp_lv_obs, roof_lv_obs);

            lv_obs[id] = tmp_lv_obs;

            // vel //
            rho_obs[id] = 1.0;
            u_obs  [id] = 0.0;
            v_obs  [id] = 0.0;
            w_obs  [id] = 0.0;
        }
    }
}


void  TestOklahoma::DigitalizeObject(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

//    const real dxc = meshValues[0].               coordinates().dx();
//    const real dxf = meshValues[DefAMR::LV_MAX-2].coordinates().dx();

#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

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

        for (int id=0; id<nn_cell; id++) {
            // lv //
            double tmp_lv_obs = std::floor( (lv_obs[id] + dx)/dx )*dx - 0.5*dx;

            lv_obs[id] = tmp_lv_obs;
        }
    }
}


void  TestOklahoma::InitObjectMapTokyo(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int  n_leaf  = tree.number_of_nodes();

    const real dx_map = 1.0;
//    const real dx_map = 2.0;
//    const real dx_map = 4.0;
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

    const real  cal_min_x = parameters.coefDomain().x_global_domain_min + dxc*(real)1.0;
    const real  cal_max_x = parameters.coefDomain().x_global_domain_max - dxc*(real)1.0;

    const real  cal_min_y = parameters.coefDomain().y_global_domain_min + dxc*(real)1.0;
    const real  cal_max_y = parameters.coefDomain().y_global_domain_max - dxc*(real)1.0;


#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real   dx = meshValues[lv].coordinates().dx();

        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        // initialize //
        real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obs  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obs    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obs    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obs    = &meshValues[lv].valueObjLS().w_obj()  [offset];

        // map //
        const real  map_offset_x = -1500.0;
        const real  map_offset_y = -1500.0;

        FuncMapData::levelset_map(lv_obs, x,y,z, DefAMR::NX_LEAF, lv, dx,
                                  height_map, mx, my, dx_map, map_offset_x, map_offset_y,
                                  cal_min_x, cal_min_y, cal_max_x, cal_max_y);
    }
    delete [] height_map;
}


void  TestOklahoma::InitObjectMapOklahoma(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int  n_leaf  = tree.number_of_nodes();

    // map //
    // 1m //
    const real dx_map = 1.0;
    const std::string fname = Foldernames::input_folder + "/oklahoma/building/buildings1m_1200_1500.txt";

//    // 2m //
//    const real dx_map = 2.0;
//    const std::string fname = Foldernames::input_folder + "/oklahoma/building/buildings2m_600_750.txt";

//    // 4m //
//    const real dx_map = 4.0;
//    const std::string fname = Foldernames::input_folder + "/oklahoma/building/buildings4m_300_375.txt";
//    const std::string fname = Foldernames::input_folder + "/oklahoma/building/test_10_10.txt";


    // read //
    std::ifstream  fin;
    fin.open(fname.c_str(), std::ios::in);
    if (!fin) { std::cout << __PRETTY_FUNCTION__ << " : error fin" << std::endl; exit(0); }

    int  mx, my;
    fin >> mx;
    fin >> my;

    if (rank_ == 0) { std::cout << "read map data = " << mx << ", " << my << std::endl; }

    // map //
    const real  domain_center_x = (parameters.coefDomain().x_global_domain_min + parameters.coefDomain().x_global_domain_max)*0.5;
    const real  domain_center_y = (parameters.coefDomain().y_global_domain_min + parameters.coefDomain().y_global_domain_max)*0.5;

//    const real  map_offset_x = 0.0;
//    const real  map_offset_y = 0.0;
//    const real  map_offset_x = domain_center_x - mx*dx_map/2 + dx_map/2;
//    const real  map_offset_y = domain_center_y - my*dx_map/2 + dx_map/2;
    const real  map_offset_x = domain_center_x - mx*dx_map/2 + dx_map/2;
    const real  map_offset_y = domain_center_y - my*dx_map/2 + dx_map/2;


    int* height_map;
    height_map = new int[mx*my];

    for (int _j=0; _j<my; _j++) {
    for (int i=0; i<mx; i++) {
        const int j  = _j;
//        const int j  = my-1 - _j;
        const int id = i + mx*j;

        int _height;

        fin >> _height;
        height_map[id] = _height;
////        std::cout << height_map[id] << ", ";
    }
//    std::cout << std::endl;
    }
    fin.close();
    // map //

    const real  cal_min_x = parameters.coefDomain().x_global_domain_min + (real)16.0 + bc_layer_min_[0];
    const real  cal_max_x = parameters.coefDomain().x_global_domain_max - (real)16.0 - bc_layer_max_[0];

    const real  cal_min_y = parameters.coefDomain().y_global_domain_min + (real)16.0 + bc_layer_min_[1];
    const real  cal_max_y = parameters.coefDomain().y_global_domain_max - (real)16.0 - bc_layer_max_[1];


#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const real   dx = meshValues[lv].coordinates().dx();

        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];

        // initialize //
        real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obs  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obs    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obs    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obs    = &meshValues[lv].valueObjLS().w_obj()  [offset];

        FuncMapData::levelset_map(lv_obs, x,y,z, DefAMR::NX_LEAF, lv, dx,
                                  height_map, mx, my, dx_map, map_offset_x, map_offset_y,
                                  cal_min_x, cal_min_y, cal_max_x, cal_max_y);
    }
    delete [] height_map;
}


void  TestOklahoma::InitDecoBoco(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

    const real dxc = meshValues[0].coordinates().dx();

//    // 2m //
//    const real width_decoboco[3] = { 8.0, 8.0, 24.0 };
//    const real pitch_decoboco[3] = { width_decoboco[0]*8, width_decoboco[1]*8, width_decoboco[2]*4 };

    // 4m //
//    const real width_decoboco[3] = { 16.0, 16.0, 8.0 };
//    const real pitch_decoboco[3] = { width_decoboco[0]*8, width_decoboco[1]*8, width_decoboco[2]*4 };

    const real width_decoboco[3] = { 8.0, 8.0, 8.0 };
    const real pitch_decoboco[3] = { width_decoboco[0]*4, width_decoboco[1]*4, width_decoboco[2]*4 };
//    const real pitch_decoboco[3] = { width_decoboco[0]*8, width_decoboco[1]*8, width_decoboco[2]*4 };

//    const real width_decoboco[3] = { 20.0, 20.0, 20.0 };
//    const real pitch_decoboco[3] = { width_decoboco[0]*6, width_decoboco[1]*6, width_decoboco[2]*6 };

    const int nx_deco = fabs(parameters.coefDomain().x_global_domain_max - parameters.coefDomain().x_global_domain_min)/pitch_decoboco[0];
    const int ny_deco = fabs(parameters.coefDomain().y_global_domain_max - parameters.coefDomain().y_global_domain_min)/pitch_decoboco[1];


//#pragma omp parallel for collapse(2)
#pragma omp parallel for collapse(2) schedule(dynamic)
    for (int jj=0; jj<ny_deco; jj++) {
    for (int ii=0; ii<nx_deco; ii++) {
        for (int l=0; l<n_leaf; l++) {
            // deco boco //
            const real _ii = ii;
            const real _jj = (ii%2==0) ? jj : jj+0.5;
//            const real _jj =  (ii%4==0) ? jj        :
//                              (ii%4==1) ? jj + (real)0.25 :
//                              (ii%4==2) ? jj + (real)0.50 :
//                                          jj + (real)0.75;

            const real box_min[3] = { parameters.coefDomain().x_global_domain_min + _ii*pitch_decoboco[0],
                                      parameters.coefDomain().y_global_domain_min + _jj*pitch_decoboco[1],
                                      parameters.coefDomain().z_global_domain_min - (real)1.0   };

            const real box_max[3] = { box_min[0] + width_decoboco[0],
                                      box_min[1] + width_decoboco[1],
                                                   width_decoboco[2]  };
            // deco boco //

            const int  lv     = tree.nodes(l)->level();
            const int  offset = tree.nodes(l)->mesh_offset();

            // coordinates xyz //
            const real*  x = &meshValues[lv].coordinates().x()[offset];
            const real*  y = &meshValues[lv].coordinates().y()[offset];
            const real*  z = &meshValues[lv].coordinates().z()[offset];

            const real   dx = meshValues[lv].coordinates().dx();

            const real _x_tmp = x[0];
            const real _y_tmp = y[0];
            const real _z_tmp = z[0];
//            if ( _x_tmp < box_min[0] - width_decoboco[0]*1.5 || _x_tmp + dx*DefAMR::NX_LEAF > box_max[0] + width_decoboco[0]*1.5 ) { continue; }
//            if ( _y_tmp < box_min[1] - width_decoboco[1]*1.5 || _y_tmp + dx*DefAMR::NX_LEAF > box_max[1] + width_decoboco[1]*1.5 ) { continue; }
//            if ( _z_tmp < box_min[2] - width_decoboco[2]*1.5 || _z_tmp + dx*DefAMR::NX_LEAF > box_max[2] + width_decoboco[2]*1.5 ) { continue; }

            if ( _z_tmp < box_min[2] - dx*(real)2.5 || _z_tmp + dx*DefAMR::NX_LEAF > box_max[2] + dx*(real)2.5 ) { continue; }
            if ( _x_tmp < box_min[0] - dx*(real)2.5 || _x_tmp + dx*DefAMR::NX_LEAF > box_max[0] + dx*(real)2.5 ) { continue; }
            if ( _y_tmp < box_min[1] - dx*(real)2.5 || _y_tmp + dx*DefAMR::NX_LEAF > box_max[1] + dx*(real)2.5 ) { continue; }

            // initialize //
            real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];
            real*  rho_obs  = &meshValues[lv].valueObjLS().rho_obj()[offset];
            real*  u_obs    = &meshValues[lv].valueObjLS().u_obj()  [offset];
            real*  v_obs    = &meshValues[lv].valueObjLS().v_obj()  [offset];
            real*  w_obs    = &meshValues[lv].valueObjLS().w_obj()  [offset];

            for (int id=0; id<nn_cell; id++) {
                // lv //
                const real x_tmp = x[id] + (real)0.5*dx;
                const real y_tmp = y[id] + (real)0.5*dx;
                const real z_tmp = z[id] + (real)0.5*dx;

                const real xyz[3] = { x_tmp, y_tmp, z_tmp };

                const real  cal_min_x1 = parameters.coefDomain().x_global_domain_min + dxc*(real)4.0 + bc_layer_min_[0];
                const real  cal_max_x1 = parameters.coefDomain().x_global_domain_max - dxc*(real)4.0 - bc_layer_max_[0];

                const real  cal_min_y1 = parameters.coefDomain().y_global_domain_min + dxc*(real)4.0 + bc_layer_min_[1];
                const real  cal_max_y1 = parameters.coefDomain().y_global_domain_max - dxc*(real)4.0 - bc_layer_max_[1];

                const real  cal_min_x2 = parameters.coefDomain().x_global_domain_min + dxc*(real)4.0 + bc_layer_min_[0] + 128.0;
                const real  cal_max_x2 = parameters.coefDomain().x_global_domain_max - dxc*(real)4.0 - bc_layer_max_[0] - 128.0;

                const real  cal_min_y2 = parameters.coefDomain().y_global_domain_min + dxc*(real)4.0 + bc_layer_min_[1] + 128.0;
                const real  cal_max_y2 = parameters.coefDomain().y_global_domain_max - dxc*(real)4.0 - bc_layer_max_[1] - 128.0;

//                real tmp_lv_obs = lv_obs[id];
//                const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);
//                lv_obs[id] = FuncObj::overlap_objects(tmp_lv_obs, box_lv_obs);

                if ( x_tmp<cal_min_x1 || x_tmp>cal_max_x1 || y_tmp<cal_min_y1 || y_tmp>cal_max_y1 ) {
                    continue;
                }
                else {
                    real  bbox_min[3] = { box_min[0], box_min[1], box_min[2] };
                    real  bbox_max[3] = { box_max[0], box_max[1], box_max[2] };

                    if (fabs(x_tmp) < (real)330.0 && fabs(y_tmp - 60.0) < (real)620.0) {
                        if ( (fabs(x_tmp + (real)220.0) < (real)110.0 && fabs(y_tmp + (real)260.0) < (real)300.0)
                          || (fabs(x_tmp - (real)110.0) < (real)220.0 && fabs(y_tmp + (real)380.0) < (real)180.0)
                          || (fabs(x_tmp - (real)  0.0) < (real)330.0 && fabs(y_tmp - (real)580.0) < (real)100.0)
                          || (fabs(x_tmp - (real)240.0) < (real) 90.0 && fabs(y_tmp - (real)400.0) < (real) 80.0) ) {
                            bbox_max[2] = (real)4.0;
                        }
                        else {
                            continue;
                        }
                    }

                    const real box_lv_obs = FuncObj::length_from_box(xyz, bbox_min, bbox_max);

                    real tmp_lv_obs = lv_obs[id];
                    lv_obs[id] = FuncObj::overlap_objects(tmp_lv_obs, box_lv_obs);
                }

            }
        }

    }
    }

}


void  TestOklahoma::InitSource(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const real source_xyz[] = { -48.0, -164.0, 2.0 };
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;

//    const float source_position = { 10.0, -170.0,

    int  num_source = 0;
    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x = &meshValues[lv].coordinates().x()[offset];
        const real*  y = &meshValues[lv].coordinates().y()[offset];
        const real*  z = &meshValues[lv].coordinates().z()[offset];
        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        const real*  lv_obs   = &meshValues[lv].valueObjLS().lv_obj() [offset];
              real*  sc_scalar   = &meshValues[lv].valueObjLS().sc_scalar()  [offset];

        const real  c_ref  = parameters.c_ref_lbm();
        for (int id=0; id<nn_cell; id++) {
            const real x_tmp = x[id] + (real)0.5*dx;
            const real y_tmp = y[id] + (real)0.5*dx;
            const real z_tmp = z[id] + (real)0.5*dx;

//            const real box_min[] = { -(real)50.0-(real)0.5*dx, -(real)150.0-(real)0.5*dx, (real)0.0 };
//            const real box_max[] = { -(real)50.0+(real)0.5*dx, -(real)150.0+(real)0.5*dx,       dx  };
            const real box_min[] = { -source_xyz[0]-(real)1.0*dx, -source_xyz[1]-(real)1.0*dx, source_xyz[2]-(real)1.0*dx };
            const real box_max[] = { -source_xyz[0]+(real)1.0*dx, -source_xyz[1]+(real)1.0*dx, source_xyz[2]+(real)1.0*dx };

//            // tokyo //
//            const real box_min[] = { -(real)2.0*dx, -(real)2.0*dx, (real)10.0 + (real)8.0+(real)0.0*dx };
//            const real box_max[] = { +(real)2.0*dx, +(real)2.0*dx, (real)10.0 + (real)8.0+(real)4.0*dx  };

            const real xyz[] = { x_tmp, y_tmp, z_tmp };
            const real box_lv_obs = FuncObj::length_from_box(xyz, box_min, box_max);

            if ( FuncObj::is_fluid(lv_obs[id]) && FuncObj::is_obj(box_lv_obs) ) {
                sc_scalar[id] = (real)1.0;

                if ( tree.nodes(i)->nodeCalFlags().Cal() ) {
                    num_source++;
                }
            }

        }

    }

    MPI_Request requ[1];
    int  num_source_g = 0;
    MPI_Iallreduce(&num_source, &num_source_g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &requ[0]);
    MPI_Waitall(1, requ, MPI_STATUSES_IGNORE);

    std::cout << rank_ << " : num_source, num_source_global = " << num_source << ", " << num_source_g << std::endl;
    if ( num_source_g == 0 ) { std::cout << "error : " << __PRETTY_FUNCTION__ << std::endl; exit(-1); }


    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

              real*  sc_scalar   = &meshValues[lv].valueObjLS().sc_scalar()  [offset];

        for (int j=0; j<nn_cell; j++) {
            sc_scalar[j] = sc_scalar[j]/num_source_g;
        }

    }
}


void  TestOklahoma::InitObservationDataOklahoma(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int time_minutes_now = start_minutes_;
    BoundaryConditions  boundaryConditions;
    boundaryConditions.ReadWRFData   (time_minutes_now+1, grids, tree, parameters, meshValues);
    boundaryConditions.ReadGroundData(time_minutes_now+1, grids, tree, parameters, meshValues);

    // swap: rho_obj <-> rhon_obj, .. etc //
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        meshValues[lv].valueObjLS().swap_objs();
    }
    boundaryConditions.ReadWRFData   (time_minutes_now+1, grids, tree, parameters, meshValues);
    boundaryConditions.ReadGroundData(time_minutes_now+1, grids, tree, parameters, meshValues);
}


void  TestOklahoma::ReadOklahoma_uvwTE(float* u, float* v, float* w, float* T, float* E, int nx, int ny, int nz, int time_num)
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


void  TestOklahoma::Refinement_uvwTE(float* uf, float* vf, float* wf, float* Tf, float* Ef,
    const float* u, const float* v, const float* w, const float* T, const float* E,
    int nx, int ny, int nz, int cx, int cy, int cz)
{
    auto fx = [nx,ny,nz](const float* f, int id){ return (f[id+1]     - f[id-1]    )*(float)0.5; };
    auto fy = [nx,ny,nz](const float* f, int id){ return (f[id+nx]    - f[id-nx]   )*(float)0.5; };
    auto fz = [nx,ny,nz](const float* f, int id){ return (f[id+nx*ny] - f[id-nx*ny])*(float)0.5; };

    const int nx_da = nx*cx;
    const int ny_da = ny*cy;
//    const int nz_da = nz*cz;


#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = i + nx*j + nx*ny*k;

        const float _u = u[id];
        const float _v = v[id];
        const float _w = w[id];
        const float _T = T[id];
        const float _E = E[id];

        const float _ux = (i == 0 || i == nx-1) ? 0.0 : fx(u, id);    const float _uy = (j == 0 || j == ny-1) ? 0.0 : fy(u, id);    const float _uz = (k == 0 || k == nz-1) ? 0.0 : fz(u, id);
        const float _vx = (i == 0 || i == nx-1) ? 0.0 : fx(v, id);    const float _vy = (j == 0 || j == ny-1) ? 0.0 : fy(v, id);    const float _vz = (k == 0 || k == nz-1) ? 0.0 : fz(v, id);
        const float _wx = (i == 0 || i == nx-1) ? 0.0 : fx(w, id);    const float _wy = (j == 0 || j == ny-1) ? 0.0 : fy(w, id);    const float _wz = (k == 0 || k == nz-1) ? 0.0 : fz(w, id);
        const float _Tx = (i == 0 || i == nx-1) ? 0.0 : fx(T, id);    const float _Ty = (j == 0 || j == ny-1) ? 0.0 : fy(T, id);    const float _Tz = (k == 0 || k == nz-1) ? 0.0 : fz(T, id);
        const float _Ex = (i == 0 || i == nx-1) ? 0.0 : fx(E, id);    const float _Ey = (j == 0 || j == ny-1) ? 0.0 : fy(E, id);    const float _Ez = (k == 0 || k == nz-1) ? 0.0 : fz(E, id);

#pragma omp simd collapse(3)
        for (int kk=0; kk<cz; kk++) {
        for (int jj=0; jj<cy; jj++) {
        for (int ii=0; ii<cx; ii++) {
            const int id_da   = (cx*i+ii) + nx_da*(cy*j+jj) + nx_da*ny_da*(cz*k+kk);

            const float ox = ((float)0.5/cx + (float)1.0/cx*ii) - (float)0.5;
            const float oy = ((float)0.5/cy + (float)1.0/cy*jj) - (float)0.5;
            const float oz = ((float)0.5/cz + (float)1.0/cz*kk) - (float)0.5;

            uf[id_da] = _u + _ux*ox + _uy*oy + _uz*oz;
            vf[id_da] = _v + _vx*ox + _vy*oy + _vz*oz;
            wf[id_da] = _w + _wx*ox + _wy*oy + _wz*oz;
            Tf[id_da] = _T + _Tx*ox + _Ty*oy + _Tz*oz;
            Ef[id_da] = _E + _Ex*ox + _Ey*oy + _Ez*oz;
        }
        }
        }
    }
    }
    }
}
