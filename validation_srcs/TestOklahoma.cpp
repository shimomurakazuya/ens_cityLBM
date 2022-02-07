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
#include "IOData.h"
#include "runtime_error.hpp"

#include <random> // ensemble


void  TestOklahoma::Flow()
{
    Field  field(opt_, comm_);

    const bool is_restart = ( opt_.restart_flags_and_step(0) != 0 );
    const auto restart_t =  opt_.restart_flags_and_step(1);

    // arbitral start_minutes for puff_release etc //
    start_minutes0_ = opt_.wrf_start();
    start_minutes_restart_ = start_minutes0_;

    // read restart parameters //
    if(is_restart) {
        field.read_field_params( restart_t );
        start_minutes_restart_ = start_minutes0_ + static_cast<int>(field.parameters().time_now()/60.0);
        if(comm_.is_rank0()) { 
            std::cout << __PRETTY_FUNCTION__ << ": " << "restart load params: "
                << "with time= " << restart_t << " (step); " << start_minutes_restart_ << " (minites)" << std::endl;
        }
    }

    // init mapmem //
    this->ReadMapFile(field.parameters());
    field.init_field_with_map(mapData_);


    // init anyway //
    Init( field.grids(),
          field.tree(),
          field.parameters(),
          field.meshValues0(),
          field.nudgingCoefAdaptation()
          );

    field.
    copy_meshValues(
        field.meshValues1(),
        field.meshValues0(),
        field.grids(),
        field.tree()
        );

    // read restart meshval //
    if (is_restart) {
        if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << " : " << " restart load meshval" << std::endl; }

        field.read_field_meshval( field.optionParser().restart_flags_and_step(1) );
    }

    // read observation again; maybe became non-consistent by restart //
    InitObservationDataOklahoma(
          field.grids(),
          field.tree(),
          field.parameters(),
          field.meshValues0(),
          field.nudgingCoefAdaptation()
          );
    InitObservationDataOklahoma(
          field.grids(),
          field.tree(),
          field.parameters(),
          field.meshValues1(),
          field.nudgingCoefAdaptation()
          );

    if (comm_.is_rank0()) { std::cout << "init_taskID_opt\n"; }
    field.init_taskID_opt(field.taskID(), field.tree(), field.meshValues0() );


#if 1
    /// submit functions ///
//    std::vector< std::function<void(int)> >  cal_funcs;
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


void  TestOklahoma::Init(
    const Grid*                  grids,
    const Tree&                  tree,
    const Parameters&            parameters,
          MeshValue*             meshValues,
    const NudgingCoefAdaptation& nudgingCoefAdaptation
    )
{
    // object //
//    InitObjectMapTokyo(grids, tree, parameters, meshValues);
//    InitObjectMapOklahoma(grids, tree, parameters, meshValues);
    InitObjectMap(grids, tree, parameters, meshValues);

    InitObjectFloarRoof(grids, tree, parameters, meshValues);
#ifndef NO_INIT_DECOBOCO
    InitDecoBoco       (grids, tree, parameters, meshValues);
#endif

    InitPlantDensity   (grids, tree, parameters, meshValues);

    DigitalizeObject   (grids, tree, parameters, meshValues);
    // object //

    InitValue (grids, tree, parameters, meshValues);
    InitValueOklahoma (grids, tree, parameters, meshValues);


    InitObservationDataOklahoma (grids, tree, parameters, meshValues, nudgingCoefAdaptation);
    InitSource (grids, tree, parameters, meshValues);

    // debug : not use (heatflux is initialized in the InitObservationDataOklahoma function) //
//    BoundaryConditions().InitializeHeatFluxObj(0.0, grids, tree, parameters, meshValues);
 
#if REGRESSION_TEST == 1
    if(rank_ == 0) { std::cout << "write binary" << std::endl; }
    IOData(rank_, 0).writeBinaries(grids, tree, meshValues, nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 0);
    MPI_Finalize();
    std::exit(0);
#endif

#if REGRESSION_TEST == 2 || REGRESSION_TEST == 4
    if(rank_ == 0) { std::cout << "write binary" << std::endl; }
    IOData(rank_, 0).writeBinaries(grids, tree, meshValues, nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank_ == 0) { std::cout << "read binary" << std::endl; }
    IOData(rank_, 0).readBinaries(grids, tree, meshValues, nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank_ == 0) { std::cout << "write binary" << std::endl; }
    IOData(rank_, 0).writeBinaries(grids, tree, meshValues, nullptr);
    MPI_Barrier(MPI_COMM_WORLD);
#if REGRESSION_TEST == 2
    MPI_Abort(MPI_COMM_WORLD, 0);
    MPI_Finalize();
    std::exit(0);
#endif
#endif
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
    BoundaryConditions  boundaryConditions(comm_);
    boundaryConditions.InitializeByWRFData   (start_minutes_restart_, grids, tree, parameters, meshValues);
    boundaryConditions.InitializeByGroundData(start_minutes_restart_, grids, tree, parameters, meshValues);
}


void  TestOklahoma::InitObjectFloarRoof(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

    const real dxc = meshValues[0].coordinates().dx();

#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();
        const auto& coordinates = meshValues[lv].coordinates();

        // coordinates xyz //
        const real   dx = coordinates.dx();

        // initialize //
        real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obj  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obj    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obj    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obj    = &meshValues[lv].valueObjLS().w_obj()  [offset];

#pragma omp simd
        for (int id=0; id<nn_cell; id++) {
            // lv //
            const real z_floar = Zmin_;
            const real z_roof  = parameters.coefDomain().z_global_domain_max - dxc*(real)3.0;
//            const real z_roof  = 1000.0;

            const real z_tmp = coordinates.zcell(offset+id);
            const real floar_lv_obj = z_floar - z_tmp;
            const real roof_lv_obj  = z_tmp - z_roof;


            real tmp_lv_obj = lv_obj[id];
            tmp_lv_obj      = FuncObj::overlap_objects(tmp_lv_obj, floar_lv_obj);
            tmp_lv_obj      = FuncObj::overlap_objects(tmp_lv_obj, roof_lv_obj);

            lv_obj[id] = tmp_lv_obj;

            // vel //
            rho_obj[id] = 1.0;
            u_obj  [id] = 0.0;
            v_obj  [id] = 0.0;
            w_obj  [id] = 0.0;
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
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real c_ref   = parameters.c_ref_lbm();

//    const real dxc = meshValues[0].               coordinates().dx();
//    const real dxf = meshValues[DefAMR::LV_MAX-2].coordinates().dx();

#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];

        for (int id=0; id<nn_cell; id++) {
            // lv //
            double tmp_lv_obj = std::floor( (lv_obj[id] + dx)/dx )*dx - 0.5*dx;

            lv_obj[id] = tmp_lv_obj;
        }
    }
}

void TestOklahoma::ReadMapFile(const Parameters& parameters) {
    // Oklahoma
    constexpr real dx_map = 1;
    mapData_.load(Foldernames::input_folder + "/msm/building/buildings1m_1200_1500.txt");
    mapData_.set_geo(
            dx_map, // dx
            parameters.coefDomain().x_global_domain_center() - mapData_.mx()*dx_map/2 + dx_map/2, // x0
            parameters.coefDomain().y_global_domain_center() - mapData_.my()*dx_map/2 + dx_map/2, // y0
            parameters.coefDomain().x_global_domain_min + (real)16.0 + bc_layer_min_[0], // west
            parameters.coefDomain().x_global_domain_max - (real)16.0 - bc_layer_max_[0], // east
            parameters.coefDomain().y_global_domain_min + (real)16.0 + bc_layer_min_[1], // south
            parameters.coefDomain().y_global_domain_max - (real)16.0 - bc_layer_max_[1]  // north
            );

//    // Tokyo
//    constexpr real dx_map = 1;
//    constexpr real dxc = 4; // coarse?
//    mapData_.load(Foldernames::input_folder + "/" + Filenames::Map_name0);
//    mapData_.set_geo(
//            dx_map, // dx
//            -1500, -1500, // x0, y0
//             dxc + parameters.coefDomain().x_global_domain_min, // west
//            -dxc + parameters.coefDomain().x_global_domain_max, // east
//             dxc + parameters.coefDomain().y_global_domain_min, // south
//            -dxc + parameters.coefDomain().y_global_domain_max  // north
//            );

//    // 国土地理院　地理院地図3D
//    constexpr real dx_map = 15.625;
//    //constexpr real dxc = 4; // coarse?
//    mapData_.load(Foldernames::input_folder + "/map/fukushima1F257x257.map");
//    mapData_.set_geo(
//            dx_map,
//            -2000, -2000, // offset
//            -1800, // west
//             1800, // east
//            -1800, // south
//             1800  // north
//    );
}

void TestOklahoma::InitObjectMap(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const auto& nodes = tree.nodes();
    const auto& n_leaf = nodes.size();
    const auto& c_ref = parameters.c_ref_lbm();

#pragma omp parallel for
    for(std::intptr_t l=0; l<n_leaf; l++) { 
        // <-- aliases
        const auto& node = nodes.at(l);
        const auto& lv = node.level();
        const auto& offset = node.mesh_offset();
        const auto& coor = meshValues[lv].coordinates();
        auto& obj = meshValues[lv].valueObjLS();
        std::function<double(int)>
            x = [&](int id) { return coor.xcell(offset + id); },
            y = [&](int id) { return coor.ycell(offset + id); },
            z = [&](int id) { return coor.zcell(offset + id); };
        // -->
        mapData_.levelset_map(
                obj.lv_obj() + offset,
                x, y, z,
                DefAMR::NX_LEAF,
                lv, coor.dx()
                );
    }
}


void  TestOklahoma::InitPlantDensity(
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
        const auto& coordinates = meshValues[lv].coordinates();

        // coordinates xyz //
        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        const real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
              real*  pad_obj  = &meshValues[lv].valueObjLS().pad_obj()[offset];

        const real  c_ref  = parameters.c_ref_lbm();
        for (int id=0; id<nn_cell; id++) {
            const real x_tmp = coordinates.xcell(offset+id);
            const real y_tmp = coordinates.ycell(offset+id);
            const real z_tmp = coordinates.zcell(offset+id);

            // park
            #ifndef ENSEMBLE_HEIGHT_PLANT_TEST
            constexpr real bbox_zmax = 16.0;
            #else
            #warning enabling ENSEMBLE_HEIGHT_PLANT_TEST
            const real bbox_zmax = 8.0 * comm_.color_ens_offseted(); // 0 -- 64 m at 9 ensembles
            #endif
//            const real box_min[] = { -305.0, -382.0,  0.0 };
            const real box_min[] = { -305.0, -402.0,  0.0 };
            const real box_max[] = {  -45.0, -182.0, bbox_zmax };

            const real xyz[] = { x_tmp, y_tmp, z_tmp };
            const real box_lv_obj = FuncObj::length_from_box(xyz, box_min, box_max);

            if ( FuncObj::is_fluid(lv_obj[id]) && FuncObj::is_obj(box_lv_obj) ) {
                pad_obj[id] = air_property::Cd_pad;
            }

        }
    }
}


void  TestOklahoma::InitDecoBoco(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;

    const double dxc = meshValues[0].coordinates().dx();
    const double dxf = meshValues[DefAMR::LV_MAX-1].coordinates().dx();

#if EXPECT_DX_MESH == 2 || EXPECT_DX_MESH == 1
    // /data/g5/a170085/research/citylbm/validation/2m_v0
    const double width_decoboco[3] = { 8.0, 8.0, 24.0 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*8, width_decoboco[1]*8, width_decoboco[2]*4 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*6, width_decoboco[1]*6, width_decoboco[2]*4 };
    const double pitch_decoboco[3] = { width_decoboco[0]*4, width_decoboco[1]*4, width_decoboco[2]*4 };
#elif EXPECT_DX_MESH == 4
    // 4m //
//    const double width_decoboco[3] = { 16.0, 16.0, 8.0 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*8, width_decoboco[1]*8, width_decoboco[2]*4 };

//    const double width_decoboco[3] = { 8.0, 8.0, 8.0 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*4, width_decoboco[1]*4, width_decoboco[2]*4 };

//    const double width_decoboco[3] = { 12.0, 12.0, 12.0 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*3.5, width_decoboco[1]*3.5, width_decoboco[2]*3.5 };

    // /data/g5/a170085/research/citylbm/validation/4m_v0
    const double width_decoboco[3] = { 8.0, 8.0, 24.0 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*8, width_decoboco[1]*8, width_decoboco[2]*4 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*6, width_decoboco[1]*6, width_decoboco[2]*4 };
    const double pitch_decoboco[3] = { width_decoboco[0]*4, width_decoboco[1]*4, width_decoboco[2]*4 };

//    // /data/g5/a170085/research/citylbm/validation/4m_v1
//    const double width_decoboco[3] = { 16.0, 16.0, 16.0 };
//    const double pitch_decoboco[3] = { width_decoboco[0]*4, width_decoboco[1]*4, width_decoboco[2]*4 };
#else
#error decoboco not defined on this EXPECT_DX_MESH
#endif

    const int nx_deco = fabs(parameters.coefDomain().x_global_domain_max - parameters.coefDomain().x_global_domain_min)/pitch_decoboco[0];
    const int ny_deco = fabs(parameters.coefDomain().y_global_domain_max - parameters.coefDomain().y_global_domain_min)/pitch_decoboco[1];

    for (int jj=0; jj<ny_deco; jj++) {
    for (int ii=0; ii<nx_deco; ii++) {
        // deco boco //
        double _di = 0;
        double _dj = 0.5 * (ii%2);
        #ifdef ENSEMBLE_SPIKE
        if(comm_.ensemble_id() != 0) {
            constexpr int ens_mv_count = 1;
            constexpr double mv_di = 0.25;
            const auto ens_seed = ii + jj*nx_deco + comm_.ensemble_id()*nx_deco*ny_deco;
            auto&& engine = std::mt19937(ens_seed);
            auto&& dist = std::uniform_int_distribution<int>(-ens_mv_count, ens_mv_count);

            volatile const double mv0 = dist(engine) * mv_di;
            volatile const double mv1 = dist(engine) * mv_di;
            _di += mv0;
            _dj += mv1;
        }
        #endif
        const double _ii = double(ii) + _di;
        const double _jj = double(jj) + _dj;
        const double box_min[3] = { parameters.coefDomain().x_global_domain_min + _ii*pitch_decoboco[0],
                                    parameters.coefDomain().y_global_domain_min + _jj*pitch_decoboco[1],
                                    parameters.coefDomain().z_global_domain_min - 1.0   };
        const double box_max[3] = { box_min[0] + width_decoboco[0],
                                    box_min[1] + width_decoboco[1],
                                                 width_decoboco[2]  };
        #pragma omp parallel for 
        for (int l=0; l<n_leaf; l++) {
            // deco boco //
            const int  lv     = tree.nodes(l)->level();
            const int  offset = tree.nodes(l)->mesh_offset();
            const auto& coordinates = meshValues[lv].coordinates();

            const double dx = meshValues[lv].coordinates().dx();

            const double _x_tmp = coordinates.xnode(offset);
            const double _y_tmp = coordinates.ynode(offset);
            const double _z_tmp = coordinates.znode(offset);

            if ( _z_tmp < box_min[2] - dx*2.5 || _z_tmp + dx*DefAMR::NX_LEAF > box_max[2] + dx*2.5 ) { continue; }
            if ( _x_tmp < box_min[0] - dx*(real)2.5 || _x_tmp + dx*DefAMR::NX_LEAF > box_max[0] + dx*(real)2.5 ) { continue; }
            if ( _y_tmp < box_min[1] - dx*(real)2.5 || _y_tmp + dx*DefAMR::NX_LEAF > box_max[1] + dx*(real)2.5 ) { continue; }

            // initialize //
            auto*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
            auto*  rho_obj  = &meshValues[lv].valueObjLS().rho_obj()[offset];
            auto*  u_obj    = &meshValues[lv].valueObjLS().u_obj()  [offset];
            auto*  v_obj    = &meshValues[lv].valueObjLS().v_obj()  [offset];
            auto*  w_obj    = &meshValues[lv].valueObjLS().w_obj()  [offset];

            for (int id=0; id<nn_cell; id++) {
                // lv //
                const double x_tmp = coordinates.xcell(offset + id);
                const double y_tmp = coordinates.ycell(offset + id);
                const double z_tmp = coordinates.zcell(offset + id);

                const double xyz[3] = { x_tmp, y_tmp, z_tmp };

                const double  cal_min_x1 = parameters.coefDomain().x_global_domain_min + dxc*4.0 + bc_layer_min_[0];
                const double  cal_max_x1 = parameters.coefDomain().x_global_domain_max - dxc*4.0 - bc_layer_max_[0];

                const double  cal_min_y1 = parameters.coefDomain().y_global_domain_min + dxc*4.0 + bc_layer_min_[1];
                const double  cal_max_y1 = parameters.coefDomain().y_global_domain_max - dxc*4.0 - bc_layer_max_[1];

                const double  cal_min_x2 = parameters.coefDomain().x_global_domain_min + dxc*4.0 + bc_layer_min_[0] + 128.0;
                const double  cal_max_x2 = parameters.coefDomain().x_global_domain_max - dxc*4.0 - bc_layer_max_[0] - 128.0;

                const double  cal_min_y2 = parameters.coefDomain().y_global_domain_min + dxc*4.0 + bc_layer_min_[1] + 128.0;
                const double  cal_max_y2 = parameters.coefDomain().y_global_domain_max - dxc*4.0 - bc_layer_max_[1] - 128.0;

                // search no spike region //
                /// 1: outer boundary ///
                if ( x_tmp<cal_min_x1 || x_tmp>cal_max_x1 || y_tmp<cal_min_y1 || y_tmp>cal_max_y1 ) {
                    continue;
                }
                /// 2: city centre ///
                constexpr double xw_bb0 = 320;
                constexpr double yw_bb0 = 620;
                constexpr double x0_bb = -10;
                constexpr double y0_bb = 60;
                #ifndef ENSEMBLE_BB_SPIKE_TEST
                constexpr double xw_bb = xw_bb0;
                constexpr double yw_bb = yw_bb0;
                #else
                #warning enabling ENSEMBLE_BB_SPIKE_TEST
                constexpr double dxbb = 40;
                const int seed = comm_.color_ens_offseted();
                const double xw_bb = xw_bb0 + seed*dxbb;
                const double yw_bb = yw_bb0 + seed*dxbb;
                #endif
                if (fabs(x_tmp - x0_bb) < xw_bb && fabs(y_tmp - y0_bb) < yw_bb) { 
                    continue;
                }

                // set spike on reamined region
                const double box_lv_obj = FuncObj::length_from_box(xyz, box_min, box_max);

                double tmp_lv_obj = lv_obj[id];
                lv_obj[id] = FuncObj::overlap_objects(tmp_lv_obj, box_lv_obj);
            } // foreach cell
        } // foreach leaf
    }
    } // foreach decoboco
}


void  TestOklahoma::InitSource(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
        if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    const real source_xyz[] = { -50.0, -176.0, 2.0 };
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;


    const int n_scalars = meshValues[0].n_scalars();
    std::vector<int>  num_source(n_scalars, 0);
    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();
        const auto& coordinates = meshValues[lv].coordinates();

        // coordinates xyz //
        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        const real*  lv_obj     = &meshValues[lv].valueObjLS().lv_obj()     [offset];
              real*  sc_scalar  = &meshValues[lv].valueObjLS().sc_scalar () [offset];
        const int nn_max = meshValues[lv].nn_max();

        const real  c_ref  = parameters.c_ref_lbm();
        for (int id=0; id<nn_cell; id++) {
            const real x_tmp = coordinates.xcell(offset + id);
            const real y_tmp = coordinates.ycell(offset + id);
            const real z_tmp = coordinates.zcell(offset + id);

//            const real box_min[] = { -(real)50.0-(real)0.5*dx, -(real)150.0-(real)0.5*dx, (real)0.0 };
//            const real box_max[] = { -(real)50.0+(real)0.5*dx, -(real)150.0+(real)0.5*dx,       dx  };

            const real box_min[] = { source_xyz[0]-(real)0.5*dx, source_xyz[1]-(real)0.5*dx, source_xyz[2]-(real)0.5*dx };
            const real box_max[] = { source_xyz[0]+(real)0.5*dx, source_xyz[1]+(real)0.5*dx, source_xyz[2]+(real)0.5*dx };
//            const real box_min[] = { source_xyz[0]-(real)1.0*dx, source_xyz[1]-(real)1.0*dx, source_xyz[2]-(real)1.0*dx };
//            const real box_max[] = { source_xyz[0]+(real)1.0*dx, source_xyz[1]+(real)1.0*dx, source_xyz[2]+(real)1.0*dx };

//            // tokyo //
//            const real box_min[] = { -(real)2.0*dx, -(real)2.0*dx, (real)10.0 + (real)8.0+(real)0.0*dx };
//            const real box_max[] = { +(real)2.0*dx, +(real)2.0*dx, (real)10.0 + (real)8.0+(real)4.0*dx  };

            const real xyz[] = { x_tmp, y_tmp, z_tmp };
            const real box_lv_obj = FuncObj::length_from_box(xyz, box_min, box_max);

            if ( FuncObj::is_fluid(lv_obj[id]) && FuncObj::is_obj(box_lv_obj) ) {
                for(int n=0; n<n_scalars; n++) {
                    sc_scalar[id + n*nn_max] = (real)1.0/(dx*dx*dx);
                    if ( tree.nodes(i)->nodeCalFlags().Cal() ) { num_source[n]++; }
                }
            }

//            for (int jjj=0; jjj<2; jjj++) {
//            for (int iii=0; iii<2; iii++) {
//                const real dx_sc = (real)4.0;
//                const real xyz_tmp[] = { x_tmp - dx_sc*(real)2.0*(iii-(real)0.5), y_tmp - dx_sc*(real)2.0*(jjj-(real)0.5), z_tmp };
//                const real box_lv_obj = FuncObj::length_from_box(xyz_tmp, box_min, box_max);
//
//                if ( FuncObj::is_fluid(lv_obj[id]) && FuncObj::is_obj(box_lv_obj) ) {
//                    sc_scalars[iii+2*jjj][id] = (real)1.0/(dx*dx*dx);
//                    if ( tree.nodes(i)->nodeCalFlags().Cal() ) { num_sources[iii+2*jjj]++; }
//                }
//            }
//            }

        }

    }

    std::vector<int>  num_source_g(n_scalars, 0);
    const auto comm_ens = comm_.col_vector().comm();
    MPI_Allreduce(num_source.data(), num_source_g.data(), n_scalars, MPI_INT, MPI_SUM, comm_ens);

//    std::cout << rank_ << " : num_source, num_source_global = " << num_source << ", " << num_source_g << std::endl;
    for(int n=0; n<n_scalars; n++) {
        std::cout << comm_.world().rank() << " : scalar id, num_source, num_source_global = " << n << ", " << num_source.at(n) << ", " << num_source_g.at(n) << std::endl;
    }

    if (comm_.is_rank0()) { 
        for(int n=0; n<n_scalars; n++) {
            std::cout << "0" << " : " << 0 << " : scalar id, num_source_global = " << n << ", " << num_source_g.at(n) << std::endl; 
        }
    }
    for(auto ni_source_g : num_source_g) {
        runtime_assert(ni_source_g > 0, "InternalError: source not found in any rank");
    }

    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();
        const int nn_max = meshValues[lv].nn_max();

        for(int n=0; n<n_scalars; n++) {
            real*  sc_scalar    = &meshValues[lv].valueObjLS().sc_scalar ()  [offset + n * nn_max];

            for (int j=0; j<nn_cell; j++) {
                sc_scalar[j] = sc_scalar[j]/(float)num_source_g.at(n);
            }
        }
//        for (int jjj=0; jjj<2; jjj++) {
//        for (int iii=0; iii<2; iii++) {
//            for (int j=0; j<nn_cell; j++) {
//                sc_scalars[iii+2*jjj][j] = sc_scalars[iii+2*jjj][j]/(float)num_sources_g[iii+2*jjj];
//            }
//        }
//        }

    }
}


void  TestOklahoma::InitObservationDataOklahoma(
    const Grid*                  grids,
    const Tree&                  tree,
    const Parameters&            parameters,
          MeshValue*             meshValues,
    const NudgingCoefAdaptation& nudgingCoefAdaptation
    )
{
    const int time_minutes_now = start_minutes_restart_;
    BoundaryConditions  boundaryConditions(comm_);
    boundaryConditions.ReadWRFData   (time_minutes_now, grids, tree, parameters, meshValues, nudgingCoefAdaptation.coef_nudging());
    boundaryConditions.ReadGroundData(time_minutes_now, grids, tree, parameters, meshValues);

    // swap: rho_obj <-> rhon_obj, .. etc //
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        meshValues[lv].valueObjLS().swap_objs();
    }
    boundaryConditions.ReadWRFData   (time_minutes_now+1, grids, tree, parameters, meshValues, nudgingCoefAdaptation.coef_nudging());
    boundaryConditions.ReadGroundData(time_minutes_now+1, grids, tree, parameters, meshValues);

    // swap: rho_obj <-> rhon_obj, .. etc //
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        meshValues[lv].valueObjLS().swap_objs();
    }
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
