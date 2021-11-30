#include "TestFlowAroundCube.h"
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


void  TestFlowAroundCube::Flow()
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
    field.init_field_wo_map();


    // init anyway //
    Init(   field.grids(),
            field.tree(),
            field.parameters(),
            field.meshValues0() );

    field.copy_meshValues(
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


void  TestFlowAroundCube::Init(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{

    InitObjectFloarRoof(grids, tree, parameters, meshValues, box_min_, box_max_);
    InitValue (grids, tree, parameters, meshValues);
    InitSource (grids, tree, parameters, meshValues);
    //initBCAroundCube(grids, tree, parameters, meshValues);
}


void  TestFlowAroundCube::InitValue(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues
    )
{
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const auto& coordinates = meshValues[lv].coordinates();

        // initialize //
        // NS //
        real*  u   = &meshValues[lv].valueNS().u()  [offset];
        real*  v   = &meshValues[lv].valueNS().v()  [offset];
        real*  w   = &meshValues[lv].valueNS().w()  [offset];
        real*  rho = &meshValues[lv].valueNS().rho()[offset];

        real*  T   = &meshValues[lv].valueNS().T()  [offset];

        const real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        const real  rho_tmp = 1.0;
        const real  c_ref  = parameters.c_ref_lbm();
        for (int id=0; id<nn_cell; id++) {
            const real zc = coordinates.zcell(offset + id);
            //const real  u_tmp   =  0.0;
            const real  u_tmp   = FuncObj::is_fluid(lv_obj[id]) ? FuncFlowAroundCube::inflow_velocity(zc) : 0.0;
            const real  v_tmp   =  0.0;
            const real  w_tmp   =  0.0;
            const real  T_tmp   = fluid_property::Temperature0;

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


void  TestFlowAroundCube::InitObjectFloarRoof(
    const Grid*         grids,
    const Tree&         tree,
    const Parameters&   parameters,
          MeshValue*    meshValues,
    const real          box_min[],
    const real          box_max[]

    )
{
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;
    const real  c_ref  = parameters.c_ref_lbm();

    const real dxc = meshValues[0].coordinates().dx();

#pragma omp parallel for
    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const auto& coordinates = meshValues[lv].coordinates();

        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        real*  lv_obj   = &meshValues[lv].valueObjLS().lv_obj() [offset];
        real*  rho_obj  = &meshValues[lv].valueObjLS().rho_obj()[offset];
        real*  u_obj    = &meshValues[lv].valueObjLS().u_obj()  [offset];
        real*  v_obj    = &meshValues[lv].valueObjLS().v_obj()  [offset];
        real*  w_obj    = &meshValues[lv].valueObjLS().w_obj()  [offset];
        auto* bcTypes_f = &meshValues[lv].valueObjLS().bcTypes_f()[offset];
        //auto* bcTypes_T = &meshValues[lv].valueObjLS().bcTypes_T()[offset];
        real* bc_weight = &meshValues[lv].valueObjLS().dirichlet_weight()[offset];

        real*  rhon_obj     = &meshValues[lv].valueObjLS().rhon_obj() [offset];
        real*  un_obj       = &meshValues[lv].valueObjLS().un_obj() [offset];
        real*  vn_obj       = &meshValues[lv].valueObjLS().vn_obj() [offset];
        real*  wn_obj       = &meshValues[lv].valueObjLS().wn_obj() [offset];
        //real*  Tn_obj       = &meshValues[lv].valueObjLS().Tn_obj() [offset];
        real*  scalarn_obj  = &meshValues[lv].valueObjLS().scalarn_obj() [offset];




#pragma omp simd
        for (int id=0; id<nn_cell; id++) {
            const real xc = coordinates.xcell(offset + id);
            const real yc = coordinates.ycell(offset + id);
            const real zc = coordinates.zcell(offset + id);
            const real  xyz[] = { xc, yc, zc};
            // lv //
            const real x_floar = parameters.coefDomain().x_global_domain_min ;
            const real x_roof  = parameters.coefDomain().x_global_domain_max ;

            const real x_floar_lv_obj = x_floar - xc;
            const real x_roof_lv_obj  = xc - x_roof;
            const real z_floar = 0;
            const real z_roof  = parameters.coefDomain().z_global_domain_max ;

            const real floar_lv_obj = z_floar - zc;
            const real roof_lv_obj  = zc - z_roof;
            const real box_lv_obj   =   FuncObj::length_from_box(xyz, box_min, box_max);
            real tmp_lv_obj = box_lv_obj;

            const real bottom_min[3] = { parameters.coefDomain().x_global_domain_min, parameters.coefDomain().y_global_domain_min, parameters.coefDomain().z_global_domain_min + (real)dxc };
            const real bottom_max[3] = { parameters.coefDomain().x_global_domain_max, parameters.coefDomain().y_global_domain_max, z_floar };
            const real bottom_lv_obj = FuncObj::length_from_box(xyz, bottom_min, bottom_max);

            tmp_lv_obj = FuncObj::overlap_objects(tmp_lv_obj, roof_lv_obj);
            tmp_lv_obj = FuncObj::overlap_objects(tmp_lv_obj, floar_lv_obj);
            tmp_lv_obj = FuncObj::overlap_objects(tmp_lv_obj, bottom_lv_obj);

            // roughness blocks
            tmp_lv_obj = levelset_from_Inflow_fin(tmp_lv_obj, xyz);

            lv_obj[id] = tmp_lv_obj;
            bcTypes_f[id] = BCTypes::CalRegion;
            // const double lv_tmp  = FuncObj::length_from_box(xyz, box_min, box_max);
            // const double nlv_obj = floor( lv_tmp/dx )*dx + 0.5*dx;

            // vel //
            u_obj  [id] = 0.0;
            v_obj  [id] = 0.0;
            w_obj  [id] = 0.0;
            un_obj  [id] = 0.0;
            vn_obj  [id] = 0.0;
            wn_obj  [id] = 0.0;
            // setting roof velocity
            if ( fabs(roof_lv_obj) <= 2.0*dx ) {
                bcTypes_f[id] = BCTypes::BCRegion;
                bc_weight[id] = 1.0;
                u_obj[id] = FuncFlowAroundCube::inflow_velocity(zc) / c_ref;
                v_obj[id] = 0.0;
                w_obj[id] = 0.0;
                un_obj[id] = FuncFlowAroundCube::inflow_velocity(zc) / c_ref;
                vn_obj[id] = 0.0;
                wn_obj[id] = 0.0;
            }

            if ( xyz[2] > xyz_max_[2] - bc_layer_max_[2] ) {
                bcTypes_f[id] = BCTypes::BCRegion;
                //bcTypes_T[id] = BCTypes::BCRegion;

                //bc_weight[j] = 0.05;
                bc_weight[id] = 0.5;
                un_obj[id] = FuncFlowAroundCube::inflow_velocity(zc)/ c_ref;
                u_obj[id] = FuncFlowAroundCube::inflow_velocity(zc) / c_ref;
            }

            // setting nudging area
            if((fabs(x_roof_lv_obj) <=  2.5 || fabs(x_floar_lv_obj) <= 2.5) && 0 <= zc ){
                bcTypes_f[id] = BCTypes::BCRegion;
                bc_weight[id] = 0.5;
                u_obj[id] = FuncFlowAroundCube::inflow_velocity(zc) / c_ref;
                v_obj[id] = 0.0;
                w_obj[id] = 0.0;
                un_obj[id] = FuncFlowAroundCube::inflow_velocity(zc) / c_ref;
                vn_obj[id] = 0.0;
                wn_obj[id] = 0.0;
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
    if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const real source_xyz[] = { -0.15, 0.0, 0.05 };
    const int  n_leaf  = tree.number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;

    int  num_source = 0;
    for (int i=0; i<n_leaf; i++) {
        const int  lv     = tree.nodes(i)->level();
        const int  offset = tree.nodes(i)->mesh_offset();

        // coordinates xyz //
        const auto& coordinates = meshValues[lv].coordinates();
        const real   dx = meshValues[lv].coordinates().dx();

        // initialize //
        const real*  lv_obj     = &meshValues[lv].valueObjLS().lv_obj()     [offset];
        real*  sc_scalar  = &meshValues[lv].valueObjLS().sc_scalar () [offset];

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


            const real xyz[] = { x_tmp, y_tmp, z_tmp };
            const real box_lv_obj = FuncObj::length_from_box(xyz, box_min, box_max);

            if ( FuncObj::is_fluid(lv_obj[id]) && FuncObj::is_obj(box_lv_obj) ) {
               // sc_scalar[id] = (real)1.0/(dx*dx*dx);
                sc_scalar[id] = sc_dscalar_ / (dx*dx*dx) * 1000.0; // ppm = mg/L : kg / m3 = 1000 * 1000 mg / 1000 L //
                if ( tree.nodes(i)->nodeCalFlags().Cal() ) { num_source++; }
            }


            for (int jjj=0; jjj<2; jjj++) {
                for (int iii=0; iii<2; iii++) {
                    const real dx_sc = (real)4.0;
                    const real xyz_tmp[] = { x_tmp - dx_sc*(real)2.0*(iii-(real)0.5), y_tmp - dx_sc*(real)2.0*(jjj-(real)0.5), z_tmp };
                    const real box_lv_obj = FuncObj::length_from_box(xyz_tmp, box_min, box_max);

                }
            }
        }
    }

    MPI_Request requ[1];
    int  num_source_g = 0;
    const auto comm_ens = comm_.col_vector().comm();
    MPI_Iallreduce(&num_source,     &num_source_g,     1, MPI_INT, MPI_SUM, comm_ens, &requ[0]);
    MPI_Waitall(1, requ, MPI_STATUSES_IGNORE);

    // std::cout << comm_.world().rank() << " : num_source, num_source_global = " << num_source << ", " << num_source_g << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_.is_rank0()) { std::cout << "0" << " : " << 0 << " :  num_source_global = " << num_source_g << std::endl; }
    runtime_assert(num_source_g > 0, "InternalError: source not found in any rank");


    for (int l=0; l<n_leaf; l++) {
        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        real*  sc_scalar    = &meshValues[lv].valueObjLS().sc_scalar ()  [offset];
            for (int id=0; id<nn_cell; id++) {
                sc_scalar[id] = sc_scalar[id]/(float)num_source_g;
            }

    }
}


real  TestFlowAroundCube::levelset_from_Inflow_fin(const real lv_obj, const real xyz[]) const
{
    const int  num_fins = num_fins_;
    const real length = (xyz_max_[1] - xyz_min_[1]) / (num_fins + 1);

    real  tmp_lv_obj = lv_obj;
    for (int i=0; i<num_fins; i++) {
        real mv_x = 0, mv_y = 0;
        #ifdef ENSEMBLE_SPIKE
        if(comm_.ensemble_id() > 0) { // except ensemble_id == 0 (control-run)
            const real ens_length_x = 10 * length;
            const auto ens_seed = i; // torima（フィン番号iは全MPIでユニークか？）
            const auto ens_seed2 = comm_.ensemble_id();
            auto&& engine = std::mt19937(std::mt19937(ens_seed)() + std::mt19937(ens_seed2)()); // torima mt19937 chain
            auto&& dist = std::uniform_real_distribution<real>(-1, 1);
            //mv_y = dist(engine);
            mv_x = std::abs(dist(engine) * ens_length_x);
        }
        #endif

        // box with ensemble shift
        const real  box_min[3] = { fin_pos_              + mv_x, xyz_min_[1] + length*(i+1) - fin_width_*(real)0.5 + mv_y, xyz_min_[2]              };
        const real  box_max[3] = { fin_pos_ + fin_width_ + mv_x, xyz_min_[1] + length*(i+1) + fin_width_*(real)0.5 + mv_y, xyz_min_[2] + fin_height_ };


        const real box_lv_obj = FuncObj::length_from_box(xyz, box_min, box_max);

        tmp_lv_obj = FuncObj::overlap_objects(tmp_lv_obj, box_lv_obj);
    }

    return  tmp_lv_obj;
}
