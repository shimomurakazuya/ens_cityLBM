#include "LBMCalculation.h"
#include "defineCal.h"
#include "defineAMR.h"
#include "defineFluidProperty.h"
#include "FuncScheduler.h"
#include "FuncLBMKernel.h"
#include "FuncLoop.h"
#include "LBMCalculationGPU.h"
#include "LBMCalculationCPU.h"

#include "Index.h"
#include "IndexLBM.h"
#include "BoundaryConditions.h"


void  LBMCalculation::
LBM_data_assimilation(int t, Field& field)
{
    const int  time_minutes_start = 360;

    const real time    = field.parameters().time_now();
    const int  time_minutes = int(time/60.0) + time_minutes_start;
    static int time_minutes_now;
    static int is_first;

//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if (rank == 0) { std::cout << "t, time, time_minutes, time_minutes_now = " << t << ", " << time << ", " << time_minutes << ", " << time_minutes_now << std::endl; }

    if (time_minutes != time_minutes_now || is_first == 0) {
        time_minutes_now = time_minutes;
        is_first = 1;

        // read: rho_obj, u_obj, ... etc. //
        BoundaryConditions  boundaryConditions;
        boundaryConditions.ReadWRFData   (time_minutes_now+1, field.grids(), field.tree(), field.parameters(), field.meshValues0());
        boundaryConditions.ReadGroundData(time_minutes_now+1, field.grids(), field.tree(), field.parameters(), field.meshValues0());

        // swap: rho_obj <-> rhon_obj, .. etc //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            field.meshValues0()[lv].valueObjLS().swap_objs();
        }
//        boundaryConditions.ReadWRFData   (time_minutes_now, field.grids(), field.tree(), field.parameters(), field.meshValues0());
//        boundaryConditions.ReadGroundData(time_minutes_now, field.grids(), field.tree(), field.parameters(), field.meshValues0());
    }
}


void  LBMCalculation::
LBM_incompressible_flow(int t, Field& field)
{
#ifdef TEMPORAL_BLOCKING_
//    LBM_incompressible_flow_w_tb (t, field);
    LBM_incompressible_flow_wo_tb (t, field);
#else
    LBM_incompressible_flow_wo_tb (t, field);
#endif
}



void  LBMCalculation::
LBM_incompressible_flow_w_tb(int t, Field& field)
{
    const int  lv_max = DefAMR::LV_MAX;

    // functions //
    std::function<void(int)>  func_cal_with_halo = [this, &field](int lv){
        field.swap_meshValue(lv);

        field.funcTimeInfo().StartTimer("LBM_to_euler_variables");
        LBM_to_euler_variables(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("LBM_to_euler_variables");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_to_euler_variables", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


        field.funcTimeInfo().StartTimer("ScalarAdvection");
        ScalarAdvection(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("ScalarAdvection");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ScalarAdvection", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


        field.funcTimeInfo().StartTimer("LBM_stream_collision");
        LBM_stream_collision_opt(lv, field,
                                 field.taskID().id_task_array_with_halo_wo_bc(lv), field.taskID().num_task_array_with_halo_wo_bc(lv),
                                 field.taskID().id_task_array_with_halo_w_bc (lv), field.taskID().num_task_array_with_halo_w_bc (lv) );
        field.funcTimeInfo().StopTimer("LBM_stream_collision");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_stream_collision", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


        // MPI //
        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        ValuesCommunications(lv, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //


        // BC //
        field.funcTimeInfo().StartTimer("boundary_conditions");
        boundary_conditions(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("boundary_conditions");
        field.funcTimeInfo().SubmitElapsedTimeInfo("boundary_conditions", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());

    };


    std::function<void(int)>  func_L2F = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("L2F");
        Val_L2F(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        LBM_L2F(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        field.funcTimeInfo().StopTimer("L2F");
        field.funcTimeInfo().SubmitElapsedTimeInfo("L2F", lv, field.taskID().id_tasks_L2F(lv).size(), field.parameters().step_now());
    };


    std::function<void(int)>  func_empty = [this, &field](int lv){};


    std::function<void(int)>  func_F2L = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("F2L");
        Val_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        LBM_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        field.funcTimeInfo().StopTimer("F2L");
        field.funcTimeInfo().SubmitElapsedTimeInfo("F2L", lv, field.taskID().id_tasks_F2L(lv).size(), field.parameters().step_now());
    };


    // cal //
    field.timerSimple().start_cal_timer(); // timer simple //
    FuncScheduler::launcher( lv_max, func_cal_with_halo, func_L2F, func_empty, func_F2L );
    field.timerSimple().stop_cal_timer();  field.timerSimple().add_cal_timer(); // timer simple //
}



void  LBMCalculation::
LBM_incompressible_flow_wo_tb(int t, Field& field)
{
    const int  lv_max = DefAMR::LV_MAX;

    // functions //
    std::function<void(int)>  func_cal = [this, &field](int lv){
        field.swap_meshValue(lv);

        field.funcTimeInfo().StartTimer("LBM_to_euler_variables");
        LBM_to_euler_variables(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("LBM_to_euler_variables");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_to_euler_variables", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


        field.funcTimeInfo().StartTimer("ScalarAdvection");
        ScalarAdvection(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("ScalarAdvection");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ScalarAdvection", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


        field.funcTimeInfo().StartTimer("LBM_stream_collision");
        LBM_stream_collision_opt(lv, field,
                                 field.taskID().id_task_array_with_halo_wo_bc(lv), field.taskID().num_task_array_with_halo_wo_bc(lv),
                                 field.taskID().id_task_array_with_halo_w_bc (lv), field.taskID().num_task_array_with_halo_w_bc (lv) );
        field.funcTimeInfo().StopTimer("LBM_stream_collision");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_stream_collision", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


        // MPI //
        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        ValuesCommunications(lv, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //


        // BC //
        field.funcTimeInfo().StartTimer("boundary_conditions");
        boundary_conditions(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("boundary_conditions");
        field.funcTimeInfo().SubmitElapsedTimeInfo("boundary_conditions", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());

    };


    std::function<void(int)>  func_L2F = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("L2F");
        Val_L2F(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        LBM_L2F(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        field.funcTimeInfo().StopTimer("L2F");
        field.funcTimeInfo().SubmitElapsedTimeInfo("L2F", lv, field.taskID().id_tasks_L2F(lv).size(), field.parameters().step_now());
    };


    std::function<void(int)>  func_L2FA = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("L2FA");
        Val_L2FA(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        LBM_L2FA(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        field.funcTimeInfo().StopTimer("F2LA");
        field.funcTimeInfo().SubmitElapsedTimeInfo("L2FA", lv, field.taskID().id_tasks_L2F(lv).size(), field.parameters().step_now());
    };


    std::function<void(int)>  func_F2L = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("F2L");
        Val_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        LBM_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        field.funcTimeInfo().StopTimer("F2L");
        field.funcTimeInfo().SubmitElapsedTimeInfo("F2L", lv, field.taskID().id_tasks_F2L(lv).size(), field.parameters().step_now());
    };


    // cal //
    field.timerSimple().start_cal_timer(); // timer simple //
    FuncScheduler::launcher   ( lv_max, func_cal, func_L2F, func_L2FA, func_F2L );
    field.timerSimple().stop_cal_timer();  field.timerSimple().add_cal_timer(); // timer simple //
}


void  LBMCalculation::
LBM_stream_collision_opt( const int lv, Field& field,
                          const int*  id_tasks_wo_bc, const int  num_tasks_wo_bc,
                          const int*  id_tasks_w_bc,  const int  num_tasks_w_bc)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();


    const real*  f_lbm  = field.meshValue    (lv).valueLBM().f_lbm();
          real*  fn_lbm = field.meshValue_new(lv).valueLBM().f_lbm();

          real*  sgs_vis = field.meshValue_new(lv).valueNS().sgs_vis();

    const real*  T      = field.meshValue    (lv).valueNS().T();
    const real*  Tn     = field.meshValue_new(lv).valueNS().T();

    const real*  lv_obs = field.meshValue(lv).valueObjLS().lv_obj();
    const real*  u_obs  = field.meshValue(lv).valueObjLS().u_obj();
    const real*  v_obs  = field.meshValue(lv).valueObjLS().v_obj();
    const real*  w_obs  = field.meshValue(lv).valueObjLS().w_obj();

    const real* u     = field.meshValue(lv).valueNS().u();
    const real* v     = field.meshValue(lv).valueNS().v();
    const real* w     = field.meshValue(lv).valueNS().w();

    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real  kvis = air_property::KViscosity;

    const real c_ref  = parameters.c_ref_lbm();
    const real dx     = field.meshValue(lv).coordinates().dx();
    const real dt     = dx/c_ref;


#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::stream_collision_sgs_cpu (
        num_tasks_w_bc,
        id_tasks_w_bc,
        mesh_offsets3x3x3,
        f_lbm,
        fn_lbm,
        sgs_vis,
        T,Tn,
        u,v,w,
        lv_obs,
        u_obs, v_obs, w_obs,
        kvis,
        c_ref, dt
        );

    LBMCalculationCPU::stream_collision_sgs_cpu_wo_wall (
        num_tasks_wo_bc,
        id_tasks_wo_bc,
        mesh_offsets3x3x3,
        f_lbm,
        fn_lbm,
        sgs_vis,
        T,Tn,
        u,v,w,
        kvis,
        c_ref, dt
        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid_wo_bc(num_tasks_wo_bc, 1, 1);
    const dim3  grid_w_bc (num_tasks_w_bc , 1, 1);

    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
    const dim3  block(block3[0], block3[1], block3[2]);

    cudaFuncSetCacheConfig( LBMCalculationGPU::stream_collision_sgs_gpu,            cudaFuncCachePreferL1 );
    cudaFuncSetCacheConfig( LBMCalculationGPU::stream_collision_sgs_gpu_wo_wall,    cudaFuncCachePreferL1 );


    LBMCalculationGPU::stream_collision_sgs_gpu <<< grid_w_bc, block >>> (
        id_tasks_w_bc,
        mesh_offsets3x3x3,
        f_lbm,
        fn_lbm,
        sgs_vis,
        T, Tn,
        u,v,w,
        lv_obs,
        u_obs, v_obs, w_obs,
        kvis,
        c_ref, dt
        );

    LBMCalculationGPU::stream_collision_sgs_gpu_wo_wall <<< grid_wo_bc, block >>> (
        id_tasks_wo_bc,
        mesh_offsets3x3x3,
        f_lbm,
        fn_lbm,
        sgs_vis,
        T, Tn,
        u,v,w,
        kvis,
        c_ref, dt
        );

    cudaDeviceSynchronize();
#endif

}


void  LBMCalculation::
LBM_to_euler_variables(const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree& tree = field.tree();

    const real* f_lbm = field.meshValue(lv).valueLBM().f_lbm();
          real* rho   = field.meshValue(lv).valueNS().rho();
          real* u     = field.meshValue(lv).valueNS().u();
          real* v     = field.meshValue(lv).valueNS().v();
          real* w     = field.meshValue(lv).valueNS().w();

    const int*  mesh_offsets = tree.mesh_offsets();

#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::to_euler_variables_cpu (
        num_tasks,
        id_tasks,
        mesh_offsets,
        f_lbm,
        rho,
        u, v, w
        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid(num_tasks, 1, 1);
    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
    const dim3  block(block3[0], block3[1], block3[2]);

    cudaFuncSetCacheConfig( LBMCalculationGPU::to_euler_variables_gpu,  cudaFuncCachePreferL1 );
    LBMCalculationGPU::to_euler_variables_gpu <<< grid, block >>> (
        id_tasks,
        mesh_offsets,
        f_lbm,
        rho,
        u, v, w
        );
    cudaDeviceSynchronize();
#endif
}


void  LBMCalculation::
ScalarAdvection(const int lv, Field& field, const int*  id_tasks, const int  num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();


    // updates //
    const real*  scalar   = field.meshValue    (lv).valueNS().scalar();
          real*  scalar_n = field.meshValue_new(lv).valueNS().scalar();

    const real*  T   = field.meshValue    (lv).valueNS().T();
          real*  T_n = field.meshValue_new(lv).valueNS().T();

    // others //
    const real*  u  = field.meshValue(lv).valueNS().u();
    const real*  v  = field.meshValue(lv).valueNS().v();
    const real*  w  = field.meshValue(lv).valueNS().w();

    const real*  sgs_vis = field.meshValue_new(lv).valueNS().sgs_vis();

    const real*  lv_obs    = field.meshValue(lv).valueObjLS().lv_obj();
//    const real*  T_obs     = field.meshValue(lv).valueObjLS().T_obj();
    const real*  T_obs     = field.meshValues0()[lv].valueObjLS().T_obj();
    const real*  Tn_obs    = field.meshValues0()[lv].valueObjLS().Tn_obj();
    const int *  bcTypes_T = field.meshValue(lv).valueObjLS().bcTypes_T();

    const real*  sc_scalar = field.meshValue(lv).valueObjLS().sc_scalar();

    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real c_ref  = parameters.c_ref_lbm();
    const real dx     = field.meshValue(lv).coordinates().dx();
    const real dt     = dx/c_ref;

    const real coef_heatf = air_property::xi;
//    const real coef_heatf = air_property::xi / (dx*dx);

    const real time_second    = parameters.time_now();
//    const int  time_minutes   = int(time_second/60.0);
    const int  time_30minutes = int(time_second/60.0/30.0);
    const bool flag_source = (time_30minutes%4 == 0) ? true : false;
//    const bool flag_source = true;

    // boundary //
    const real time_minutes  = parameters.time_now()/(real)60.0;
    const int  _time_minutes = std::floor(time_minutes);
    const real time_weight0  = (real)1.0 - (time_minutes - _time_minutes);


#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::scalar_advection_cpu (
        num_tasks,
        id_tasks,
        mesh_offsets3x3x3,
        scalar,
        scalar_n,
        u,v,w,
        sgs_vis,
        lv_obs,
        sc_scalar,
        dx, dt
        );

    LBMCalculationCPU::temperature_advection_cpu(
        num_tasks,
        id_tasks,
        mesh_offsets3x3x3,
        T, T_n,
        u,v,w,
        sgs_vis,
        lv_obs,
        T_obs, Tn_obs, bcTypes_T, time_weight0,
        coef_heatf,
        dx, dt
        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid(num_tasks, 1, 1);
    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
    const dim3  block(block3[0], block3[1], block3[2]);

    cudaFuncSetCacheConfig( LBMCalculationGPU::scalar_advection_gpu,        cudaFuncCachePreferL1 );
    cudaFuncSetCacheConfig( LBMCalculationGPU::temperature_advection_gpu,   cudaFuncCachePreferL1 );

    if (flag_source) {
        LBMCalculationGPU::scalar_advection_w_source_gpu <<< grid, block >>> (
            id_tasks,
            mesh_offsets3x3x3,
            scalar,
            scalar_n,
            u,v,w,
            sgs_vis,
            lv_obs,
            sc_scalar,
            dx, dt
            );
    }
    else {
        LBMCalculationGPU::scalar_advection_gpu <<< grid, block >>> (
            id_tasks,
            mesh_offsets3x3x3,
            scalar,
            scalar_n,
            u,v,w,
            sgs_vis,
            lv_obs,
            dx, dt
            );
    }

    LBMCalculationGPU::temperature_advection_gpu <<< grid, block >>> (
        id_tasks,
        mesh_offsets3x3x3,
        T, T_n,
        u,v,w,
        sgs_vis,
        lv_obs,
        T_obs, Tn_obs, bcTypes_T, time_weight0,
        coef_heatf,
        dx, dt
        );

    cudaDeviceSynchronize();
#endif
}


void  LBMCalculation::
LBM_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const real* f_lbm  = field.meshValue_new(lv).  valueLBM().f_lbm();
          real* f_lbmF = field.meshValue_new(lv+1).valueLBM().f_lbm();

    const int* id_child          = tree.id_child();
    const int* mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real kvis  = air_property::KViscosity;
    const real c_ref = parameters.c_ref_lbm();
    const real dt0   = parameters.dt0();


#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::LBM_L2F_cpu(
        num_tasks,
        id_tasks,
        id_child,
        mesh_offsets3x3x3,
        lv,
        f_lbmF,
        f_lbm,
        kvis, c_ref, dt0
        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid(num_tasks, 8, 1);
    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
    const dim3  block(block3[0], block3[1], block3[2]);

    cudaFuncSetCacheConfig( LBMCalculationGPU::LBM_L2F_gpu, cudaFuncCachePreferL1 );
    LBMCalculationGPU::LBM_L2F_gpu <<< grid, block >>> (
        id_tasks,
        id_child,
        mesh_offsets3x3x3,
        lv,
        f_lbmF,
        f_lbm,
        kvis, c_ref, dt0
        );
    cudaDeviceSynchronize();
#endif
}


void  LBMCalculation::
LBM_L2FA(const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const real*  f_lbm  = field.meshValue    (lv).  valueLBM().f_lbm();
    const real*  fn_lbm = field.meshValue_new(lv).  valueLBM().f_lbm();
          real*  f_lbmF = field.meshValue_new(lv+1).valueLBM().f_lbm();

    const int*   id_child          = tree.id_child();
    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real kvis = air_property::KViscosity;
    const real c_ref  = parameters.c_ref_lbm();
    const real dt0    = parameters.dt0();


#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::LBM_L2FA_cpu(
        num_tasks,
        id_tasks,
        id_child,
        mesh_offsets3x3x3,
        lv,
        f_lbmF,
        f_lbm, fn_lbm,
        kvis, c_ref, dt0
        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid(num_tasks, 8, 1);
    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
    const dim3  block(block3[0], block3[1], block3[2]);

    cudaFuncSetCacheConfig( LBMCalculationGPU::LBM_L2FA_gpu, cudaFuncCachePreferL1 );
    LBMCalculationGPU::LBM_L2FA_gpu <<< grid, block >>> (
        id_tasks,
        id_child,
        mesh_offsets3x3x3,
        lv,
        f_lbmF,
        f_lbm, fn_lbm,
        kvis, c_ref, dt0
        );
    cudaDeviceSynchronize();
#endif
}


void  LBMCalculation::
Val_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int*   id_child          = tree.id_child();
    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

//    const real*  scalar  = field.meshValue_new(lv).  valueNS().scalar();
//          real*  scalarF = field.meshValue_new(lv+1).valueNS().scalar();

    std::vector<real*> Vals;
    std::vector<real*> ValFs;

    Vals .emplace_back( field.meshValue_new(lv).  valueNS().scalar() );
    ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().scalar() );

    Vals .emplace_back( field.meshValue_new(lv).  valueNS().T() );
    ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().T() );


    for (int ind=0; ind<Vals.size(); ind++) {
        const real* val  = Vals [ind];
              real* valF = ValFs[ind];

#ifdef CPU_CALCULATION__ // CPU : OMP //
        LBMCalculationCPU::Val_L2F_cpu(
            num_tasks,
            id_tasks,
            id_child,
            mesh_offsets3x3x3,
            lv,
            valF,
            val
//            scalarF,
//            scalar
            );
#endif

#ifdef GPU_CALCULATION__ // GPU //
        const dim3  grid(num_tasks, 8, 1);
        int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
        const dim3  block(block3[0], block3[1], block3[2]);

        cudaFuncSetCacheConfig( LBMCalculationGPU::Val_L2F_gpu, cudaFuncCachePreferL1 );
        LBMCalculationGPU::Val_L2F_gpu <<< grid, block >>> (
            id_tasks,
            id_child,
            mesh_offsets3x3x3,
            lv,
            valF,
            val
//            scalarF,
//            scalar
            );
        cudaDeviceSynchronize();
#endif
    }

}


void  LBMCalculation::
Val_L2FA (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int*   id_child          = tree.id_child();
    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

//    const real*  scalar  = field.meshValue    (lv).  valueNS().scalar();
//    const real*  scalarn = field.meshValue_new(lv).  valueNS().scalar();
//          real*  scalarF = field.meshValue_new(lv+1).valueNS().scalar();

    std::vector<real*> Vals;
    std::vector<real*> Valns;
    std::vector<real*> ValFs;

    Vals .emplace_back( field.meshValue    (lv).  valueNS().scalar() );
    Valns.emplace_back( field.meshValue_new(lv).  valueNS().scalar() );
    ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().scalar() );

    Vals .emplace_back( field.meshValue    (lv).  valueNS().T() );
    Valns.emplace_back( field.meshValue_new(lv).  valueNS().T() );
    ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().T() );

    for (int ind=0; ind<Vals.size(); ind++) {
        const real* val  = Vals [ind];
        const real* valn = Valns[ind];
              real* valF = ValFs[ind];


#ifdef CPU_CALCULATION__ // CPU : OMP //
        LBMCalculationCPU::Val_L2FA_cpu(
            num_tasks,
            id_tasks,
            id_child,
            mesh_offsets3x3x3,
            lv,
            valF,
            val,
            valn
//            scalarF,
//            scalar,
//            scalarn
            );
#endif

#ifdef GPU_CALCULATION__ // GPU //
        const dim3  grid(num_tasks, 8, 1);
        int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
        const dim3  block(block3[0], block3[1], block3[2]);

        cudaFuncSetCacheConfig( LBMCalculationGPU::Val_L2FA_gpu, cudaFuncCachePreferL1 );
        LBMCalculationGPU::Val_L2FA_gpu <<< grid, block >>> (
            id_tasks,
            id_child,
            mesh_offsets3x3x3,
            lv,
            valF,
            val,
            valn
//            scalarF,
//            scalar,
//            scalarn
            );
        cudaDeviceSynchronize();
#endif
    }

}


void  LBMCalculation::
LBM_F2L (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int    lvF    = lv+1;
    const real*  f_lbm  = field.meshValue_new(lvF).  valueLBM().f_lbm();
          real*  f_lbmC = field.meshValue_new(lvF-1).valueLBM().f_lbm();


    const int*  id_parent = tree.id_parent();
    const int*  mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const int*  node_pos_x = tree.node_pos_x();
    const int*  node_pos_y = tree.node_pos_y();
    const int*  node_pos_z = tree.node_pos_z();

    // F2L : lauch job at lv //
    //     : cal   job at lv+1 //
    const real kvis = air_property::KViscosity;
    const real c_ref  = parameters.c_ref_lbm();
    const real dt0    = parameters.dt0();


#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::LBM_F2L_cpu(
        num_tasks,
        id_tasks,
        id_parent,
        mesh_offsets3x3x3,
        lvF,
        node_pos_x, node_pos_y, node_pos_z,
        f_lbmC,
        f_lbm,
        kvis, c_ref, dt0
        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid(num_tasks, 1, 1);
    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF/2, DefAMR::NX_LEAF/2, DefAMR::NX_LEAF/2);
    const dim3  block(block3[0], block3[1], block3[2]);

    cudaFuncSetCacheConfig( LBMCalculationGPU::LBM_F2L_gpu, cudaFuncCachePreferL1 );
    LBMCalculationGPU::LBM_F2L_gpu <<< grid, block >>> (
        id_tasks,
        id_parent,
        mesh_offsets3x3x3,
        lvF,
        node_pos_x, node_pos_y, node_pos_z,
        f_lbmC,
        f_lbm,
        kvis, c_ref, dt0
        );
    cudaDeviceSynchronize();
#endif
}


void  LBMCalculation::
Val_F2L (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int    lvF    = lv+1;

    const int*  id_parent = tree.id_parent();
    const int*  mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const int*  node_pos_x = tree.node_pos_x();
    const int*  node_pos_y = tree.node_pos_y();
    const int*  node_pos_z = tree.node_pos_z();

    // F2L : lauch job at lv //
    //     : cal   job at lv+1 //


//    const real*  scalar  = field.meshValue_new(lvF).  valueNS().scalar();
//          real*  scalarC = field.meshValue_new(lvF-1).valueNS().scalar();

    std::vector<real*> Vals;
    std::vector<real*> ValCs;

    Vals .emplace_back( field.meshValue_new(lvF).  valueNS().scalar() );
    ValCs.emplace_back( field.meshValue_new(lvF-1).valueNS().scalar() );

    Vals .emplace_back( field.meshValue_new(lvF).  valueNS().T() );
    ValCs.emplace_back( field.meshValue_new(lvF-1).valueNS().T() );


    for (int ind=0; ind<Vals.size(); ind++) {
        const real* val  = Vals [ind];
              real* valC = ValCs[ind];

#ifdef CPU_CALCULATION__ // CPU : OMP //
        LBMCalculationCPU::Val_F2L_cpu(
            num_tasks,
            id_tasks,
            id_parent,
            mesh_offsets3x3x3,
            lvF,
            node_pos_x, node_pos_y, node_pos_z,
            valC,
            val
//            scalarC,
//            scalar
            );
#endif

#ifdef GPU_CALCULATION__ // GPU //
        const dim3  grid(num_tasks, 1, 1);
        int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF/2, DefAMR::NX_LEAF/2, DefAMR::NX_LEAF/2);
        const dim3  block(block3[0], block3[1], block3[2]);

        cudaFuncSetCacheConfig( LBMCalculationGPU::Val_F2L_gpu, cudaFuncCachePreferL1 );
        LBMCalculationGPU::Val_F2L_gpu <<< grid, block >>> (
            id_tasks,
            id_parent,
            mesh_offsets3x3x3,
            lvF,
            node_pos_x, node_pos_y, node_pos_z,
            valC,
            val
//            scalarC,
//            scalar
            );
        cudaDeviceSynchronize();
#endif
    }

}


void  LBMCalculation::
boundary_conditions(const int lv, Field& field, const int*  id_tasks, const int  num_tasks)
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();


    const real*  f_lbm  = field.meshValue    (lv).valueLBM().f_lbm();
          real*  fn_lbm = field.meshValue_new(lv).valueLBM().f_lbm();

    const real*  scalar   = field.meshValue    (lv).valueNS().scalar();
          real*  scalar_n = field.meshValue_new(lv).valueNS().scalar();

    const real*  rho = field.meshValue(lv).valueNS().rho();
    const real*  u   = field.meshValue(lv).valueNS().u();
    const real*  v   = field.meshValue(lv).valueNS().v();
    const real*  w   = field.meshValue(lv).valueNS().w();

          real*  Tn      = field.meshValue_new(lv).valueNS().T();
          real*  scalarn = field.meshValue_new(lv).valueNS().scalar();

    // solid wall //
    const real*  lv_obs      = field.meshValues0()[lv].valueObjLS().lv_obj();

    const real*  rho_obs     = field.meshValues0()[lv].valueObjLS().rho_obj();
    const real*  u_obs       = field.meshValues0()[lv].valueObjLS().u_obj();
    const real*  v_obs       = field.meshValues0()[lv].valueObjLS().v_obj();
    const real*  w_obs       = field.meshValues0()[lv].valueObjLS().w_obj();
    const real*  T_obs       = field.meshValues0()[lv].valueObjLS().T_obj();
    const real*  scalar_obs  = field.meshValues0()[lv].valueObjLS().scalar_obj();

    const real*  rhon_obs    = field.meshValues0()[lv].valueObjLS().rhon_obj();
    const real*  un_obs      = field.meshValues0()[lv].valueObjLS().un_obj();
    const real*  vn_obs      = field.meshValues0()[lv].valueObjLS().vn_obj();
    const real*  wn_obs      = field.meshValues0()[lv].valueObjLS().wn_obj();
    const real*  Tn_obs      = field.meshValues0()[lv].valueObjLS().Tn_obj();
    const real*  scalarn_obs = field.meshValues0()[lv].valueObjLS().scalarn_obj();

    // inflow & outflow bc //
    const int*   bcTypes_f = field.meshValue(lv).valueObjLS().bcTypes_f();

    const real*  dirichlet_weight = field.meshValue(lv).valueObjLS().dirichlet_weight();
    const real*  viscosity_weight = field.meshValue(lv).valueObjLS().viscosity_weight();

//    // solid wall //
//    const real*  lv_obs  = field.meshValue(lv).valueObjLS().lv_obj();
//    const real*  rho_obs = field.meshValue(lv).valueObjLS().rho_obj();
//    const real*  u_obs   = field.meshValue(lv).valueObjLS().u_obj();
//    const real*  v_obs   = field.meshValue(lv).valueObjLS().v_obj();
//    const real*  w_obs   = field.meshValue(lv).valueObjLS().w_obj();
//    const real*  T_obs   = field.meshValue(lv).valueObjLS().T_obj();
//
//    const real*  rhon_obs = field.meshValue(lv).valueObjLS().rhon_obj();
//    const real*  un_obs   = field.meshValue(lv).valueObjLS().un_obj();
//    const real*  vn_obs   = field.meshValue(lv).valueObjLS().vn_obj();
//    const real*  wn_obs   = field.meshValue(lv).valueObjLS().wn_obj();
//    const real*  Tn_obs   = field.meshValue(lv).valueObjLS().Tn_obj();
//
//    // inflow & outflow bc //
//    const int*   bcTypes_f = field.meshValue(lv).valueObjLS().bcTypes_f();
//
//    const real*  dirichlet_weight = field.meshValue(lv).valueObjLS().dirichlet_weight();
//    const real*  viscosity_weight = field.meshValue(lv).valueObjLS().viscosity_weight();


    const int*  mesh_offsets = tree.mesh_offsets();
    const int*  mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real c_ref  = parameters.c_ref_lbm();
    const real dx     = field.meshValue(lv).coordinates().dx();
    const real dt     = dx/c_ref;


//    const real time_second  = parameters.time_now();
//    const int  time_minutes = int(time_second/60.0);
//    const real time_weight0 = (real)1.0 - (time_second - time_minutes*(real)60.0)/(real)60.0;

    const real time_minutes  = parameters.time_now()/(real)60.0;
    const int  _time_minutes = std::floor(time_minutes);
    const real time_weight0  = (real)1.0 - (time_minutes - _time_minutes);

//    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    if (rank == 0) { std::cout << "time_second, time_minutes, time_weight0 = " << time_second << ", " << time_minutes << ", " << time_weight0 << std::endl; }


#ifdef CPU_CALCULATION__ // CPU : OMP //
    LBMCalculationCPU::boundary_conditions_in_wall_cpu (
        num_tasks,
        id_tasks,
        mesh_offsets,
        f_lbm, fn_lbm,
        rho,
        u, v, w,
        lv_obs,
        u_obs, v_obs, w_obs,
        c_ref, dx, dt
        );

//    LBMCalculationCPU::boundary_conditions_inflow_outflow_cpu (
//        num_tasks,
//        id_tasks,
//        mesh_offsets3x3x3,
//        f_lbm, fn_lbm,
//        rho,
//        u, v, w,
//        // solid wall //
//        lv_obs,
//        rho_obs,
//        u_obs, v_obs, w_obs,
//        // inflow & outflow bc //
//        bcTypes_f,
//        dirichlet_weight, viscosity_weight,
//        c_ref, dx, dt
//        );
#endif

#ifdef GPU_CALCULATION__ // GPU //
    const dim3  grid(num_tasks, 1, 1);
    int  block3[3]; set_cuda_block(block3, DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
    const dim3  block(block3[0], block3[1], block3[2]);

    LBMCalculationGPU::boundary_conditions_in_wall_gpu <<< grid, block >>> (
        id_tasks,
        mesh_offsets,
        f_lbm, fn_lbm,
        rho,
        u, v, w,
        lv_obs,
        u_obs, v_obs, w_obs,
        c_ref, dx, dt
        );

    cudaDeviceSynchronize();

//    LBMCalculationGPU::boundary_conditions_inflow_outflow_gpu <<< grid, block >>> (
//        id_tasks,
//        mesh_offsets3x3x3,
//        f_lbm, fn_lbm,
//        rho,
//        u, v, w,
//        // solid wall //
//        lv_obs,
//        rho_obs,
//        u_obs, v_obs, w_obs,
//        // inflow & outflow bc //
//        bcTypes_f,
//        dirichlet_weight, viscosity_weight,
//        c_ref, dx, dt
//        );

    LBMCalculationGPU::boundary_condition_data_assimilation_gpu <<< grid, block >>> (
        time_weight0,
        id_tasks,
        mesh_offsets3x3x3,
        f_lbm, fn_lbm,
        rho,
        u, v, w,
        Tn, scalarn,
        // solid wall //
        lv_obs,
        rho_obs,  u_obs,  v_obs,  w_obs,  T_obs,  scalar_obs,
        rhon_obs, un_obs, vn_obs, wn_obs, Tn_obs, scalarn_obs,
        // inflow & outflow bc //
        bcTypes_f,
        dirichlet_weight, viscosity_weight,
        c_ref, dx, dt
        );

    cudaDeviceSynchronize();
#endif

}


void  LBMCalculation::
ValuesCommunications(const int lv, Field& field)
{
    const int  n_leaf      = field.tree().number_of_nodes();
    const int  n_leaf_gmax = field.tree().number_of_nodes_global_max();

    bool  enable_comm = false;
    bool  tmp_enable_comm;


    // pack & unpack //
    tmp_enable_comm = mpiCommunication_.ValCommPackUnpack( lv, field, field.tree().mpiPackUnpackInfo(lv), field.tree(), field.meshValue_new(lv).valueNS(), field.meshValue_new(lv).valueObjLS(), field.meshValue_new(lv).valueBuff() );
    enable_comm = enable_comm || tmp_enable_comm;


    // pack & unpack //
    tmp_enable_comm = mpiCommunication_.LBMCommPackUnpack( lv, field, field.tree().lbmPutInfoSrcDst(lv), field.tree().mpiPackUnpackInfo(lv), n_leaf*DefAMR::NN_LEAF, n_leaf_gmax*DefAMR::NN_LEAF, field.tree(), field.meshValue_new(lv).valueLBM(), field.meshValue_new(lv).valueBuff() );
    enable_comm = enable_comm || tmp_enable_comm;


//    if (enable_comm) { MPI_Barrier(MPI_COMM_WORLD); }
}


void  LBMCalculation::
set_cuda_block(int block3[], int blockDimX, int blockDimY, int blockDimZ)
{
    int   nthread = CUDA_THREAD_MAX;
    block3[0] = std::min(blockDimX, nthread);  nthread = nthread/block3[0];
    block3[1] = std::min(blockDimY, nthread);  nthread = nthread/block3[1];
    block3[2] = std::min(blockDimZ, nthread);  nthread = nthread/block3[2];

    if (block3[0] == 0 || block3[1] == 0 || block3[2] == 0) {
        std::cout << __PRETTY_FUNCTION__ << " : error block3" << std::endl;
        exit(0);
    }

//    std::cout << block3[0] << ", " << block3[1] << ", " << block3[2] << std::endl;
//    exit(0);
}
