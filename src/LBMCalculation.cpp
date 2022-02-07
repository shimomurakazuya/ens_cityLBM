#include <climits>
#include "LBMCalculation.h"
#include "defineCal.h"
#include "defineAMR.h"
#include "defineFluidProperty.h"
#include "FuncScheduler.h"
#include "FuncLBMKernel.h"
#include "FuncLoop.h"
#include "FuncAllocate.h"
#include "FuncMath.h"

#include "Index.h"
#include "IndexLBM.h"
#include "BoundaryConditions.h"

#include "NSKernel.h"
#include "LBMKernel.h"
#include "AMRKernel.h"
#include "foreach.h"

#include "mpi_wrapper.hpp"
#include "range.hpp"

void  LBMCalculation::
LBM_data_assimilation(int t, Field& field)
{
    const int  time_minutes_start = field.optionParser().wrf_start();
#ifdef USE_WRF0
#warning USE_WRF0
    return;;
#endif

    const real time    = field.parameters().time_now();
    const int  time_minutes = int(time/60.0) + time_minutes_start;
    static int time_minutes_now = time_minutes;

    if (time_minutes != time_minutes_now) {
        if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": t, time, time_minutes, time_minutes_now = " << t << ", " << time << ", " << time_minutes << ", " << time_minutes_now << std::endl; }
        time_minutes_now = time_minutes;

        // particle filter
        #if defined(USE_PF_NUD5)
        constexpr int pf_time_minutes = 5;
        #elif defined(USE_PF_NUD15)
        constexpr int pf_time_minutes = 15;
        #else
        constexpr int pf_time_minutes = INT_MAX;
        #endif

        auto& nudgingCoefAdaptation = field.nudgingCoefAdaptation();
        if (time_minutes_now%pf_time_minutes == 0) { // default
            const float coef_ave_nudging_old = nudgingCoefAdaptation.coef_ave_nudging();
            const float coef_min_nudging_old = nudgingCoefAdaptation.coef_min_nudging();
            const float coef_max_nudging_old = nudgingCoefAdaptation.coef_max_nudging();
            nudgingCoefAdaptation.UpdateCoef(time_minutes_now, field.grids(), field.tree(), field.parameters(), field.meshValues0(), field.valueTimeAverage1min(), true); // meshValues0 ?????
            const float coef_ave_nudging = nudgingCoefAdaptation.coef_ave_nudging();
            const float coef_min_nudging = nudgingCoefAdaptation.coef_min_nudging();
            const float coef_max_nudging = nudgingCoefAdaptation.coef_max_nudging();

            if (comm_.is_rank0()) {
                std::cout << t << " : Particle filter (ave, min, max) (new <- old) : " << ", "
                    << coef_ave_nudging     << ", "  << coef_min_nudging     << ", "  << coef_max_nudging     << ", " << " <- " << ", "
                    << coef_ave_nudging_old << ", "  << coef_min_nudging_old << ", "  << coef_max_nudging_old << std::endl;
            }
        }

        // read: rho_obj, u_obj, ... etc. //
        BoundaryConditions  boundaryConditions(comm_);
        boundaryConditions.ReadWRFData   (time_minutes_now+1, field.grids(), field.tree(), field.parameters(), field.meshValues0(), nudgingCoefAdaptation.coef_nudging());
        boundaryConditions.ReadGroundData(time_minutes_now+1, field.grids(), field.tree(), field.parameters(), field.meshValues0());

        // swap: rho_obj <-> rhon_obj, .. etc //
        for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
            field.meshValues0()[lv].valueObjLS().swap_objs();
        }
        if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": load next data: fin" <<std::endl; }
    }
}


void  LBMCalculation::
LBM_incompressible_flow(int t, Field& field)
{
//    LBM_incompressible_flow_w_tb (t, field);
    LBM_incompressible_flow_wo_tb (t, field);
}



void  LBMCalculation::
LBM_incompressible_flow_w_tb(int t, Field& field)
{
    const int  lv_max = DefAMR::LV_MAX;

    // functions //
    std::function<void(int)>  func_cal_with_halo = [this, &field](int lv){
        field.swap_meshValue(lv);

        field.funcTimeInfo().StartTimer("LBM_to_euler_variables");
//        LBM_to_euler_variables(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        LBM_to_euler_variables_with_time_average(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("LBM_to_euler_variables");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_to_euler_variables", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


//        field.timerSimple().start_comm_timer(); // timer simple //
//        field.funcTimeInfo().StartTimer("ValuesCommunications");
//        TimeAverageCommunications(lv, field);
//        field.funcTimeInfo().StopTimer("ValuesCommunications");
//        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
//        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //

        // NS //
        field.funcTimeInfo().StartTimer("NSUpdate");
        NSUpdate(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("NSUpdate");
        field.funcTimeInfo().SubmitElapsedTimeInfo("NSUpdate", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());

        // MPI //
        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //

        // LBM //
        field.funcTimeInfo().StartTimer("LBM_stream_collision");
        LBM_stream_collision_opt(lv, field,
                                 field.taskID().id_task_array_with_halo_wo_bc(lv), field.taskID().num_task_array_with_halo_wo_bc(lv),
                                 field.taskID().id_task_array_with_halo_w_bc (lv), field.taskID().num_task_array_with_halo_w_bc (lv) );
        field.funcTimeInfo().StopTimer("LBM_stream_collision");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_stream_collision", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());

        // MPI //
        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        LBMCommunications(lv, field);
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


        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv+1, field);
        LBMCommunications(lv+1, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv+1, field.tree().mpiPutInfo(lv+1).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //
    };


    std::function<void(int)>  func_empty = [this, &field](int lv){};


    std::function<void(int)>  func_F2L = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("F2L");
        Val_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        LBM_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        field.funcTimeInfo().StopTimer("F2L");
        field.funcTimeInfo().SubmitElapsedTimeInfo("F2L", lv, field.taskID().id_tasks_F2L(lv).size(), field.parameters().step_now());

        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv, field);
        LBMCommunications(lv, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //
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
//        LBM_to_euler_variables(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        LBM_to_euler_variables_with_time_average(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("LBM_to_euler_variables");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_to_euler_variables", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());


//        field.timerSimple().start_comm_timer(); // timer simple //
//        field.funcTimeInfo().StartTimer("ValuesCommunications");
//        TimeAverageCommunications(lv, field);
//        field.funcTimeInfo().StopTimer("ValuesCommunications");
//        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
//        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //


        // NS
        field.funcTimeInfo().StartTimer("NSUpdate");
        NSUpdate(lv, field, field.taskID().id_task_array_with_halo(lv), field.taskID().num_task_array_with_halo(lv));
        field.funcTimeInfo().StopTimer("NSUpdate");
        field.funcTimeInfo().SubmitElapsedTimeInfo("NSUpdate", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());

        // MPI //
        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //

        // LBM //
        field.funcTimeInfo().StartTimer("LBM_stream_collision");
        LBM_stream_collision_opt(lv, field,
                                 field.taskID().id_task_array_with_halo_wo_bc(lv), field.taskID().num_task_array_with_halo_wo_bc(lv),
                                 field.taskID().id_task_array_with_halo_w_bc (lv), field.taskID().num_task_array_with_halo_w_bc (lv) );
        field.funcTimeInfo().StopTimer("LBM_stream_collision");
        field.funcTimeInfo().SubmitElapsedTimeInfo("LBM_stream_collision", lv, field.taskID().id_tasks_with_halo(lv).size(), field.parameters().step_now());

        // MPI //
        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        LBMCommunications(lv, field);
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


        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv+1, field);
        LBMCommunications(lv+1, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv+1, field.tree().mpiPutInfo(lv+1).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //
    };


    std::function<void(int)>  func_L2FA = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("L2FA");
        Val_L2FA(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        LBM_L2FA(lv, field, field.taskID().id_task_array_L2F(lv), field.taskID().num_task_array_L2F(lv));
        field.funcTimeInfo().StopTimer("F2LA");
        field.funcTimeInfo().SubmitElapsedTimeInfo("L2FA", lv, field.taskID().id_tasks_L2F(lv).size(), field.parameters().step_now());


        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv+1, field);
        LBMCommunications(lv+1, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv+1, field.tree().mpiPutInfo(lv+1).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //
    };


    std::function<void(int)>  func_F2L = [this, &field](int lv){
        field.funcTimeInfo().StartTimer("F2L");
        Val_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        LBM_F2L(lv, field, field.taskID().id_task_array_F2L(lv), field.taskID().num_task_array_F2L(lv));
        field.funcTimeInfo().StopTimer("F2L");
        field.funcTimeInfo().SubmitElapsedTimeInfo("F2L", lv, field.taskID().id_tasks_F2L(lv).size(), field.parameters().step_now());


        field.timerSimple().start_comm_timer(); // timer simple //
        field.funcTimeInfo().StartTimer("ValuesCommunications");
        NSCommunications(lv, field);
        LBMCommunications(lv, field);
        field.funcTimeInfo().StopTimer("ValuesCommunications");
        field.funcTimeInfo().SubmitElapsedTimeInfo("ValuesCommunications", lv, field.tree().mpiPutInfo(lv).size(), field.parameters().step_now());
        field.timerSimple().stop_comm_timer();  field.timerSimple().add_comm_timer(); // timer simple //
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
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();


    const real*  f_lbm  = field.meshValue    (lv).valueLBM().f_lbm();
          real*  fn_lbm = field.meshValue_new(lv).valueLBM().f_lbm();

          real*  sgs_vis = field.meshValue_new(lv).valueNS().sgs_vis();

    const real*  T      = field.meshValue    (lv).valueNS().T();
    const real*  Tn     = field.meshValue_new(lv).valueNS().T();

    const real*  lv_obj = field.meshValue(lv).valueObjLS().lv_obj();
    const real*  u_obj  = field.meshValue(lv).valueObjLS().u_obj();
    const real*  v_obj  = field.meshValue(lv).valueObjLS().v_obj();
    const real*  w_obj  = field.meshValue(lv).valueObjLS().w_obj();

    const real*  pad_obj = field.meshValue(lv).valueObjLS().pad_obj();

    const real* u     = field.meshValue(lv).valueNS().u();
    const real* v     = field.meshValue(lv).valueNS().v();
    const real* w     = field.meshValue(lv).valueNS().w();

    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real  kvis = air_property::KViscosity;

    const real c_ref  = parameters.c_ref_lbm();
    const real dx     = field.meshValue(lv).coordinates().dx();
    const real dt     = dx/c_ref;


    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks_w_bc,
        [=] __HD__ () {
            LBMKernel::stream_collision_sgs<NX_LEAF>(
                id_tasks_w_bc, num_tasks_w_bc, mesh_offsets3x3x3,
                f_lbm, fn_lbm,
                sgs_vis,
                T, Tn,
                u,v,w,
                lv_obj,
                u_obj, v_obj, w_obj,
                pad_obj,
                kvis,
                c_ref, dt
                );
        }
    );


    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks_wo_bc,
        [=] __HD__ () {
            LBMKernel::stream_collision_sgs_wo_wall<NX_LEAF>(
                id_tasks_wo_bc, num_tasks_wo_bc, mesh_offsets3x3x3,
                f_lbm, fn_lbm,
                sgs_vis,
                T, Tn,
                u,v,w,
                pad_obj,
                kvis,
                c_ref, dt
                );
        }
    );
}


void  LBMCalculation::
LBM_to_euler_variables(const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const Tree& tree = field.tree();

    const real* f_lbm = field.meshValue(lv).valueLBM().f_lbm();
          real* rho   = field.meshValue(lv).valueNS().rho();
          real* u     = field.meshValue(lv).valueNS().u();
          real* v     = field.meshValue(lv).valueNS().v();
          real* w     = field.meshValue(lv).valueNS().w();

    const int*  mesh_offsets = tree.mesh_offsets();

    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            LBMKernel::to_euler_variables<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets,
                f_lbm,
                rho,
                u, v, w
                );
        }
    );
}


void  LBMCalculation::
LBM_to_euler_variables_with_time_average(const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const Tree& tree = field.tree();

    const real* f_lbm = field.meshValue(lv).valueLBM().f_lbm();
          real* rho   = field.meshValue(lv).valueNS().rho();
          real* u     = field.meshValue(lv).valueNS().u();
          real* v     = field.meshValue(lv).valueNS().v();
          real* w     = field.meshValue(lv).valueNS().w();
    const real* T     = field.meshValue(lv).valueNS().T();

          real* u_mean_1min  = field.valueTimeAverage1min(lv).u_mean();
          real* v_mean_1min  = field.valueTimeAverage1min(lv).v_mean();
          real* w_mean_1min  = field.valueTimeAverage1min(lv).w_mean();
          real* T_mean_1min  = field.valueTimeAverage1min(lv).T_mean();

          real* vel2_fluc_1min  = field.valueTimeAverage1min(lv).vel2_fluc();
//          real* uu_fluc_1min  = field.valueTimeAverage1min(lv).uu_fluc();
//          real* vv_fluc_1min  = field.valueTimeAverage1min(lv).vv_fluc();
//          real* ww_fluc_1min  = field.valueTimeAverage1min(lv).ww_fluc();

    const real  time_1min =  field.valueTimeAverage1min(lv).time_min_av();

    const int*  mesh_offsets = tree.mesh_offsets();

    const real c_ref  = field.parameters().c_ref_lbm();
    const real dx     = field.meshValue(lv).coordinates().dx();
    const real dt     = dx/c_ref;

    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            LBMKernel::to_euler_variables_with_time_average<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets,
                f_lbm,
                rho, u, v, w, T,
                u_mean_1min,  v_mean_1min,  w_mean_1min,  T_mean_1min,
                vel2_fluc_1min,
//                uu_fluc_1min, vv_fluc_1min, ww_fluc_1min, 
                time_1min,
                dt
                );
        }
    );
}


void  LBMCalculation::
NSUpdate(const int lv, Field& field, const int*  id_tasks, const int  num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();


    // updates //
    const real*  scalar   = field.meshValue    (lv).valueNS().scalar();
          real*  scalar_n = field.meshValue_new(lv).valueNS().scalar();

    const real*  T   = field.meshValue    (lv).valueNS().T();
          real*  Tn = field.meshValue_new(lv).valueNS().T();

    // others //
    const real*  u  = field.meshValue(lv).valueNS().u();
    const real*  v  = field.meshValue(lv).valueNS().v();
    const real*  w  = field.meshValue(lv).valueNS().w();

    const real*  sgs_vis = field.meshValue(lv).valueNS().sgs_vis();

    const real*  lv_obj    = field.meshValue(lv).valueObjLS().lv_obj();
//    const real*  T_obj     = field.meshValue(lv).valueObjLS().T_obj();
    const real*  T_obj     = field.meshValues0()[lv].valueObjLS().T_obj();
    const real*  Tn_obj    = field.meshValues0()[lv].valueObjLS().Tn_obj();

    const real*  hflux_z_obj  = field.meshValues0()[lv].valueObjLS().hflux_z_obj();
    const real*  hfluxn_z_obj = field.meshValues0()[lv].valueObjLS().hfluxn_z_obj();

    const real*  sc_scalar = field.meshValue(lv).valueObjLS().sc_scalar();

    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
    const int    n_scalars = field.meshValue(lv).n_scalars();

    const real c_ref  = parameters.c_ref_lbm();
    const real dx     = field.meshValue(lv).coordinates().dx();
    const real dt     = dx/c_ref;

    const real coef_heatf = air_property::xi;
//    const real coef_heatf = air_property::xi / (dx*dx);

    const real time_second    = parameters.time_now();

    // boundary //
    const real time_minutes  = parameters.time_now()/(real)60.0;
    const int  _time_minutes = std::floor(time_minutes);
    const real time_weight0  = (real)1.0 - (time_minutes - _time_minutes);

    // avearage //
    const real* u_mean  = field.valueTimeAverage1min(lv).u_mean();
    const real* v_mean  = field.valueTimeAverage1min(lv).v_mean();
    const real* w_mean  = field.valueTimeAverage1min(lv).w_mean();
    const real* T_mean  = field.valueTimeAverage1min(lv).T_mean();

    #ifndef SC_SCALAR_ALWAYS
    /// JU2003, IOP6: July 16, 2003
    /// cont30min:: 9:00,11:00,13:00 CDT, 
    /// puff_once:: 15:00,15:20,15:40,16:00 CDT 
    ///// 7:00 CDT (420 min CDT) == 360 min @ WRF
    /////  --> cdt_minute = wrf_minute + 60
    const auto wrf_minute = time_minutes + field.optionParser().wrf_start();
    const auto cdt_minute = wrf_minute + 60;
    // flag source @ continuous release //
    const bool flag_source = false
        || ( 9*60    <= cdt_minute && cdt_minute <  9*60+30) // cont
        || (11*60    <= cdt_minute && cdt_minute < 11*60+30) //
        || (13*60    <= cdt_minute && cdt_minute < 13*60+30) //
    ; // flag_source
    //if(flag_source && comm_.is_rank0()) {
    //    std::cout << __PRETTY_FUNCTION__ << ": flag_source get true at:" 
    //        << parameters.step_now() << ", "
    //        << time_minutes << ", "
    //        << wrf_minute << ", "
    //        << int(cdt_minute)/60 << ":" << int(cdt_minute)%60
    //        << " (step, minute, minute@wrf, cdt)"
    //        << std::endl;
    //}
    // flag puff @ puff release //
    const bool flag_puff = false
        || (15*60    - dt/60./2. <= cdt_minute && cdt_minute < 15*60    + dt/60./2.) // puff
        || (15*60+20 - dt/60./2. <= cdt_minute && cdt_minute < 15*60+20 + dt/60./2.) // 
        || (15*60+40 - dt/60./2. <= cdt_minute && cdt_minute < 15*60+40 + dt/60./2.) //
        || (16*60    - dt/60./2. <= cdt_minute && cdt_minute < 16*60    + dt/60./2.) //
    ; // flag_puff
    if(flag_puff && comm_.is_rank0()) {
        std::cout << __PRETTY_FUNCTION__ << ": flag_puff get true at:" 
            << parameters.step_now() << ", "
            << time_minutes << ", "
            << wrf_minute << ", "
            << int(cdt_minute)/60 << ":" << int(cdt_minute)%60
            << " (step, minute, minute@wrf, cdt)"
            << std::endl;
    }
    #else
    constexpr bool flag_source = true, flag_puff = false;
    #endif // ifndef SC_SCALAR_ALWAYS

    // scalars //
    std::vector<real*> scalars;


    // puff release //
    const auto nn_max = field.meshValue(lv).nn_max();
    if(flag_puff) {
        const auto memType = field.meshValue(lv).valueNS().memType();

        real* scalar = field.meshValue(lv).valueNS().scalar();
        real* scalarn = field.meshValue_new(lv).valueNS().scalar();
        const real* sc_scalar = field.meshValue(lv).valueObjLS().sc_scalar();
        FuncAllocate::copy_values(scalar , sc_scalar, nn_max * n_scalars, memType);
        FuncAllocate::copy_values(scalarn, sc_scalar, nn_max * n_scalars, memType);
    }


    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
#if 0 // hflux : WRF
            NSKernel::TemperatureAdvDiff1st_hflux<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets3x3x3,
                Tn, T,
                u,v,w,
                sgs_vis,
                lv_obj,
                T_obj, Tn_obj,
                hflux_z_obj, hfluxn_z_obj,
                time_weight0,
                coef_heatf,
                dx, dt
                );
#else // hflux : LBM
            NSKernel::TemperatureAdvDiff1st_hflux_les<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets3x3x3,
                Tn, T,
                u,v,w,
                sgs_vis,
                lv_obj,
                T_obj, Tn_obj,
                time_weight0,
                coef_heatf,
                dx, dt,
                c_ref, 
                u_mean, v_mean, w_mean, T_mean
                );
#endif
        }
    );

    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            #if defined(SCALAR_WENO5)
            NSKernel::ScalarConvectionweno<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets3x3x3,
                nn_max, n_scalars,
                scalar_n,
                scalar,
                u,v,w,
                lv_obj,
                sc_scalar, flag_source,
                dx, dt
                );
            #else
            NSKernel::ScalarConvectionorg<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets3x3x3,
                nn_max, n_scalars,
                scalar_n,
                scalar,
                u,v,w,
                lv_obj,
                sc_scalar, flag_source,
                dx, dt
                );
            #endif
        }
    );

    #ifdef SCALAR_NON_PRIODIC
    const auto& coordinates = field.meshValue(lv).coordinates();
    //const auto nn_max = field.meshValue(lv).nn_max();
    const auto xmin = parameters.coefDomain().x_global_domain_min;
    const auto ymin = parameters.coefDomain().y_global_domain_min;
    const auto zmin = parameters.coefDomain().z_global_domain_min;
    const auto xmax = parameters.coefDomain().x_global_domain_max;
    const auto ymax = parameters.coefDomain().y_global_domain_max;
    const auto zmax = parameters.coefDomain().z_global_domain_max;

    foreach::exec_1d<foreach::opti>(
        nn_max,
        [=] __HD__ () {
            FOR_EACH1D_XX(i, nn_max) {
                const auto xc = coordinates.xnode(i);
                const auto yc = coordinates.ynode(i);
                //const auto zc = z[i];
                if(xc <= xmin + 2*dx || xmax - 2*dx <= xc
                || yc <= ymin + 2*dx || ymax - 2*dx <= yc) {
                    for(int ns=0; ns<n_scalars; ns++) {
                        scalar_n[i + ns*nn_max] = 0;
                    }
                }
            }
        }
    );

    #endif // SCALAR_NON_PRIODIC

}


void  LBMCalculation::
LBM_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const real* f_lbm  = field.meshValue_new(lv).  valueLBM().f_lbm();
          real* f_lbmF = field.meshValue_new(lv+1).valueLBM().f_lbm();

    const int* id_child          = tree.id_child();
    const int* mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const real kvis  = air_property::KViscosity;
    const real c_ref = parameters.c_ref_lbm();
    const real dt0   = parameters.dt0();

    foreach::exec_block_amr_L2F< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            AMRKernel::LBM_L2F<NX_LEAF>(
                id_tasks,
                num_tasks,
                id_child,
                mesh_offsets3x3x3,
                lv,
                f_lbmF,
                f_lbm,
                kvis, c_ref, dt0
                );
        }
    );
}


void  LBMCalculation::
LBM_L2FA(const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;
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

    foreach::exec_block_amr_L2F< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            AMRKernel::LBM_L2FA<NX_LEAF>(
                id_tasks,
                num_tasks,
                id_child,
                mesh_offsets3x3x3,
                lv,
                f_lbmF,
                f_lbm, fn_lbm,
                kvis, c_ref, dt0
                );
        }
    );
}


void  LBMCalculation::
Val_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int*   id_child          = tree.id_child();
    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
    const int    n_scalars = field.meshValue(lv).n_scalars();

    std::vector<real*> Vals;
    std::vector<real*> ValFs;

    for(int n=0; n<n_scalars; n++) {
        Vals .emplace_back( field.meshValue_new(lv).  valueNS().scalar(n) );
        ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().scalar(n) );
    }

    Vals .emplace_back( field.meshValue_new(lv).  valueNS().T() );
    ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().T() );


    for (int ind=0; ind<(int)Vals.size(); ind++) {
        const real* val  = Vals [ind];
              real* valF = ValFs[ind];

        foreach::exec_block_amr_L2F< foreach::opti, NX_LEAF >(
            num_tasks,
            [=] __HD__ () {
                AMRKernel::Val_L2F<NX_LEAF>(
                    id_tasks,
                    num_tasks,
                    id_child,
                    mesh_offsets3x3x3,
                    lv,
                    valF,
                    val
                    );
            }
        );
    }
}


void  LBMCalculation::
Val_L2FA (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int*   id_child          = tree.id_child();
    const int*   mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
    const int    n_scalars = field.meshValue(lv).n_scalars();

    std::vector<real*> Vals;
    std::vector<real*> Valns;
    std::vector<real*> ValFs;

    for(int n=0; n<n_scalars; n++) {
        Vals .emplace_back( field.meshValue    (lv).  valueNS().scalar(n) );
        Valns.emplace_back( field.meshValue_new(lv).  valueNS().scalar(n) );
        ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().scalar(n) );
    }

    Vals .emplace_back( field.meshValue    (lv).  valueNS().T() );
    Valns.emplace_back( field.meshValue_new(lv).  valueNS().T() );
    ValFs.emplace_back( field.meshValue_new(lv+1).valueNS().T() );


    for (int ind=0; ind<(int)Vals.size(); ind++) {
        const real* val  = Vals [ind];
        const real* valn = Valns[ind];
              real* valF = ValFs[ind];

            foreach::exec_block_amr_L2F< foreach::opti, NX_LEAF >(
                num_tasks,
                [=] __HD__ () {
                    AMRKernel::Val_L2FA<NX_LEAF>(
                        id_tasks,
                        num_tasks,
                        id_child,
                        mesh_offsets3x3x3,
                        lv,
                        valF,
                        val, valn
                        );
                }
            );
    }
}


void  LBMCalculation::
LBM_F2L (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

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


    foreach::exec_block_amr_F2L< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            AMRKernel::LBM_F2L<NX_LEAF>(
                id_tasks,
                num_tasks,
                id_parent,
                mesh_offsets3x3x3,
                lvF,
                node_pos_x, node_pos_y, node_pos_z,
                f_lbmC,
                f_lbm,
                kvis, c_ref, dt0
                );
        }
    );
}


void  LBMCalculation::
Val_F2L (const int lv, Field& field, const int*  id_tasks, const int num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    const int    lvF    = lv+1;

    const int*  id_parent = tree.id_parent();
    const int*  mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
    const int    nn_max = field.meshValue(lv).nn_max();
    const int    n_scalars = field.meshValue(lv).n_scalars();

    const int*  node_pos_x = tree.node_pos_x();
    const int*  node_pos_y = tree.node_pos_y();
    const int*  node_pos_z = tree.node_pos_z();

    // F2L : lauch job at lv //
    //     : cal   job at lv+1 //


    std::vector<real*> Vals;
    std::vector<real*> ValCs;

    for(int n=0; n<n_scalars; n++) {
        Vals .emplace_back( field.meshValue_new(lvF).  valueNS().scalar(n) );
        ValCs.emplace_back( field.meshValue_new(lvF-1).valueNS().scalar(n) );
    }

    Vals .emplace_back( field.meshValue_new(lvF).  valueNS().T() );
    ValCs.emplace_back( field.meshValue_new(lvF-1).valueNS().T() );

    for (int ind=0; ind<(int)Vals.size(); ind++) {
        const real* val  = Vals [ind];
              real* valC = ValCs[ind];

        foreach::exec_block_amr_F2L< foreach::opti, NX_LEAF >(
            num_tasks,
            [=] __HD__ () {
                AMRKernel::Val_F2L<NX_LEAF>(
                    id_tasks,
                    num_tasks,
                    id_parent,
                    mesh_offsets3x3x3,
                    lvF,
                    node_pos_x, node_pos_y, node_pos_z,
                    valC,
                    val
                    );
            }
        );
    }
}


void  LBMCalculation::
boundary_conditions(const int lv, Field& field, const int*  id_tasks, const int  num_tasks)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();


    const real*  f_lbm  = field.meshValue    (lv).valueLBM().f_lbm();
          real*  fn_lbm = field.meshValue_new(lv).valueLBM().f_lbm();

    const real*  scalar   = field.meshValue    (lv).valueNS().scalar();

    const real*  rho = field.meshValue(lv).valueNS().rho();
    const real*  u   = field.meshValue(lv).valueNS().u();
    const real*  v   = field.meshValue(lv).valueNS().v();
    const real*  w   = field.meshValue(lv).valueNS().w();

          real*  Tn      = field.meshValue_new(lv).valueNS().T();
          real*  scalarn = field.meshValue_new(lv).valueNS().scalar();

//    const real*  u_mean   = field.meshValue(lv).valueNS().u();
//    const real*  v_mean   = field.meshValue(lv).valueNS().v();
//    const real*  w_mean   = field.meshValue(lv).valueNS().w();
    const real*  u_mean   = field.valueTimeAverage1min(lv).u_mean();
    const real*  v_mean   = field.valueTimeAverage1min(lv).v_mean();
    const real*  w_mean   = field.valueTimeAverage1min(lv).w_mean();

    // solid wall //
    const real*  lv_obj      = field.meshValues0()[lv].valueObjLS().lv_obj();

    const real*  rho_obj     = field.meshValues0()[lv].valueObjLS().rho_obj();
    const real*  u_obj       = field.meshValues0()[lv].valueObjLS().u_obj();
    const real*  v_obj       = field.meshValues0()[lv].valueObjLS().v_obj();
    const real*  w_obj       = field.meshValues0()[lv].valueObjLS().w_obj();
    const real*  T_obj       = field.meshValues0()[lv].valueObjLS().T_obj();
    const real*  scalar_obj  = field.meshValues0()[lv].valueObjLS().scalar_obj();
//    const real*  hflux_z_obj = field.meshValues0()[lv].valueObjLS().hflux_z_obj();

    const real*  rhon_obj     = field.meshValues0()[lv].valueObjLS().rhon_obj();
    const real*  un_obj       = field.meshValues0()[lv].valueObjLS().un_obj();
    const real*  vn_obj       = field.meshValues0()[lv].valueObjLS().vn_obj();
    const real*  wn_obj       = field.meshValues0()[lv].valueObjLS().wn_obj();
    const real*  Tn_obj       = field.meshValues0()[lv].valueObjLS().Tn_obj();
    const real*  scalarn_obj  = field.meshValues0()[lv].valueObjLS().scalarn_obj();
//    const real*  hfluxn_z_obj = field.meshValues0()[lv].valueObjLS().hfluxn_z_obj();

    // inflow & outflow bc //
    const int*   bcTypes_f = field.meshValue(lv).valueObjLS().bcTypes_f();

//    const real*  dirichlet_weight = field.meshValue(lv).valueObjLS().dirichlet_weight(); // bug???
    const real*  dirichlet_weight = field.meshValues0()[lv].valueObjLS().dirichlet_weight();

//    // solid wall //
//    const real*  lv_obj  = field.meshValue(lv).valueObjLS().lv_obj();
//    const real*  rho_obj = field.meshValue(lv).valueObjLS().rho_obj();
//    const real*  u_obj   = field.meshValue(lv).valueObjLS().u_obj();
//    const real*  v_obj   = field.meshValue(lv).valueObjLS().v_obj();
//    const real*  w_obj   = field.meshValue(lv).valueObjLS().w_obj();
//    const real*  T_obj   = field.meshValue(lv).valueObjLS().T_obj();
//
//    const real*  rhon_obj = field.meshValue(lv).valueObjLS().rhon_obj();
//    const real*  un_obj   = field.meshValue(lv).valueObjLS().un_obj();
//    const real*  vn_obj   = field.meshValue(lv).valueObjLS().vn_obj();
//    const real*  wn_obj   = field.meshValue(lv).valueObjLS().wn_obj();
//    const real*  Tn_obj   = field.meshValue(lv).valueObjLS().Tn_obj();
//
//    // inflow & outflow bc //
//    const int*   bcTypes_f = field.meshValue(lv).valueObjLS().bcTypes_f();
//
//    const real*  dirichlet_weight = field.meshValue(lv).valueObjLS().dirichlet_weight();


    const int*  mesh_offsets = tree.mesh_offsets();
//    const int*  mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
    const int   nn_max    = field.meshValue(lv).nn_max();
    const int   n_scalars = field.meshValue(lv).n_scalars();

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


    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            LBMKernel::boundary_conditions_in_wall<NX_LEAF>(
                id_tasks, num_tasks, mesh_offsets,
                fn_lbm,
                lv_obj,
                u_obj, v_obj, w_obj,
                c_ref, dx, dt
                );
        }
    );

    foreach::exec_block_amr< foreach::opti, NX_LEAF >(
        num_tasks,
        [=] __HD__ () {
            LBMKernel::boundary_conditions_data_assimilation<NX_LEAF>(
                nn_max,
                n_scalars,
                time_weight0,
                id_tasks,
                num_tasks,
                mesh_offsets,
                f_lbm, fn_lbm,
                u_mean, v_mean, w_mean,
                Tn, scalarn,
                // solid wall //
                lv_obj,
                rho_obj,  u_obj,  v_obj,  w_obj,  T_obj,  scalar_obj,
                rhon_obj, un_obj, vn_obj, wn_obj, Tn_obj, scalarn_obj,
                // inflow & outflow bc //
                bcTypes_f,
                dirichlet_weight,
                c_ref, dx, dt
                );
        }
    );

}


void  LBMCalculation::
NSCommunications(const int lv, Field& field)
{
    TCommunications(lv, field);
    ScalarCommunications(lv, field);
}


void  LBMCalculation::
TCommunications(const int lv, Field& field)
{
#ifdef USE_TEMPORAL_BLOCKING
    const int  n_leaf      = field.tree().number_of_nodes();
    const int  n_leaf_gmax = field.tree().number_of_nodes_global_max();

    bool  enable_comm = false;
    bool  tmp_enable_comm;


    tmp_enable_comm = mpiCommunication_.TCommPackUnpack( lv, field, field.tree().mpiPackUnpackInfo(lv), field.tree(), field.meshValue_new(lv).valueNS(), field.meshValue_new(lv).valueObjLS(), field.valueBuff() );
    enable_comm = enable_comm || tmp_enable_comm;


//    if (enable_comm) { MPI_Barrier(MPI_COMM_WORLD); }
#else
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const MPIPackUnpackInfo&  mpiPackUnpackInfo = field.tree().mpiPackUnpackInfo(lv);
          ValueBuff&          valueBuff         = field.valueBuff();
    real* _T = field.meshValue_new(lv).valueNS().T();

    mpiCommunicationMG_.ValCommPackUnpackMG<NX_LEAF>(
        mpiPackUnpackInfo.sendPackUnpackCommInfo, 
        mpiPackUnpackInfo.recvPackUnpackCommInfo, 
//        mpiPackUnpackInfo.all, 
//        mpiPackUnpackInfo.slice3, 
        mpiPackUnpackInfo.slice1, 
        std::vector<real*>{ _T  }, valueBuff
        );
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}


void  LBMCalculation::
ScalarCommunications(const int lv, Field& field)
{
#ifdef USE_TEMPORAL_BLOCKING
    const int  n_leaf      = field.tree().number_of_nodes();
    const int  n_leaf_gmax = field.tree().number_of_nodes_global_max();

    bool  enable_comm = false;
    bool  tmp_enable_comm;


    tmp_enable_comm = mpiCommunication_.ScalarCommPackUnpack( lv, field, field.tree().mpiPackUnpackInfo(lv), field.tree(), field.meshValue_new(lv).valueNS(), field.meshValue_new(lv).valueObjLS(), field.valueBuff() );
    enable_comm = enable_comm || tmp_enable_comm;


//    if (enable_comm) { MPI_Barrier(MPI_COMM_WORLD); }
#else
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const MPIPackUnpackInfo&  mpiPackUnpackInfo = field.tree().mpiPackUnpackInfo(lv);
          ValueBuff&          valueBuff         = field.valueBuff();

    const int nn_max = field.meshValue(lv).nn_max();
    const int n_scalars = field.meshValue(lv).n_scalars();

    std::vector<real*>  vals;
    for(int n=0; n<n_scalars; n++) {
        int idx = n * nn_max;
        vals.push_back( &(field.meshValue_new(lv).valueNS().scalar()[idx]) );
    }

    #if defined(SCALAR_WENO5)
    mpiCommunicationMG_.ValCommPackUnpackMG<NX_LEAF>(
        mpiPackUnpackInfo.sendPackUnpackCommInfo, 
        mpiPackUnpackInfo.recvPackUnpackCommInfo, 
        mpiPackUnpackInfo.slice3, 
        vals, valueBuff
        );
    #else
    mpiCommunicationMG_.ValCommPackUnpackMG<NX_LEAF>(
        mpiPackUnpackInfo.sendPackUnpackCommInfo, 
        mpiPackUnpackInfo.recvPackUnpackCommInfo, 
        mpiPackUnpackInfo.slice1, 
        vals, valueBuff
        );
    #endif
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}


void  LBMCalculation::
TimeAverageCommunications(const int lv, Field& field)
{
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const MPIPackUnpackInfo&  mpiPackUnpackInfo = field.tree().mpiPackUnpackInfo(lv);
          ValueBuff&          valueBuff         = field.valueBuff();
    real* _u_mean_1min  = field.valueTimeAverage1min(lv).u_mean();
    real* _v_mean_1min  = field.valueTimeAverage1min(lv).v_mean();
    real* _w_mean_1min  = field.valueTimeAverage1min(lv).w_mean();
    real* _T_mean_1min  = field.valueTimeAverage1min(lv).T_mean();

    real* _vel2_fluc_1min  = field.valueTimeAverage1min(lv).vel2_fluc();
//    real* _uu_fluc_1min  = field.valueTimeAverage1min(lv).uu_fluc();
//    real* _vv_fluc_1min  = field.valueTimeAverage1min(lv).vv_fluc();
//    real* _ww_fluc_1min  = field.valueTimeAverage1min(lv).ww_fluc();

    mpiCommunicationMG_.ValCommPackUnpackMG<NX_LEAF>(
        mpiPackUnpackInfo.sendPackUnpackCommInfo, 
        mpiPackUnpackInfo.recvPackUnpackCommInfo, 
//        mpiPackUnpackInfo.all, 
//        mpiPackUnpackInfo.slice3, 
        mpiPackUnpackInfo.slice1, 
        std::vector<real*>{
            _u_mean_1min,  _v_mean_1min,  _w_mean_1min,  _T_mean_1min,
            _vel2_fluc_1min
//            _uu_fluc_1min, _vv_fluc_1min, _ww_fluc_1min
        },
        valueBuff
        );
    MPI_Barrier(MPI_COMM_WORLD);
}


void  LBMCalculation::
LBMCommunications(const int lv, Field& field)
{
#ifdef USE_TEMPORAL_BLOCKING
    const int  n_leaf      = field.tree().number_of_nodes();
    const int  n_leaf_gmax = field.tree().number_of_nodes_global_max();

    bool  enable_comm = false;
    bool  tmp_enable_comm;


    tmp_enable_comm = mpiCommunication_.LBMCommPackUnpack( lv, field, field.tree().lbmPutInfoSrcDst(lv), field.tree().mpiPackUnpackInfo(lv), n_leaf*DefAMR::NN_LEAF, n_leaf_gmax*DefAMR::NN_LEAF, field.tree(), field.meshValue_new(lv).valueNS(), field.meshValue_new(lv).valueLBM(), field.valueBuff() );
    enable_comm = enable_comm || tmp_enable_comm;
#else
    constexpr int NX_LEAF = DefAMR::NX_LEAF;

    const MPIPackUnpackInfo&  mpiPackUnpackInfo = field.tree().mpiPackUnpackInfo(lv);
          ValueBuff&          valueBuff         = field.valueBuff();
    real* _lbm = field.meshValue_new(lv).valueLBM().f_lbm();

    mpiCommunicationMG_.ValCommPackUnpackMG<NX_LEAF>(
        mpiPackUnpackInfo.sendPackUnpackCommInfo, 
        mpiPackUnpackInfo.recvPackUnpackCommInfo, 
//        mpiPackUnpackInfo.all_lbm, 
//        mpiPackUnpackInfo.slice3_lbm, 
        mpiPackUnpackInfo.slice1_lbm, 
        std::vector<real*>{ _lbm }, valueBuff
        );
    MPI_Barrier(MPI_COMM_WORLD);
#endif
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
