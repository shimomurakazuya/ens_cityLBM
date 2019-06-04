#include "OutputFunc.h"
#include <mpi.h>
#include <cmath>
#include "defineFluidProperty.h"
#include "FuncMath.h"
#include "FuncObj.h"
#include "PostprocessChannelFlow.h"
#include "PostprocessTaylorGreen.h"
#include "PostprocessNaturalConvection2d.h"
#include "PostprocessNaturalConvection3d.h"


void OutputFunc::
OutputFluidData(int t, Field& field, WorkerThread& iothread)
{
#if 1
    field.parameter_time(t);

    const Parameters& parameters = field.parameters();
    if ( parameters.coefCFROutput().is_cout_step(t) || parameters.coefCFROutput().is_fout_step(t) ) {
        field.update_meshValues(field.grids(), field.tree());
    }

    // output //
    cout_step(t, parameters);
    cout_monitor(t, field);
    cout_timer(t, field);


//    // channel flow //
//    output_channel_flow_data(t, field);

//    // taylor green //
//    output_taylor_green_data(t, field);

//    // natural convection //
//    output_natural_convection2d_data(t, field);
//    output_natural_convection3d_data(t, field);


    // hdf, vtu //
    output_data(t, field, iothread);

    // monitor //
    output_monitor_data(t, field);


    // elapsed time //
    field.allTimeInfo(). OutputElapsedTimeInfo(t);
    field.funcTimeInfo().OutputElapsedTimeInfo(t);
    field.mpiTimeInfo(). OutputElapsedTimeInfo(t);
#endif
}


void OutputFunc::
cout_step(int t, const Parameters& parameters)
const
{
    if ( !parameters.coefCFROutput().is_cout_step(t) ) { return; }
    if ( !rank_ == 0 ) { return; }

    std::cout << std::endl;
    std::cout << "t / step_end = " << t << " / " << parameters.step_end() << " ( time = " << parameters.time_now() << " / " << parameters.time_end() << " ) " << std::endl;
//    if ( t%20 == 0 ) { std::cout << std::flush; }
}


void OutputFunc::
cout_monitor(int t, const Field& field)
const
{
    const Tree&       tree       = field.tree();
    const Parameters& parameters = field.parameters();

    if ( !parameters.coefCFROutput().is_cout_step(t) ) { return; }

    const int  n_leaf = tree.number_of_nodes();
    std::vector<int>  id_tasks;
    for (int i=0; i<n_leaf; i++) {
        if ( tree.nodes(i)->nodeCalFlags().Cal() ) {
            id_tasks.push_back(i);
        }
    }

    const double  kvis = air_property::KViscosity;


    const double dx_coarse = field.meshValue(0).coordinates().dx();
    const double c_ref     = parameters.c_ref_lbm();
    const double dt_coarse = dx_coarse/c_ref;

    double  num = 0.0;
    double  energy  = 0.0;
    double  rho_all = 0.0;
    double  vel_max = 0.0;
    double  scalar_all = 0.0;
    double  scalar_min =  1.0e8;
    double  scalar_max = -1.0e8;
    double  T_average = 0.0;
    double  T_min =  1.0e8;
    double  T_max = -1.0e8;
    double  sgs_vis_ave = 0.0;
    double  sgs_vis_max = 0.0;

#pragma omp parallel for reduction(+:num,energy,rho_all,scalar_all,T_average,sgs_vis_ave) reduction(min:scalar_min,T_min) reduction(max:vel_max,scalar_max,T_max,sgs_vis_max)
    for (int it=0; it<(int)id_tasks.size(); it++) {
        const int  i = id_tasks[it];

        const int  lv     = tree.nodes(i)->level();

        const double dx     = field.meshValue(lv).coordinates().dx();
//        const double c_ref  = parameters.c_ref_lbm();

        const real*  rho = field.meshValue(lv).valueNS().rho();
        const real*  u   = field.meshValue(lv).valueNS().u();
        const real*  v   = field.meshValue(lv).valueNS().v();
        const real*  w   = field.meshValue(lv).valueNS().w();
        const real*  T   = field.meshValue(lv).valueNS().T();

        const real*  sgs_vis = field.meshValue(lv).valueNS().sgs_vis();

        const real*  scalar = field.meshValue(lv).valueNS().scalar();

        const real*  lv_obs = field.meshValue(lv).valueObjLS().lv_obj();
//        const real*  u_obs  = field.meshValue(lv).valueObjLS().u_obj();
//        const real*  v_obs  = field.meshValue(lv).valueObjLS().v_obj();
//        const real*  w_obs  = field.meshValue(lv).valueObjLS().w_obj();
        const Array3D<int>  offsets = tree.nodes(i)->neighbor_mesh_offsets();

        // monitor //
        const int  offset0 = offsets.offset0();

        for (int j=0; j<DefAMR::NN_LEAF; j++) {
            const int  id = j + offset0;

            if ( FuncObj::is_obj(lv_obs[id]) ) { continue; }

            num        += 1.0;
            vel_max     = std::fmax( (float)vel_max, (float)sqrt(u[id]*u[id] + v[id]*v[id] + w[id]*w[id]) );
            rho_all    +=     pow(dx, 3) * rho[id];

            scalar_all +=     pow(dx, 3) * scalar[id];
            scalar_min = std::fmin((float)scalar_min, (float)scalar[id]);
            scalar_max = std::fmax((float)scalar_max, (float)scalar[id]);

            energy  += 0.5*pow(dx, 3) * FuncMath::norm(u[id], v[id], w[id]) * pow(c_ref, 2) * rho[id];

            sgs_vis_ave += sgs_vis[id];
            sgs_vis_max  = std::fmax( (float)sgs_vis_max, (float)sgs_vis[id] );

            T_average += T[id];
            T_min = std::fmin((float)T_min, T[id]);
            T_max = std::fmax((float)T_max, T[id]);
        }
    }
//    std::cout << "scalar_all = " << scalar_all << std::endl;

    double num_g     = 0.0;
    double vel_max_g = 0.0;
    double rho_all_g = 0.0;
    double energy_g  = 0.0;
    double scalar_all_g = 0.0;
    double scalar_min_g = 0.0;
    double scalar_max_g = 0.0;
    double T_average_g  = 0.0;
    double T_min_g  = 0.0;
    double T_max_g  = 0.0;
    double sgs_vis_ave_g  = 0.0;
    double sgs_vis_max_g  = 0.0;

    MPI_Request requ[12];

    MPI_Iallreduce(&num,          &num_g,          1, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, &requ[0]);
    MPI_Iallreduce(&vel_max,      &vel_max_g,      1, MPI_DOUBLE,  MPI_MAX, MPI_COMM_WORLD, &requ[1]);
    MPI_Iallreduce(&rho_all,      &rho_all_g,      1, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, &requ[2]);
    MPI_Iallreduce(&energy,       &energy_g,       1, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, &requ[3]);
    MPI_Iallreduce(&scalar_all,   &scalar_all_g,   1, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, &requ[4]);
    MPI_Iallreduce(&scalar_min,   &scalar_min_g,   1, MPI_DOUBLE,  MPI_MIN, MPI_COMM_WORLD, &requ[5]);
    MPI_Iallreduce(&scalar_max,   &scalar_max_g,   1, MPI_DOUBLE,  MPI_MAX, MPI_COMM_WORLD, &requ[6]);
    MPI_Iallreduce(&T_average,    &T_average_g,    1, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, &requ[7]);
    MPI_Iallreduce(&T_min,        &T_min_g,        1, MPI_DOUBLE,  MPI_MIN, MPI_COMM_WORLD, &requ[8]);
    MPI_Iallreduce(&T_max,        &T_max_g,        1, MPI_DOUBLE,  MPI_MAX, MPI_COMM_WORLD, &requ[9]);
    MPI_Iallreduce(&sgs_vis_ave,  &sgs_vis_ave_g,  1, MPI_DOUBLE,  MPI_SUM, MPI_COMM_WORLD, &requ[10]);
    MPI_Iallreduce(&sgs_vis_max,  &sgs_vis_max_g,  1, MPI_DOUBLE,  MPI_MAX, MPI_COMM_WORLD, &requ[11]);

    MPI_Waitall(12, requ, MPI_STATUSES_IGNORE);

    T_average_g   /= num_g;
    sgs_vis_ave   /= num_g;
    sgs_vis_ave_g /= num_g;

    // number of leaves //
    const int  n_leaf_g     = tree.number_of_nodes_global_total();
    const int  n_leaf_g_max = tree.number_of_nodes_global_max();
    const int  n_leaf_g_min = tree.number_of_nodes_global_min();

    int  n_leaf_lv_g    [DefAMR::LV_MAX];
    int  n_leaf_lv_g_max[DefAMR::LV_MAX];
    int  n_leaf_lv_g_min[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        n_leaf_lv_g    [i] = tree.number_of_nodes_lv_global_total(i);
        n_leaf_lv_g_max[i] = tree.number_of_nodes_lv_global_max  (i);
        n_leaf_lv_g_min[i] = tree.number_of_nodes_lv_global_min  (i);
    }

    // number of time-integrated leaves //
    int  n_upleaf[DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        n_upleaf[i] = field.taskID().num_task_array_lbm(i);
    }

    int  n_upleaf_g_sum[DefAMR::LV_MAX];
    int  n_upleaf_g_max[DefAMR::LV_MAX];
    int  n_upleaf_g_min[DefAMR::LV_MAX];

    MPI_Request requ_leaf[3*DefAMR::LV_MAX];
    for (int i=0; i<DefAMR::LV_MAX; i++) {
        MPI_Iallreduce(&n_upleaf[i],    &n_upleaf_g_sum[i],  1, MPI_INT,  MPI_SUM, MPI_COMM_WORLD, &requ_leaf[0 + 3*i]);
        MPI_Iallreduce(&n_upleaf[i],    &n_upleaf_g_max[i],  1, MPI_INT,  MPI_MAX, MPI_COMM_WORLD, &requ_leaf[1 + 3*i]);
        MPI_Iallreduce(&n_upleaf[i],    &n_upleaf_g_min[i],  1, MPI_INT,  MPI_MIN, MPI_COMM_WORLD, &requ_leaf[2 + 3*i]);
    }
    MPI_Waitall(3*DefAMR::LV_MAX,  requ_leaf,  MPI_STATUSES_IGNORE);


    MPI_Barrier(MPI_COMM_WORLD);
    if (rank_ == 0) {
        // allocate //
        std::cout << "n_leaf(total)       = " << n_leaf_g << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_leaf_lv_g[i] << ", "; }
        std::cout << " ) ";
        std::cout << " : nn_global (10^6) = " << n_leaf_g * DefAMR::NN_LEAF / 1.0e6 << std::endl;

        std::cout << "n_leaf(max@node)    = " << n_leaf_g_max << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_leaf_lv_g_max[i] << ", "; }
        std::cout << " ) ";
        std::cout << " : nn_max    (10^6) = " << n_leaf_g_max * DefAMR::NN_LEAF / 1.0e6 << std::endl;

        std::cout << "n_leaf(min@node)    = " << n_leaf_g_min << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_leaf_lv_g_min[i] << ", "; }
        std::cout << " ) ";
        std::cout << " : nn_min    (10^6) = " << n_leaf_g_min * DefAMR::NN_LEAF / 1.0e6 << std::endl;

        // updated //
        std::cout << "n_leaf(total, lv)    : updated           = "<< " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_sum[i] << ", "; }
        std::cout << " ) ";
        std::cout << "grid points(total, lv) 10^6 : updated    = "<< " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_sum[i] * DefAMR::NN_LEAF/1.0e6 << ", "; }
        std::cout << " ) " << std::endl;

        std::cout << "n_leaf(max@node, lv) : updated           = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_max[i] << ", "; }
        std::cout << " ) ";
        std::cout << "grid points(max@node, lv) 10^6 : updated = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_max[i] * DefAMR::NN_LEAF/1.0e6 << ", "; }
        std::cout << " ) " << std::endl;

        std::cout << "n_leaf(min@node, lv) : updated           = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_min[i] << ", "; }
        std::cout << " ) ";
        std::cout << "grid points(min@node, lv) 10^6 : updated = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_min[i] * DefAMR::NN_LEAF/1.0e6 << ", "; }
        std::cout << " ) " << std::endl;


        std::cout << "vel_max             = " << vel_max_g    << " :  ( " << vel_max_g*parameters.c_ref_lbm() << " m/s )"<< std::endl;
        std::cout << "rho_all             = " << rho_all_g    << std::endl;
        std::cout << "energy              = " << energy_g     << std::endl;
        std::cout << "scalar(all,min,max) = " << scalar_all_g << " ( " << scalar_min_g << ", " << scalar_max_g << " ) " << std::endl;
        std::cout << "T(ave,min,max)      = " << T_average_g  << " ( " << T_min_g << ", " << T_max_g << " ) " << std::endl;
        std::cout << "sgs_vis/vis         = " << sgs_vis_ave_g/kvis << ", " << sgs_vis_max_g/kvis << " (average, max)"  << std::endl;

        std::cout << "dt(coarse)          = " << dt_coarse << " (sec.)" << std::endl;
        std::cout << "dx                  = ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << field.meshValue(i).coordinates().dx() << ", "; }
        std::cout << std::endl;
    }
}


void OutputFunc::
cout_timer(int t, Field& field)
{
    static int  is_first_cal;
    static int  is_first_mlups;
    const Parameters& parameters = field.parameters();

    if ( !parameters.coefCFROutput().is_cout_step(t) ) { return; }

    if (is_first_cal == 0) {
    }
    else {
        gettimeofday(&t_end_,   NULL);

        const int     cout_step      = parameters.coefCFROutput().cout_step;
        const double  mupdated_grids = field.taskID().num_task_lbm_per_step() * DefAMR::NN_LEAF/1.0e6 * cout_step;
        const double  elapsed_time_cout = cal_elapsed_time_msec(t_mid_, t_end_) / 1000.0;
        const double  mlups = mupdated_grids / elapsed_time_cout;

        double cal_time  = field.timerSimple().t_cal_total();
        double comm_time = field.timerSimple().t_comm_total();

        MPI_Request requ[8];
        double elapsed_time_cout_g[3]; // average, min, max //
        double mlups_g[3];

        double cal_time_g;
        double comm_time_g;

        MPI_Iallreduce(&elapsed_time_cout, &elapsed_time_cout_g[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &requ[0]);
        MPI_Iallreduce(&elapsed_time_cout, &elapsed_time_cout_g[1], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &requ[1]);
        MPI_Iallreduce(&elapsed_time_cout, &elapsed_time_cout_g[2], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &requ[2]);

        MPI_Iallreduce(&mlups, &mlups_g[0], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &requ[3]);
        MPI_Iallreduce(&mlups, &mlups_g[1], 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &requ[4]);
        MPI_Iallreduce(&mlups, &mlups_g[2], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &requ[5]);

        MPI_Iallreduce(&cal_time , &cal_time_g , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &requ[6]);
        MPI_Iallreduce(&comm_time, &comm_time_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &requ[7]);

        MPI_Waitall(8, requ, MPI_STATUSES_IGNORE);

        int ncpu; MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
        elapsed_time_cout_g[0] /= ncpu;
        mlups_g[0]             /= ncpu;

        cal_time  /= ncpu;
        comm_time /= ncpu;


        if (rank_ == 0) {
            std::cout << "*****" << std::endl;
            std::cout << "MLUPS w/o halo,  total, mean = " << mlups_g[0]*ncpu << ", " << mlups_g[0] << " :  (min, max) = ( " << mlups_g[1] << ", " << mlups_g[2] << " ) " << std::endl;
            std::cout << "sec/step = " << elapsed_time_cout_g[2]/cout_step << "  ( total elapsed_time (minites) = " << (double)cal_elapsed_time_msec(t_begin_, t_end_)/1000.0/60.0 << " ) " << std::endl;
            std::cout << "time  ratio (gpu, comm) = " << (cal_time - comm_time)/cal_time << ", " << comm_time/cal_time << std::endl;
            std::cout << "*****" << std::endl;

            std::ofstream fout_mlups;
            if (is_first_mlups == 0)  { fout_mlups.open("../io/output/MLUPS.txt"); }
            else                      { fout_mlups.open("../io/output/MLUPS.txt", std::ios_base::app); }

            if (is_first_mlups == 0) { fout_mlups    << "MLUPS_total, MLUPS_per_GPU, total_time_min, cal_total, gpu, comm" << std::endl; }
            fout_mlups  << mlups_g[0]*ncpu  << "," << mlups_g[0] << "," << (double)cal_elapsed_time_msec(t_begin_, t_end_)/1000.0/60.0 << ", " << cal_time_g << "," << cal_time_g - comm_time_g << ", " << comm_time_g << std::endl;

            fout_mlups.close();

            is_first_mlups = 1;
        }
        field.timerSimple().reset_cal_timer();
        field.timerSimple().reset_comm_timer();

        gettimeofday(&t_mid_,   NULL);
    }
    is_first_cal = 1;
}


void OutputFunc::
output_data(int t, Field& field, WorkerThread& iothread)
const
{
//    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    const Parameters& parameters = field.parameters();

#if 0
    // w/o overlap //
    if ( parameters.coefCFROutput().is_fout_step(t) ) {
        const int  t_fout = parameters.coefCFROutput().t_fout(t);

        // copy_io //
        field.copy_io_field();
        field.write_field(t_fout);
    }
#else // w overlap //
    // thread : work //
    if ( parameters.coefCFROutput().is_fout_step(t) ) {
        const int  t_fout = parameters.coefCFROutput().t_fout(t);

        // copy_io //
        field.copy_io_field();


        // thread output //
        if (rank_ == 0) { std::cout << "io work : t fout = " << t_fout << std::endl; }
        iothread.work(true, [t_fout, &parameters, &field](){ field.write_field(t_fout); } );
    }

    // thread : join //
    if ( parameters.coefCFROutput().is_just_before_fout_step(t) || parameters.is_step_end(t) ) { // overlap //
        if (rank_ == 0) {
            std::cout << std::endl;
            std::cout << "io join : t = " << t << std::endl;
            }
        iothread.join(true);
    }
#endif
}


void OutputFunc::
output_monitor_data(int t, Field& field)
const
{
//    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    const Parameters& parameters = field.parameters();

    if ( parameters.coefCFROutput().is_fout_step(t) ) {
        const int  t_fout = parameters.coefCFROutput().t_fout(t);

        field.write_monitor(t_fout);
    }
}


void OutputFunc::
output_channel_flow_data(int t, Field& field)
const
{
    const Parameters& parameters = field.parameters();
    if ( parameters.coefCFROutput().is_cout_step(t) || parameters.is_step_end(t) ) {
//        field.copy_io_field();

        const int  t_cout = parameters.coefCFROutput().t_cout(t);

        PostprocessChannelFlow postprocess;
        postprocess.OutputStatistics(t_cout, field, field.meshValues0());
    }
}


void OutputFunc::
output_taylor_green_data(int t, Field& field)
const
{
    const Parameters& parameters = field.parameters();
    if ( parameters.coefCFROutput().is_cout_step(t) ) {
        const int  t_cout = parameters.coefCFROutput().t_cout(t);
        const real time   = parameters.time_now();

        PostprocessTaylorGreen postprocess;
        postprocess.OutputStatistics(t_cout, time, field, field.meshValues0());
    }
}


void OutputFunc::
output_natural_convection2d_data(int t, Field& field)
const
{
    const Parameters& parameters = field.parameters();
    if ( parameters.coefCFROutput().is_cout_step(t) ) {
        const int  t_cout = parameters.coefCFROutput().t_cout(t);
        const real time   = parameters.time_now();

        PostprocessNaturalConvection2d postprocess;
        postprocess.OutputStatistics(t_cout, field, field.meshValues0());
    }
}


void OutputFunc::
output_natural_convection3d_data(int t, Field& field)
const
{
    const Parameters& parameters = field.parameters();
    if ( parameters.coefCFROutput().is_cout_step(t) ) {
        const int  t_cout = parameters.coefCFROutput().t_cout(t);
        const real time   = parameters.time_now();

        PostprocessNaturalConvection3d postprocess;
        postprocess.OutputStatistics(t_cout, field, field.meshValues0());
    }
}


double OutputFunc::
cal_elapsed_time_msec (
    const struct timeval&   begin,
    const struct timeval&   end
    )
const
{
    // msec //
    return    (end.tv_sec  - begin.tv_sec) * 1000
            + (end.tv_usec - begin.tv_usec) / 1000.0;
}
