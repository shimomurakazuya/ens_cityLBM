#include "OutputFunc.h"
#include <mpi.h>
#include <cmath>
#include <iomanip>
#include "defineFluidProperty.h"
#include "FuncMath.h"
#include "FuncObj.h"
#include "PostprocessChannelFlow.h"
#include "PostprocessTaylorGreen.h"
#include "PostprocessNaturalConvection2d.h"
#include "PostprocessNaturalConvection3d.h"
#include "MemoryUsage.hpp"

// include visualization
#include "kvs_wrapper.h" // kawamura
#include "FuncAMRMesh.h" // kawamura
#include <chrono> //kawamura
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
// end add include

#include <kvs/KVSMLObjectUnstructuredVolume>
#include <kvs/UnstructuredVolumeExporter>

void OutputFunc::
OutputFluidData(int t, Field& field, WorkerThread& iothread)
{
#if 1
    const Parameters& parameters = field.parameters();
    
    // output (cout_step) //
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


    // vtu, restart (fout_step, rout_step) //
    output_data(t, field, iothread);

    // postprocess monitor (mout_step) //
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
    if ( !comm_.is_rank0() ) { return; }

    std::cout << std::endl;
    std::cout << "t / step_end = " << t << " / " << parameters.step_end() << " ( time = " << parameters.time_now() << " / " << parameters.time_end() << " ) " << std::endl;
//    if ( t%20 == 0 ) { std::cout << std::flush; }
}


void OutputFunc::
cout_monitor(int t, const Field& field)
const
{
    const auto comm_target = comm_.col_vector().comm();
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


    const double dx_coarse = field.meshValue_new(0).coordinates().dx();
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

        const double dx     = field.meshValue_new(lv).coordinates().dx();
//        const double c_ref  = parameters.c_ref_lbm();

        const real*  rho = field.meshValue_new(lv).valueNS().rho();
        const real*  u   = field.meshValue_new(lv).valueNS().u();
        const real*  v   = field.meshValue_new(lv).valueNS().v();
        const real*  w   = field.meshValue_new(lv).valueNS().w();
        const real*  T   = field.meshValue_new(lv).valueNS().T();

        const real*  sgs_vis = field.meshValue_new(lv).valueNS().sgs_vis();

        const real*  scalar = field.meshValue_new(lv).valueNS().scalar();

        const real*  lv_obj = field.meshValue_new(lv).valueObjLS().lv_obj();
//        const real*  u_obj  = field.meshValue_new(lv).valueObjLS().u_obj();
//        const real*  v_obj  = field.meshValue_new(lv).valueObjLS().v_obj();
//        const real*  w_obj  = field.meshValue_new(lv).valueObjLS().w_obj();
        const Array3D<int>  offsets = tree.nodes(i)->neighbor_mesh_offsets();

        // monitor //
        const int  offset0 = offsets.offset0();

        for (int j=0; j<DefAMR::NN_LEAF; j++) {
            const int  id = j + offset0;

            if ( FuncObj::is_obj(lv_obj[id]) ) { continue; }

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

    MPI_Iallreduce(&num,          &num_g,          1, MPI_DOUBLE,  MPI_SUM, comm_target, &requ[0]);
    MPI_Iallreduce(&vel_max,      &vel_max_g,      1, MPI_DOUBLE,  MPI_MAX, comm_target, &requ[1]);
    MPI_Iallreduce(&rho_all,      &rho_all_g,      1, MPI_DOUBLE,  MPI_SUM, comm_target, &requ[2]);
    MPI_Iallreduce(&energy,       &energy_g,       1, MPI_DOUBLE,  MPI_SUM, comm_target, &requ[3]);
    MPI_Iallreduce(&scalar_all,   &scalar_all_g,   1, MPI_DOUBLE,  MPI_SUM, comm_target, &requ[4]);
    MPI_Iallreduce(&scalar_min,   &scalar_min_g,   1, MPI_DOUBLE,  MPI_MIN, comm_target, &requ[5]);
    MPI_Iallreduce(&scalar_max,   &scalar_max_g,   1, MPI_DOUBLE,  MPI_MAX, comm_target, &requ[6]);
    MPI_Iallreduce(&T_average,    &T_average_g,    1, MPI_DOUBLE,  MPI_SUM, comm_target, &requ[7]);
    MPI_Iallreduce(&T_min,        &T_min_g,        1, MPI_DOUBLE,  MPI_MIN, comm_target, &requ[8]);
    MPI_Iallreduce(&T_max,        &T_max_g,        1, MPI_DOUBLE,  MPI_MAX, comm_target, &requ[9]);
    MPI_Iallreduce(&sgs_vis_ave,  &sgs_vis_ave_g,  1, MPI_DOUBLE,  MPI_SUM, comm_target, &requ[10]);
    MPI_Iallreduce(&sgs_vis_max,  &sgs_vis_max_g,  1, MPI_DOUBLE,  MPI_MAX, comm_target, &requ[11]);

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
        MPI_Iallreduce(&n_upleaf[i],    &n_upleaf_g_sum[i],  1, MPI_INT,  MPI_SUM, comm_target, &requ_leaf[0 + 3*i]);
        MPI_Iallreduce(&n_upleaf[i],    &n_upleaf_g_max[i],  1, MPI_INT,  MPI_MAX, comm_target, &requ_leaf[1 + 3*i]);
        MPI_Iallreduce(&n_upleaf[i],    &n_upleaf_g_min[i],  1, MPI_INT,  MPI_MIN, comm_target, &requ_leaf[2 + 3*i]);
    }
    MPI_Waitall(3*DefAMR::LV_MAX,  requ_leaf,  MPI_STATUSES_IGNORE);


    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_.is_rank0()) {
        // memory //
        print_memory_info();
        
        // allocate //
        std::cout << "n_leaf(total)       = " << n_leaf_g << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_leaf_lv_g[i] << ", "; }
        std::cout << " ) ";
        std::cout << " : nn_global (10^6) = " << static_cast<float>(n_leaf_g) * DefAMR::NN_LEAF / 1.0e6 << std::endl;

        std::cout << "n_leaf(max@node)    = " << n_leaf_g_max << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_leaf_lv_g_max[i] << ", "; }
        std::cout << " ) ";
        std::cout << " : nn_max    (10^6) = " << static_cast<float>(n_leaf_g_max) * DefAMR::NN_LEAF / 1.0e6 << std::endl;

        std::cout << "n_leaf(min@node)    = " << n_leaf_g_min << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_leaf_lv_g_min[i] << ", "; }
        std::cout << " ) ";
        std::cout << " : nn_min    (10^6) = " << static_cast<float>(n_leaf_g_min) * DefAMR::NN_LEAF / 1.0e6 << std::endl;

        // updated //
        std::cout << "n_leaf(total, lv)    : updated           = "<< " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_sum[i] << ", "; }
        std::cout << " ) ";
        std::cout << "grid points(total, lv) 10^6 : updated    = "<< " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << static_cast<float>(n_upleaf_g_sum[i]) * DefAMR::NN_LEAF/1.0e6 << ", "; }
        std::cout << " ) " << std::endl;

        std::cout << "n_leaf(max@node, lv) : updated           = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_max[i] << ", "; }
        std::cout << " ) ";
        std::cout << "grid points(max@node, lv) 10^6 : updated = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << static_cast<float>(n_upleaf_g_max[i]) * DefAMR::NN_LEAF/1.0e6 << ", "; }
        std::cout << " ) " << std::endl;

        std::cout << "n_leaf(min@node, lv) : updated           = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << n_upleaf_g_min[i] << ", "; }
        std::cout << " ) ";
        std::cout << "grid points(min@node, lv) 10^6 : updated = " << " ( ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << static_cast<float>(n_upleaf_g_min[i]) * DefAMR::NN_LEAF/1.0e6 << ", "; }
        std::cout << " ) " << std::endl;


        auto default_precision = std::cout.precision();
        auto high_precision = 
            #ifdef COUT_MONITOR_PRECISION
            COUT_MONITOR_PRECISION;
            #else
            default_precision;
            #endif
        std::cout << "vel_max             = " << vel_max_g    << " :  ( " << vel_max_g*parameters.c_ref_lbm() << " m/s )"<< std::endl;
        std::cout << "rho_all             = " 
            << std::setprecision(high_precision) << rho_all_g << std::setprecision(default_precision)
            << std::endl;
        std::cout << "energy              = "
            << std::setprecision(high_precision) << energy_g << std::setprecision(default_precision)
            << std::endl;
        std::cout << "scalar(all,min,max) = " 
            << std::setprecision(high_precision) << scalar_all_g << std::setprecision(default_precision)
            << " ( " << scalar_min_g << ", " << scalar_max_g << " ) " << std::endl;
        std::cout << "T(ave,min,max)      = " 
            << std::setprecision(high_precision) << T_average_g << std::setprecision(default_precision)
            << " ( " << T_min_g << ", " << T_max_g << " ) " << std::endl;
        std::cout << "sgs_vis/vis         = " 
            << std::setprecision(high_precision) << sgs_vis_ave_g/kvis << std::setprecision(default_precision)
            << ", " << sgs_vis_max_g/kvis << " (average, max)"  << std::endl;

        std::cout << "dt(coarse)          = " << dt_coarse << " (sec.)" << std::endl;
        std::cout << "dx                  = ";
        for (int i=0; i<DefAMR::LV_MAX; i++) { std::cout << field.meshValue_new(i).coordinates().dx() << ", "; }
        std::cout << std::endl;
    }
}


void OutputFunc::
cout_timer(int t, Field& field)
{
    const auto comm_target = comm_.col_vector().comm();
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

        MPI_Iallreduce(&elapsed_time_cout, &elapsed_time_cout_g[0], 1, MPI_DOUBLE, MPI_SUM, comm_target, &requ[0]);
        MPI_Iallreduce(&elapsed_time_cout, &elapsed_time_cout_g[1], 1, MPI_DOUBLE, MPI_MIN, comm_target, &requ[1]);
        MPI_Iallreduce(&elapsed_time_cout, &elapsed_time_cout_g[2], 1, MPI_DOUBLE, MPI_MAX, comm_target, &requ[2]);

        MPI_Iallreduce(&mlups, &mlups_g[0], 1, MPI_DOUBLE, MPI_SUM, comm_target, &requ[3]);
        MPI_Iallreduce(&mlups, &mlups_g[1], 1, MPI_DOUBLE, MPI_MIN, comm_target, &requ[4]);
        MPI_Iallreduce(&mlups, &mlups_g[2], 1, MPI_DOUBLE, MPI_MAX, comm_target, &requ[5]);

        MPI_Iallreduce(&cal_time , &cal_time_g , 1, MPI_DOUBLE, MPI_SUM, comm_target, &requ[6]);
        MPI_Iallreduce(&comm_time, &comm_time_g, 1, MPI_DOUBLE, MPI_SUM, comm_target, &requ[7]);

        MPI_Waitall(8, requ, MPI_STATUSES_IGNORE);

        const auto ncpu = comm_.col_vector().size();
        elapsed_time_cout_g[0] /= ncpu;
        mlups_g[0]             /= ncpu;

        cal_time  /= ncpu;
        comm_time /= ncpu;


        if (comm_.is_rank0()) {
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
    const auto& coefCFROutput = parameters.coefCFROutput();

    const bool is_fout_step = coefCFROutput.is_fout_step(t);
    const bool is_rout_step = coefCFROutput.is_rout_step(t);

#if 0
    // w/o overlap //
    if( is_fout_step || is_rout_step ) {
        // copy_io //
        field.copy_io_field();

        const int  t_fout = coefCFROutput.t_fout(t);
        const int  t_rout = coefCFROutput.t_rout(t);
        if( is_fout_step ) { field.write_field_vtk(t_fout); }
        if( is_rout_step ) { field.write_field(t_rout); }
    }
#else // w overlap //
    if ( parameters.coefCFROutput().is_fout_step(t) ) 
    {
        const int  t_fout = parameters.coefCFROutput().t_fout(t);

        //-------------------------//
        //---Start Visualization---//
        //-------------------------//
        int mpi_rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

        //
        //全体の最大最小値を示すpfiファイルを生成
        //
        char filename[256];
        sprintf(filename,"./particle_out/t_pfi_coords_minmax.txt");

        const OptionParser& optionParser = field.optionParser();

////        //cityLBMのMin/Maxに合わせる
//        const float x_global_domain_min = -2048;
//        const float y_global_domain_min = -2048;
//        const float z_global_domain_min = -8;
//
//        const float x_global_domain_max =  2048;
//        const float y_global_domain_max =  2048;
//        const float z_global_domain_max =  2552;


//        // 建物周辺のMin/Maxに合わせる
        const float x_global_domain_min = -620;
        const float y_global_domain_min = -770;
        const float z_global_domain_min = -8;

        const float x_global_domain_max =  620;
        const float y_global_domain_max =  770;
        const float z_global_domain_max =  320;
//        const float z_global_domain_max =  150;


        // hasegawa 2020; mpirank==0のみfopen(). 全mpirankで書き込みのfopenは危険かも。
        if(mpi_rank == 0) {
            FILE* fp = fopen( filename, "w" );

            if( fp )
            {
                fprintf( fp, "%f %f %f %f %f %f\n",
                        x_global_domain_min,
                        y_global_domain_min,
                        z_global_domain_min,
                        x_global_domain_max,
                        y_global_domain_max,
                        z_global_domain_max );
                fclose( fp );
            }
        }

        //
        // Calculate Non sleeve leaves.
        //

        // Start measure time
        std::chrono::system_clock::time_point  start, end;
        if(mpi_rank == 0) {
            start = std::chrono::system_clock::now();
        }

        const Tree& tree = field.tree();
        const int   n_leaves = tree.number_of_nodes(); //袖領域を含めたリーフ数
              int   nl = 0;               //num. of leaves without sleeve.
              int*  non_sleeve_leaf_index = new int[ n_leaves ];
              int*  leaf_around_building  = new int[ n_leaves ];

        //袖領域を飛ばすために、重複の無いleafを数え上げる
        for (int l=0; l<n_leaves; l++)
        {
            non_sleeve_leaf_index[ l ] = -1;
            if( !tree.nodes(l)->nodeCalFlags().Cal() )
            {
                continue;
            }
            non_sleeve_leaf_index[ l ] = nl;
            nl++;
        }

        //
        // re-allocate from cell centered to vertex values
        //
        //float* cell_length = new float [ nl ];
        std::vector<float> cell_length( nl );
        //float* leaf_min_coord = new float [ 3*nl ];
        std::vector<float> leaf_min_coord( 3*nl );
        const int tmp_nl = nl; //num. of leaves without sleeve.を保存
        if(comm_.is_rank0())std::cout <<  mpi_rank << ", num. of leaves without sleeve = " << tmp_nl << std::endl;

        const int  nx = DefAMR::NX_LEAF;
        const int  ny = DefAMR::NX_LEAF;
        const int  nz = DefAMR::NX_LEAF;
        const int  nvariables = 5;

        //　建物周辺のリーフの数を数え上げる
        int n_leaves_around_building = 0;

        for (int l=0; l<n_leaves; l++)
        {
            leaf_around_building[ l ] = -1;
            if( non_sleeve_leaf_index[ l ] == -1 ){ continue; }

            const int index = non_sleeve_leaf_index[ l ];

            const int lv = tree.nodes(l)->level();

            const MeshValue& meshValue = field.meshValue( lv );

            const Array3D<int> offsets = tree.nodes( l )->neighbor_mesh_offsets();

            const int  offset     = offsets.offset0();
            const float leaf_min_coord_x = meshValue.coordinates().xnode(offset);
            const float leaf_min_coord_y = meshValue.coordinates().ynode(offset);
            const float leaf_min_coord_z = meshValue.coordinates().znode(offset);
            if(     -620 < leaf_min_coord_x  && leaf_min_coord_x < 620 
                &&  -770 < leaf_min_coord_y  && leaf_min_coord_y < 770
                &&  -8   < leaf_min_coord_z  && leaf_min_coord_z < 320 )
//                &&  -8   < leaf_min_coord_z  && leaf_min_coord_z < 150 )
            {
                leaf_around_building[ l ] = n_leaves_around_building;
                n_leaves_around_building++;
            }


        }//end of for l
        if(comm_.is_rank0())std::cout <<  mpi_rank << ", n_leaves_around_building = " << n_leaves_around_building << std::endl;
        

        // リーフ数を建物周辺のもので上書き
        nl = n_leaves_around_building;
        cell_length.resize(nl);
        leaf_min_coord.resize(3*nl);

        float** values = new float*[nvariables];
        for( int i=0; i<nvariables; i++ )
        {
            values[i] = new float[ (nx+1) * (ny+1) * (nz+1) * nl ];
        }

        const real c_ref = parameters.c_ref_lbm();

        // 4 unstruct kvs_wrapper
        int ncells = nl*nx*ny*nz;
        int nnodes = nl*(nx+1)*(ny+1)*(nz+1);
//        int nnodes = coords.size()/3;
        std::vector<float> coords;
        std::vector<unsigned int> connections;
        connections.resize(ncells*8);
        coords.resize(3*nl*(nx+1)*(ny+1)*(nz+1));

        for (int l=0; l<n_leaves; l++)
        {
//            std::cout << "leaf_around_building[ l ] = " << leaf_around_building[ l ] << std::endl;
//            if( non_sleeve_leaf_index[ l ] == -1 ){ continue; }
            if( leaf_around_building[ l ] == -1 ){ continue; }

//            const int index = non_sleeve_leaf_index[ l ];
            const int index = leaf_around_building[ l ];

            const int lv = tree.nodes(l)->level();

            const MeshValue& meshValue = field.meshValue( lv );

            const Array3D<int> offsets = tree.nodes( l )->neighbor_mesh_offsets();

            const int  offset     = offsets.offset0();

            cell_length   [   index   ] = meshValue.coordinates().dx();
            //leaf_min_coord[ 3*index   ] = meshValue.coordinates().x(offset);
            //leaf_min_coord[ 3*index+1 ] = meshValue.coordinates().y(offset);
            //leaf_min_coord[ 3*index+2 ] = meshValue.coordinates().z(offset);
            leaf_min_coord[ 3*index   ] = meshValue.coordinates().xnode(offset);
            leaf_min_coord[ 3*index+1 ] = meshValue.coordinates().ynode(offset);
            leaf_min_coord[ 3*index+2 ] = meshValue.coordinates().znode(offset);
            for (int k=0; k<nz+1; k++ )
            {
                for (int j=0; j<ny+1; j++ )
                {
                    for (int i=0; i<nx+1; i++)
                    {
                        const int cnt =
                            i + j*(nx+1) + k*(nx+1)*(ny+1) + index*(nx+1)*(ny+1)*(nz+1);
//                        values[0][cnt] =
                        values[1][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().u(),
                            i,j,k, nx, offsets ) * c_ref;
                        values[2][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().v(),
                            i,j,k, nx, offsets ) * c_ref;
                        values[3][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().w(),
                            i,j,k, nx, offsets ) * c_ref;
                        values[4][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().T(),
                            i,j,k, nx, offsets );
                        values[0][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().scalar(),
                            i,j,k, nx, offsets ) ;
//                            i,j,k, nx, offsets ) * 1000;
//                        速度の絶対値
//                        values[0][cnt] = sqrt(values[4][cnt]*values[4][cnt] + values[1][cnt]*values[1][cnt] + values[2][cnt]*values[2][cnt] ) ;
                            
                        coords[3*cnt+0]=leaf_min_coord[ 3*index   ] + i* meshValue.coordinates().dx() ;
                        coords[3*cnt+1]=leaf_min_coord[ 3*index+1 ] + j* meshValue.coordinates().dx() ;
                        coords[3*cnt+2]=leaf_min_coord[ 3*index+2 ] + k* meshValue.coordinates().dx() ;

                    }//end of for i
                }//end of for j
            }//end of for k
        }//end of for l

//        int nnodes = coords.size()/3;
//        std::cout << "coords.size = " << coords.size() << std::endl;
//        std::cout << __LINE__ <<std::endl;
        int line_size  = static_cast<int>( nx+1 );
        int slice_size = static_cast<int>( (nx+1) * (ny+1) );
        long vertex_index = 0;
        long connection_index = 0;
        for (int l=0; l<n_leaves; l++)
        {
//            if( non_sleeve_leaf_index[ l ] == -1 ){ continue; }
//            const int index = non_sleeve_leaf_index[ l ];
            if( leaf_around_building[ l ] == -1 ){ continue; }
            const int index = leaf_around_building[ l ];
            //long connection_index = index*(nx)*(ny)*(nz);
            //connection
            vertex_index= index * (nx+1) * (ny+1) * (nz+1);
            for (int k=0; k<nz; k++ )
            {
                for (int j=0; j<ny; j++ )
                {
                    for (int i=0; i<nx; i++)
                    {
                        const int local_vertex_index[8] =
                        {
                            vertex_index,
                            vertex_index + 1,
                            vertex_index + line_size,
                            vertex_index + line_size + 1,
                            vertex_index + slice_size,
                            vertex_index + slice_size + 1,
                            vertex_index + slice_size + line_size,
                            vertex_index + slice_size + line_size + 1
                        };
                            vertex_index++;

                            // hexahedra-1
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 0 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 1 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 3 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 2 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 4 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 5 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 7 ] );
                            connections[ connection_index++ ] = static_cast<unsigned int>( local_vertex_index[ 6 ] );

                    }//end of for i
                                    vertex_index++;
                }//end of for j
                                vertex_index+= line_size;
            }//end of for k
//                               vertex_index+= slice_size;
        }//end of for l

         // Start measure time
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@Reconstruction time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }

        int resolution[4] = { nx+1, ny+1, nz+1, nl };
        domain_parameters dom_unstruct = {
            x_global_domain_min,
            y_global_domain_min,
            z_global_domain_min,
            x_global_domain_max,
            y_global_domain_max,
            z_global_domain_max
        };

       static int time_step = 0;

//:      unstruct
//        ensemble_generate_particles( time_step, dom_unstruct,
//                            values, nvariables,
//                            coords.data(), nnodes,
//                            connections.data(), ncells, pbvr::VolumeObjectBase::CellType::Hexahedra );

//        // 各ノード(アンサンブルデータを出力)
        generate_particles( time_step, dom_unstruct,
                            values, nvariables,
                            coords.data(), nnodes,
                            connections.data(), ncells, pbvr::VolumeObjectBase::CellType::Hexahedra );


        time_step++;

//        delete [] non_sleeve_leaf_index;
//        delete [] leaf_around_building;
        //delete [] cell_length;
        //delete [] leaf_min_coord;

        // Start measure time
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@Generation time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }

# if 0
         float max; 
         float min; 
//
         max = -1000000;
         min =  1000000;
//
              const int N =(nx+1) * (ny+1) * (nz+1) * nl; 
              for (int i=0; i<N; i++)
              {
                 max =  std::fmax(values[0][i],max ); // vel_mag
                 min =  std::fmin(values[0][i],min ); // vel_mag

              }//end of for i
//
              std::cout << mpi_rank << ", max[0] = " << max << std::endl;
//              std::cout << "max[1] = " << max[1] << std::endl;
//              std::cout << "max[2] = " << max[2] << std::endl;
//              std::cout << "max[3] = " << max[3] << std::endl;
//              std::cout << "max[4] = " << max[4] << std::endl;
              std::cout << mpi_rank << ", min[0] = " << min << std::endl;
//              std::cout << "min[1] = " << min[1] << std::endl;
//              std::cout << "min[2] = " << min[2] << std::endl;
//              std::cout << "min[3] = " << min[3] << std::endl;
//              std::cout << "min[4] = " << min[4] << std::endl;


        kvs::ValueArray<float> tmp_coords(coords);
        kvs::ValueArray<unsigned int> tmp_connections(connections);
        std::vector<float> values_ar(nnodes);
        for (int i=0; i< nnodes; i++)
        {
            values_ar[i] = values[0][i];
        }
        kvs::ValueArray<float> tmp_values(values_ar);
		kvs::AnyValueArray array_values(tmp_values);
		//kvs::AnyValueArray array_values(values[0], nnodes);
        std::cout << "tmp_values.size() = " << tmp_values.size() << ", nnodes = " << nnodes << ", tmp_coords.size() = " << tmp_coords.size() << std::endl;
        for (int i= nnodes-10 ; i < nnodes ; i++)
        {
            std::cout << "array_values[i] = " << array_values.at<float>(i) << std::endl; 
        }

        kvs::UnstructuredVolumeObject* ave_volume = new kvs::UnstructuredVolumeObject(kvs::UnstructuredVolumeObject::CellType::Hexahedra
                ,nnodes,ncells,nvariables
                ,tmp_coords, tmp_connections, array_values);

        ave_volume -> updateMinMaxCoords(); 
        std::cout << *ave_volume << std::endl;
        if(comm_.ensemble_id() == 0)
        {
            std::cout << mpi_rank << ", min:" << ave_volume->minObjectCoord()   << ", max:" << ave_volume->maxObjectCoord() << std::endl;
            std::cout << mpi_rank << ", min:" << ave_volume->minExternalCoord() << ", max:" << ave_volume->maxExternalCoord() << std::endl;
        }

        //kvsml ファイル出力     
     std::stringstream ss;
     ss << "cityLBM_ens" << std::setfill('0') << std::setw(2) << comm_.ensemble_id(); 
     ss << "_";
     ss << std::setfill('0') << std::setw(5) << time_step;
     ss << "_";
     ss << std::setfill('0') << std::setw(7) << comm_.row_id();
     ss << "_";
     ss << std::setfill('0') << std::setw(7) << comm_.n_rows();
     ss << ".kvsml";

     if(comm_.is_rank0()) std::cout << "kvsml_file_name = " << ss.str() << std::endl;
//        const std::string& suffix = "cityLBM_" + std::to_string(i) + "_rank" + std::to_string(comm_.world().rank())  ".dat";
//   if(mpi_rank ==0)
//   {
       kvs::KVSMLObjectUnstructuredVolume* kvsml_object = new kvs::UnstructuredVolumeExporter<kvs::KVSMLObjectUnstructuredVolume>( ave_volume );
       kvsml_object->setWritingDataType( kvs::KVSMLObjectUnstructuredVolume::ExternalBinary );
       kvsml_object->write( ss.str() );
       delete kvsml_object;
//   }
        time_step++;
#endif
//         float max[5]; 
//         float min[5]; 
//
//            for (int i=0; i<5; i++)
//            {
//                max[i] = -1000000;
//                min[i] =  1000000;
//            }
//
//              const int N =(nx+1) * (ny+1) * (nz+1) * nl; 
//              for (int i=0; i<N; i++)
//              {
//                 max[0] =  std::fmax(values[0][i],max[0] ); // vel_mag
//                 max[1] =  std::fmax(values[1][i],max[1] ); // v
//                 max[2] =  std::fmax(values[2][i],max[2] ); // w 
//                 max[3] =  std::fmax(values[3][i],max[3] ); // T
//                 max[4] =  std::fmax(values[4][i],max[4] ); // u
//                 min[0] =  std::fmin(values[0][i],min[0] ); // vel_mag
//                 min[1] =  std::fmin(values[1][i],min[1] ); // v
//                 min[2] =  std::fmin(values[2][i],min[2] ); // w 
//                 min[3] =  std::fmin(values[3][i],min[3] ); // T
//                 min[4] =  std::fmin(values[4][i],min[4] ); // u
//
//              }//end of for i
//
//              std::cout << "max[0] = " << max[0] << std::endl;
//              std::cout << "max[1] = " << max[1] << std::endl;
//              std::cout << "max[2] = " << max[2] << std::endl;
//              std::cout << "max[3] = " << max[3] << std::endl;
//              std::cout << "max[4] = " << max[4] << std::endl;
//              std::cout << "min[0] = " << min[0] << std::endl;
//              std::cout << "min[1] = " << min[1] << std::endl;
//              std::cout << "min[2] = " << min[2] << std::endl;
//              std::cout << "min[3] = " << min[3] << std::endl;
//              std::cout << "min[4] = " << min[4] << std::endl;

        //-------------------------//
        //----End Visualization----//
        //-------------------------//
#if 0 
        // Start measure time
        if(mpi_rank == 0) {
            start = std::chrono::system_clock::now();
        }

        // 従来のデータ集約 
        // リーフ数を袖領域以外の値で上書き
        nl = tmp_nl;
        cell_length.resize(nl);
        leaf_min_coord.resize(3*nl);

        values = new float*[nvariables];
        for( int i=0; i<nvariables; i++ )
        {
            values[i] = new float[ (nx+1) * (ny+1) * (nz+1) * nl ];
        }

        ncells = nl*nx*ny*nz;
        nnodes = nl*(nx+1)*(ny+1)*(nz+1);
//        std::vector<float> coords;
//        std::vector<unsigned int> connections;
        connections.resize(ncells*8);
        coords.resize(3*nl*(nx+1)*(ny+1)*(nz+1));
        for (int l=0; l<n_leaves; l++)
        {
            if( non_sleeve_leaf_index[ l ] == -1 ){ continue; }

            const int index = non_sleeve_leaf_index[ l ];

            const int lv = tree.nodes(l)->level();

            const MeshValue& meshValue = field.meshValue( lv );

            const Array3D<int> offsets = tree.nodes( l )->neighbor_mesh_offsets();

            const int  offset     = offsets.offset0();

            cell_length   [   index   ] = meshValue.coordinates().dx();
            leaf_min_coord[ 3*index   ] = meshValue.coordinates().xnode(offset);
            leaf_min_coord[ 3*index+1 ] = meshValue.coordinates().ynode(offset);
            leaf_min_coord[ 3*index+2 ] = meshValue.coordinates().znode(offset);
            for (int k=0; k<nz+1; k++ )
            {
                for (int j=0; j<ny+1; j++ )
                {
                    for (int i=0; i<nx+1; i++)
                    {
                        const int cnt =
                            i + j*(nx+1) + k*(nx+1)*(ny+1) + index*(nx+1)*(ny+1)*(nz+1);
                        values[4][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().u(),
                            i,j,k, nx, offsets ) * c_ref;
                        values[1][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().v(),
                            i,j,k, nx, offsets ) * c_ref;
                        values[2][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().w(),
                            i,j,k, nx, offsets ) * c_ref;
                        values[3][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().T(),
                            i,j,k, nx, offsets );
                        values[0][cnt] =
                            FuncAMRMesh::CellToNode( meshValue.valueNS().scalar(),
                            i,j,k, nx, offsets );

                        coords[3*cnt+0]=leaf_min_coord[ 3*index   ] + i* meshValue.coordinates().dx() ;
                        coords[3*cnt+1]=leaf_min_coord[ 3*index+1 ] + j* meshValue.coordinates().dx() ;
                        coords[3*cnt+2]=leaf_min_coord[ 3*index+2 ] + k* meshValue.coordinates().dx() ;

                    }//end of for i
                }//end of for j
            }//end of for k
        }//end of for l
        // 集約処理
        const int N = (nx+1) * (ny+1) * (nz+1) * nl;
        std::vector<float> reduce_values_average(N);
        std::vector<float> reduce_values_average1(N);
        std::vector<float> reduce_values_average2(N);
        std::vector<float> reduce_values_average3(N);
        std::vector<float> reduce_values_average4(N);
        std::vector<float> reduce_values_test(N);
        std::vector<float> reduce_values_test1(N);
        std::vector<float> reduce_values_test2(N);
        std::vector<float> reduce_values_test3(N);
        std::vector<float> reduce_values_test4(N);
        std::vector<float> reduce_values_varience(N);
        std::vector<float> varience_values(N);
        std::vector<float> varience_values1(N);
        std::vector<float> varience_values2(N);
        std::vector<float> varience_values3(N);
        std::vector<float> varience_values4(N);
        int Nv = 5;
        int mpi_size = 36;
 
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@Reconstruction time 4 allreduce = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }
 
        comm_.world().reduce_sum_array(values[0], reduce_values_test, N );
        comm_.world().reduce_sum_array(values[1], reduce_values_test1, N );
        comm_.world().reduce_sum_array(values[2], reduce_values_test2, N );
        comm_.world().reduce_sum_array(values[3], reduce_values_test3, N );
        comm_.world().reduce_sum_array(values[4], reduce_values_test4, N );

        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@allnode_all_reduce = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }

        // アンサンブル方向の集約 
        comm_.row_vector().reduce_sum_array(values[0], reduce_values_test, N );
        comm_.row_vector().reduce_sum_array(values[1], reduce_values_test1, N );
        comm_.row_vector().reduce_sum_array(values[2], reduce_values_test2, N );
        comm_.row_vector().reduce_sum_array(values[3], reduce_values_test3, N );
        comm_.row_vector().reduce_sum_array(values[4], reduce_values_test4, N );
        
        // allgather
        std::vector<float> recv(N * comm_.n_rows());
        comm_.col_vector().reduce_gather_array(values[0], recv, N ); 
        comm_.col_vector().reduce_gather_array(values[1], recv, N ); 
        comm_.col_vector().reduce_gather_array(values[2], recv, N ); 
        comm_.col_vector().reduce_gather_array(values[3], recv, N ); 
        comm_.col_vector().reduce_gather_array(values[4], recv, N ); 
        

        

        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@allnode_reduce_gather = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }
 
#if 0
        // 集約用配列
     
         if(comm_.ensemble_id()==0 ) std::cout << "mpi_rank = " << mpi_rank  << ", N = " << N << std::endl;

        // 空間とアンサンブル方向について、集約のためのコミュニケータ分割が必要
        // アンサンブル数の確認
        const int n_col  = comm_.n_cols(); 
        if(comm_.is_rank0())std::cout << "n_cols() = " << comm_.n_cols() << std::endl; 
        comm_.row_vector().reduce_sum_array(values[0], reduce_values_average, N );
        comm_.row_vector().reduce_sum_array(values[1], reduce_values_average1, N );
        comm_.row_vector().reduce_sum_array(values[2], reduce_values_average2, N );
        comm_.row_vector().reduce_sum_array(values[3], reduce_values_average3, N );
        comm_.row_vector().reduce_sum_array(values[4], reduce_values_average4, N );


        // 集約処理
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@average_all_reduce time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }

//        for (int i=0; i< 5; i++)
//        { 
//            //for (auto& x : reduce_values_average) {
//            for (auto& x : recv) {
//                x /= n_col;
//            }
//        }

        for (int i=0; i< N; i++)
        { 
                reduce_values_average[i]  /= mpi_size;
                reduce_values_average1[i] /= mpi_size;
                reduce_values_average2[i] /= mpi_size;
                reduce_values_average3[i] /= mpi_size;
                reduce_values_average4[i] /= mpi_size;
        }


        // 集約処理
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@average_calc time = " << elapsed << " [msec]" << std::endl;
//            std::cout << "@@average_all_reduce time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }
 
#if 0
//        for (int i = 0;i < Nv; i++)
//        {
                for (int i=0; i< N; i++)
                {
                    varience_values[i]  = (reduce_values_average[i]  - values[0][i] )*(reduce_values_average[i]  - values[0][i] );
                    varience_values1[i] = (reduce_values_average1[i] - values[1][i] )*(reduce_values_average1[i] - values[1][i] );
                    varience_values2[i] = (reduce_values_average2[i] - values[2][i] )*(reduce_values_average2[i] - values[2][i] );
                    varience_values3[i] = (reduce_values_average3[i] - values[3][i] )*(reduce_values_average3[i] - values[3][i] );
                    varience_values4[i] = (reduce_values_average4[i] - values[4][i] )*(reduce_values_average4[i] - values[4][i] );
                } 
//        }

//        for (int i = 0;i < Nv; i++)
//        {
//            for (int i=0; i< recv.size(); i++)
//            {
//                varience_values[i] = (reduce_values_average[i] - values[0][i] )*(reduce_values_average[i] - values[0][i]);
//            } 
//        }


        for (int i = 0;i < Nv; i++)
        {
//            for (auto& x : varience_values) {
//                x /= n_col;
//            }
                for (int i=0; i< N; i++)
                {
                    varience_values[i]  /=mpi_size ;
                    varience_values1[i] /=mpi_size ;
                    varience_values2[i] /=mpi_size ;
                    varience_values3[i] /=mpi_size ;
                    varience_values4[i] /=mpi_size ;
                } 

        }

        // 集約処理
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@varience_calc time = " << elapsed << " [msec]" << std::endl;
//            std::cout << "@@average_all_reduce time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }
 
//        comm_.row_vector().reduce_sum_array(varience_values.data(), reduce_values_varience, N ); 
 
        comm_.world().reduce_sum_array(varience_values.data(), reduce_values_test, N );
        comm_.world().reduce_sum_array(varience_values.data(), reduce_values_test, N );
        comm_.world().reduce_sum_array(varience_values.data(), reduce_values_test, N );
        comm_.world().reduce_sum_array(varience_values.data(), reduce_values_test, N );
        comm_.world().reduce_sum_array(varience_values.data(), reduce_values_test, N );


        // 集約処理
        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@varience_allreduce time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }

#endif
//#pragma omp parallel
//        for (int i=0; i< N; i++)
//        {
//            values[0][i] = reduce_values_average[i] ; 
//        }
//    
//        generate_particles( time_step, dom_unstruct,
//                            values, nvariables,
//                            coords.data(), nnodes,
//                            connections.data(), ncells, pbvr::VolumeObjectBase::CellType::Hexahedra );

        if(mpi_rank == 0) {
            end = std::chrono::system_clock::now();
            double elapsed =
                std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
            std::cout << "@@Generation_allreduce time = " << elapsed << " [msec]" << std::endl;
            start = std::chrono::system_clock::now();
        }

#endif

        time_step++;

//        if(comm_.row_id() == 0 )
//        {
//            for (int i =N-10; i< N; i++ )
//            {
//                std::cout  << "values["<< i << "]  =" << values[0][i]  << std::endl;
//            }
//        }
//        if(comm_.is_rank0() ) 
//        {
//             for (int i =N-10; i< N; i++ )
//            {
//                std::cout  << "reduce_values_average["<< i << "]  =" <<  reduce_values_average[i]  << std::endl;
//            }
//        }
        for( int i=0; i<nvariables; i++ )
        {
            delete [] values[i];
        }
        delete [] values;
#else
        for( int i=0; i<nvariables; i++ )
        {
            delete [] values[i];
        }
        delete [] values;

#endif
        delete [] non_sleeve_leaf_index;
        delete [] leaf_around_building;
    
    }


    // thread : work //
    if( is_fout_step || is_rout_step ) {
        // copy_io //
        field.copy_io_field();

        const int  t_fout = coefCFROutput.t_fout(t);
        const int  t_rout = coefCFROutput.t_rout(t);

        // thread output //
        if(comm_.is_rank0()) { std::cout << "io work : t fout, rout = " << t_fout << ", " << t_rout << std::endl; }

#ifndef NO_IOTHREAD
        iothread.work(true, [is_fout_step, is_rout_step, t_fout, t_rout, &field](){
            if( is_fout_step ) { field.write_field_vtk(t_fout); }
            if( is_rout_step ) { field.write_field(t_rout); }
            } );
#else
        if( is_fout_step ) { field.write_field_vtk(t_fout); }
        if( is_rout_step ) { field.write_field(t_rout); }
#endif
    }

    // thread : join //
    if( coefCFROutput.is_just_before_fout_step(t) || coefCFROutput.is_just_before_rout_step(t) || parameters.is_step_end(t) ) { // overlap //
        if(comm_.is_rank0()) {
            std::cout << std::endl;
            std::cout << "io join : t = " << t << std::endl;
            }
        iothread.join(true);
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}


void OutputFunc::
output_monitor_data(int t, Field& field)
const
{
//    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    const Parameters& parameters = field.parameters();

    if ( parameters.coefCFROutput().is_mout_step(t) ) {
        const int  t_mout = parameters.coefCFROutput().t_mout(t);

        field.write_monitor(t_mout);
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
