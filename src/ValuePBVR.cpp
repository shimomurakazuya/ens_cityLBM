#ifdef USE_PBVR
#include <cuda.h>
#include <chrono>
#include <sstream>

#include "ValuePBVR.h"
#include "invoker.hpp"
#include "range.hpp"
#include "Index.h"
#include "IndexLBM.h"
#include "Field.h"
#include "FuncAllocate.h"
#include "FuncAMRMesh.h"
#include "FuncObj.h"
#include "cuda_safe_call.hpp"
#include "runtime_error.hpp"
#include "mpi_wrapper.hpp"

#include "kvs_wrapper.h"

void ValuePBVR::init(const Field* p_field) {
    auto& field = *p_field;
    memType_ = MemType::Managed;

    for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
        auto num_tasks = field.taskID().num_task_array_lbm(lv); // task_lbm == 袖のぞくリーフ群
        n_leaves_levels_[lv] = num_tasks;
    }
    auto n_leaves_ = n_leaves();
    allocate(n_leaves_);

    init_stations(p_field);
}

void ValuePBVR::init_stations(const Field* p_field) {
    auto& field = *p_field;
    auto mpi = util::mpi(MPI_COMM_WORLD);

    if(mpi.rank() == 0) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }
    const auto n_leaf = field.tree().number_of_nodes();
    const auto nx = DefAMR::NX_LEAF;

    // set xyz //
    auto& st0 = stations_.at( 0); st0.id =  "1"; st0.x = -170.; st0.y =  340.; st0.z = 2.;
    auto& st1 = stations_.at( 1); st1.id =  "2"; st1.x =  -40.; st1.y =  330.; st1.z = 2.;
    auto& st2 = stations_.at( 2); st2.id =  "3"; st2.x =  140.; st2.y =  330.; st2.z = 2.;
    auto& st3 = stations_.at( 3); st3.id =  "4"; st3.x = -170.; st3.y =  220.; st3.z = 2.;
    auto& st4 = stations_.at( 4); st4.id =  "5"; st4.x =  -40.; st4.y =  220.; st4.z = 2.;
    auto& st5 = stations_.at( 5); st5.id =  "6"; st5.x =  100.; st5.y =  180.; st5.z = 2.;
    auto& st6 = stations_.at( 6); st6.id =  "7"; st6.x = -190.; st6.y =  110.; st6.z = 2.;
    auto& st7 = stations_.at( 7); st7.id =  "8"; st7.x =  -40.; st7.y =  110.; st7.z = 2.;
    auto& st8 = stations_.at( 8); st8.id =  "9"; st8.x =  120.; st8.y =  110.; st8.z = 2.;
    auto& st9 = stations_.at( 9); st9.id = "10"; st9.x =  -50.; st9.y =  -75.; st9.z = 2.;
    auto& sta = stations_.at(10); sta.id = "11"; sta.x =-1010.; sta.y =-1750.; sta.z = 2.;
    auto& stb = stations_.at(11); stb.id = "12"; stb.x = -196.; stb.y = -344.; stb.z = 2.;
    // xyz to rank_l_ijk_iijjkk
    for(auto&& st: stations_) {
        for(std::intptr_t l=0; l<n_leaf; l++) {
            if(st.rank >= 0) { continue; } // skip if station is already found //
            const auto& node = *field.tree().nodes(l);
            if(! node.nodeCalFlags().Cal()) { continue; } // skip halo //
            const auto level = node.level();
            const auto& meshValue = field.meshValue(level);
            const auto& coord = meshValue.coordinates();
            const auto offset = node.neighbor_mesh_offsets().offset0();
            const auto dx0 = coord.dx();
            const auto x0 = coord.x(offset);
            const auto y0 = coord.y(offset);
            const auto z0 = coord.z(offset);
            // skip nochance leaf //
            if(! (
                   x0 <= st.x && st.x <= x0 + dx0*nx
                && y0 <= st.y && st.y <= y0 + dx0*nx
                && z0 <= st.z && st.z <= z0 + dx0*nx
            )) { continue; }
            // seek this leaf //
            for(int k=0; k<nx; k++) for(int j=0; j<nx; j++) for(int i=0; i<nx; i++) {
                const auto x = x0 + dx0 * (0.5 + i); // cell center
                const auto y = y0 + dx0 * (0.5 + j); // cell center
                const auto z = z0 + dx0 * (0.5 + k); // cell center
                // if in the cell, decide xyz --> lijk //
                if( (
                      x-0.5*dx0 <= st.x && st.x <= x+0.5*dx0
                   && y-0.5*dx0 <= st.y && st.y <= y+0.5*dx0
                   && z-0.5*dx0 <= st.z && st.z <= z+0.5*dx0
                )) { // nearest neightbor torima
                    st.rank = mpi.rank();
                    st.l = l;
                    st.i = i;
                    st.j = j;
                    st.k = k;
                    st.ii = 0;
                    st.jj = 0;
                    st.kk = 0;
                    std::cout << " station " << st.id << " was found on "
                        << "rank " << st.rank
                        << ", leaf " << st.l
                        << ", ijk (" << st.i << ", " << st.j << ", " << st.k << std::endl;
                }
            }
        }
    }
    
    // allreduce
    for(auto&& st: stations_) {
        st.rank = mpi.reduce_max(st.rank);
        if(st.rank == mpi.rank()) {
            std::cout << " station " << st.id << " was found on "
                << "rank " << st.rank
                << ", leaf " << st.l
                << ", ijk (" << st.i << ", " << st.j << ", " << st.k << std::endl;
        }
        if(st.rank < 0 && mpi.rank() == 0) {
            std::cerr << " station " << st.id << " not found" << std::endl;
        }
    }
}

void ValuePBVR::allocate(const size_t n_leaves) {
    std::cout << __PRETTY_FUNCTION__ << ": n_leaves = " << n_leaves << std::endl;
    for ( auto&& v: variables_ ) {
        const auto nn_max = n_leaves*nx*ny*nz;
        FuncAllocate::allocate_value(&v, nn_max, memType_);
    }
    FuncAllocate::allocate_value(&cell_length_, n_leaves, memType_);
    FuncAllocate::allocate_value(&leaf_min_coord_, 3*n_leaves, memType_);

    //// cpu
    //for ( auto&& v: variables_cpu_ ) {
    //    const auto nn_max = n_leaves*nx*ny*nz;
    //    v = new float[nn_max];
    //}
    //cell_length_cpu_ = new float[n_leaves];
    //leaf_min_coord_cpu_ = new float[n_leaves];
}

void ValuePBVR::send_gpu_to_cpu() {
    throw STD_RUNTIME_ERROR("not-implemented function");
    //const auto nlb = n_leaves() * sizeof(float);
    //for(auto lv: util::irange(nvariables)) {
    //    CUDA_SAFE_CALL(cudaMemcpy(variables_cpu_[lv], variables_[lv], nlb*nx*ny*nz, cudaMemcpyDeviceToHost));
    //}
    //CUDA_SAFE_CALL(cudaMemcpy(cell_length_cpu_, cell_length_, nlb, cudaMemcpyDeviceToHost));
    //CUDA_SAFE_CALL(cudaMemcpy(leaf_min_coord_cpu_, leaf_min_coord_, 3*nlb, cudaMemcpyDeviceToHost));
}

void ValuePBVR::release() {
    for ( auto&& v: variables_ ) {
        FuncAllocate::release_value(v, memType_);
    }
    FuncAllocate::release_value(cell_length_, memType_);
    FuncAllocate::release_value(leaf_min_coord_, memType_);

    //// cpu
    //for ( auto&& v: variables_cpu_ ) {
    //    delete [] v;
    //}
    //delete [] cell_length_cpu_;
    //delete [] leaf_min_coord_cpu_;
}

void ValuePBVR::update(const Field* p_field) {
    auto& field = *p_field;

    size_t n_leaves_ofs = 0;
    for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
        const auto num_tasks = field.taskID().num_task_array_lbm(lv);
        const auto id_tasks = field.taskID().id_task_array_lbm(lv);
        const auto& meshValue = field.meshValue(lv);

        const auto* mesh_offsets3x3x3 = field.tree().mesh_offsets3x3x3();
        const auto* u_cell = meshValue.valueNS().u();
        const auto* v_cell = meshValue.valueNS().v();
        const auto* w_cell = meshValue.valueNS().w();
        const auto* T_cell = meshValue.valueNS().T();
        const auto* scalar_cell = meshValue.valueNS().scalar();
        const auto* x = meshValue.coordinates().x();
        const auto* y = meshValue.coordinates().y();
        const auto* z = meshValue.coordinates().z();
        const auto dx = meshValue.coordinates().dx();
        const auto c_ref = field.parameters().c_ref_lbm();

        auto* u_node      = variables_[0] + n_leaves_ofs * nx*ny*nz;
        auto* v_node      = variables_[1] + n_leaves_ofs * nx*ny*nz;
        auto* w_node      = variables_[2] + n_leaves_ofs * nx*ny*nz;
        auto* T_node      = variables_[3] + n_leaves_ofs * nx*ny*nz;
        auto* scalar_node = variables_[4] + n_leaves_ofs * nx*ny*nz;
        auto* cell_length = cell_length_;
        auto* leaf_min_coord = leaf_min_coord_;

        util::invoke_device<<<dim3(num_tasks), dim3(nx, ny, nz)>>>(
            [=]__device__() {
                auto i = threadIdx.x;
                auto j = threadIdx.y;
                auto k = threadIdx.z;
                auto l = blockIdx.x;
                auto idl = id_tasks[l];

                // get offset3d, ids
                int  offset3d[27];
                for (int idv=0; idv<27; idv++) {
                    offset3d    [idv] = mesh_offsets3x3x3[idl*27 + idv];
                }
                int ids[27];
                for (int kv=-1; kv<=1; kv++) { for (int jv=-1; jv<=1; jv++) { for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                    ids[idv_leaf] = Index::id(i+iv,j+jv,k+kv, offset3d);
                } } }

                // leaf info
                if(i == 0 && j == 0 && k == 0) {
                    const auto id0 = Index::id(ids, 0, 0, 0);
                    cell_length[l + n_leaves_ofs] = dx;
                    leaf_min_coord[0+  3*(l + n_leaves_ofs)] = x[id0];
                    leaf_min_coord[1+  3*(l + n_leaves_ofs)] = y[id0];
                    leaf_min_coord[2+  3*(l + n_leaves_ofs)] = z[id0];
                }

                // variables on 5x5x5 nodes
                double uu = 0, vv = 0, ww = 0, TT = 0, ss = 0;
                for(int kk=-1; kk<1; kk++) { for(int jj=-1; jj<1; jj++) { for(int ii=-1; ii<2; ii++) {
                    const auto id_cell = Index::id(ids, ii, jj, kk);
                    uu += u_cell[id_cell];
                    vv += v_cell[id_cell];
                    ww += w_cell[id_cell];
                    TT += T_cell[id_cell];
                    ss += scalar_cell[id_cell];
                } } }
                const auto id_node = i + nx*j + nx*ny*k + nx*ny*nz*l;
                u_node[id_node] = uu/8 * c_ref;
                v_node[id_node] = vv/8 * c_ref;
                w_node[id_node] = ww/8 * c_ref;
                T_node[id_node] = TT/8;
                scalar_node[id_node] = ss/8;
            }
        ); // invoke_device
        CUDA_SAFE_CALL(cudaDeviceSynchronize());

        // count concat
        n_leaves_ofs += num_tasks;
    } // for lv
    //send_gpu_to_cpu();
}

void ValuePBVR::output_pbvr_particles(const Field* p_field, int step) {
    int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
    auto& field = *p_field;
    auto& optionParser = field.optionParser();
    auto& parameters = field.parameters();

    //-------------------------//
    //---Start Visualization---//
    //-------------------------//

    //
    //全体の最大最小値を示すpfiファイルを生成
    //
    char filename[256];
    sprintf(filename,"./jupiter_particle_out/t_pfi_coords_minmax.txt");


    // 都市標高データのMin/Maxに合わせる
    const float x_global_domain_min = -598;
    const float y_global_domain_min = -748;
    const float z_global_domain_min =    0;

    const float x_global_domain_max =  598;
    const float y_global_domain_max =  748;
    const float z_global_domain_max =  150;


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
    auto tp_start = std::chrono::high_resolution_clock::now();
    update(&field);
    const int n_leaves = this->n_leaves();
    float* cell_length = this->cell_length();
    float* leaf_min_coord = this->leaf_min_coord();
    const float** variables = this->variables();

    const real c_ref = parameters.c_ref_lbm();

    auto tp_end = std::chrono::high_resolution_clock::now();

    if(mpi_rank == 0) {
        auto dur_nsec = std::chrono::duration_cast<std::chrono::nanoseconds >( tp_end - tp_start ).count();
        std::cout << "conv. CityLBM-->PBVR: duration = " << dur_nsec*1e-9 << " sec" << std::endl;
    }

    int resolution[4] = { nx, ny, nz, n_leaves };
    domain_parameters dom = {
        x_global_domain_min,
        y_global_domain_min,
        z_global_domain_min,
        x_global_domain_max,
        y_global_domain_max,
        z_global_domain_max,
        resolution
    };

    static int time_step = 0;
    static int count_step = 0;
    const  int iteration_number = 10;

    // check leaf_min_coord
    for( int l = 0; l < n_leaves; l+= 100 )
    {
        const float x = leaf_min_coord[ 3*l   ];
        const float y = leaf_min_coord[ 3*l+1 ];
        const float z = leaf_min_coord[ 3*l+2 ];

        if( -350 < x && x < 350 && 
            -400 < y && y < 700 &&
               0 < z && z < 300 )
        {
            std::cout << "(x,y,z)= ( "
                      << x << ", "
                      << y << ", "
                      << z << " )" <<std::endl;
        }
    }

    static MultiVariableScalars comb_mvscalars( nvariables );

    MultiVariableScalars mvscalars =
        generate_samples( count_step, dom,
                          cell_length, leaf_min_coord,
                          nvariables, variables );

    comb_mvscalars += mvscalars;

    count_step++;

    if( count_step == iteration_number )
    {
        output_samples( time_step, comb_mvscalars );

        generate_particles( time_step, dom,
                            cell_length, leaf_min_coord,
                            nvariables, variables );

        time_step++;
        count_step = 0;
        comb_mvscalars.reset();
    }

    //-------------------------//
    //----End Visualization----//
    //-------------------------//
}

void ValuePBVR::
output_scalar_timeseries(Field* p_field, int step) {
    auto& field = *p_field;
    auto& optionParser = field.optionParser();
    auto& parameters = field.parameters();

    auto mpi = util::mpi(MPI_COMM_WORLD);
    
    mpi.barrier();

    float Scalars[12];
    float Averages[12];
    // hasegawa 2020; include/OutpuFunc.hで設定した観測点の値を得るサンプルコード //
    for (int i=0; i<stations_.size(); i++) {
        const auto& st = stations_.at(i);
        // 瞬時値getter //
        real scalar = -99999; // mpi_max取るので、初期値はNANではなく十分小さい負の数
        if(st.rank == mpi.rank()) {
            const auto& node = *field.tree().nodes(st.l);
            const auto& level = node.level();
            const auto& offsets = node.neighbor_mesh_offsets();
            const auto* val_ptr = field.meshValue(level).valueNS().scalar();
            scalar = FuncAMRMesh::RawData(val_ptr, st.i, st.j, st.k, DefAMR::NX_LEAF, offsets);
        }
        scalar = mpi.reduce_max(scalar);
        Scalars[i] = scalar;
        // 時間平均値getter //
        real scalar_ave = -99999; // mpi_max取るので、初期値はNANではなく十分小さい負の数
        if(st.rank == mpi.rank()) {
            const auto& node = *field.tree().nodes(st.l);
            const auto& level = node.level();
            const auto& offsets = node.neighbor_mesh_offsets();
            const auto* val_ptr = field.valueStat(level).sc_sum();
            scalar_ave = FuncAMRMesh::RawData(val_ptr, st.i, st.j, st.k, DefAMR::NX_LEAF, offsets);
            scalar_ave /= field.valueStat(level).t_count();
        }
        scalar_ave = mpi.reduce_max(scalar_ave);
        Averages[i] = scalar_ave;
        // for check
        if(mpi.rank() == 0) {
            std::cout << " station " << st.id << ": scalar = " << scalar << ", scalar_ave = " << scalar_ave << std::endl;
        }
    }
    if(mpi.rank() == 0) {
        std::ofstream ofs( "station_values.txt" );
        for( int i=0; i<12; i++ ){
            ofs << Scalars[i] << "," << Averages[i] << std::endl;
        }
        ofs.close();
    }

    for(int level=0; level<DefAMR::LV_MAX; level++) {
        field.valueStat(level).zeroset(); // 時間平均値をリセットして取り直し始める
    }
}

void ValuePBVR::
steer_source_posi(Field* p_field, int step) {
    auto& field = *p_field;
    auto mpi = util::mpi(MPI_COMM_WORLD);

    //**********************************//
    // start Visualization and Steering //
    //**********************************//
           volatile real          source_posi[3] = { -99999, -99999, -99999 };
    static volatile real previous_source_posi[3] = { -99999, -99999, -99999 };

    std::string   filename = "source_posi.txt";
    std::ifstream fin( filename.c_str(), std::ios::in );

    if( !fin.is_open() )
    {
        std::cerr<<"Cannot open "<<filename<<std::endl;
        fin.close();
    }
    else
    {
        std::string line;
        int i = 0;
        while( std::getline( fin, line ) )
        {
            if(i==3) break;
            source_posi[i] = std::atof( line.c_str() );
            i++;
        }

        std::cout << "source_posi[0]=" << source_posi[0] << std::endl;
        std::cout << "source_posi[1]=" << source_posi[1] << std::endl;
        std::cout << "source_posi[2]=" << source_posi[2] << std::endl;

    }

    // hasegawa 2020; Steeringスキップ時、全MPIプロセスで更新判定をallreduceしてからreturn。
    // float の eq は怪しい？ 同じstd::atof()で呼んでいるから完全一致することを期待？
    int to_return_this_rank = 0;
    if( previous_source_posi[0] == source_posi[0] &&
        previous_source_posi[1] == source_posi[1] &&
        previous_source_posi[2] == source_posi[2]   )
    {
        to_return_this_rank = 1;
        // return;; // 運が悪いと一部のrankのみreturnしてしまい、デッドロックになる。のでコメントアウト。
    }
    int to_return_world = -1;
    MPI_Allreduce(&to_return_this_rank, &to_return_world, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    if(to_return_world > 0) 
    {
        std::stringstream ss;
        ss<<"XX @ rank " << mpi.rank() << " XX " <<std::endl;
        ss<<" Skip InitSource() "<<std::endl;
        ss<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<std::endl;
        std::cout << ss.str() << std::flush;
        return;
    }

    previous_source_posi[0] = source_posi[0];
    previous_source_posi[1] = source_posi[1];
    previous_source_posi[2] = source_posi[2];

    std::cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<std::endl;
    std::cout<<" CALLED InitSource FUNCTION"<<std::endl;
    std::cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<std::endl;


    //**********************************//
    //   end Visualization and Steering //
    //**********************************//

    const int  n_leaf  = field.tree().number_of_nodes();
    const int  nn_cell = DefAMR::NN_LEAF;

    int  num_source = 0;
    int  num_sources[4] = { 0, 0, 0, 0 };
    for (int i=0; i<n_leaf; i++) {
        const int  lv     = field.tree().nodes(i)->level();
        const int  offset = field.tree().nodes(i)->mesh_offset();

        // coordinates xyz //
        const real*  x  = &field.meshValue(lv).coordinates().x()[offset];
        const real*  y  = &field.meshValue(lv).coordinates().y()[offset];
        const real*  z  = &field.meshValue(lv).coordinates().z()[offset];
        const real   dx =  field.meshValue(lv).coordinates().dx();

        // initialize //
        const real*  lv_obj     = &field.meshValue(lv).valueObjLS().lv_obj()     [offset];
              real*  sc_scalar  = &field.meshValue(lv).valueObjLS().sc_scalar () [offset];
              real*     scalar  = &field.meshValue(lv).valueNS   ().   scalar () [offset];

              real*  sc_scalarn  = &field.meshValue_new(lv).valueObjLS().sc_scalar () [offset];
              real*     scalarn  = &field.meshValue_new(lv).valueNS   ().   scalar () [offset];

        const real  c_ref  = field.parameters().c_ref_lbm();
        for (int id=0; id<nn_cell; id++) {
            // kawamura 2020, hasegawa mod.
            // reset scalar //
            scalar [id] = 0;

            scalarn [id] = 0;

            // seek scalar source //
            const real x_tmp = x[id] + (real)0.5*dx;
            const real y_tmp = y[id] + (real)0.5*dx;
            const real z_tmp = z[id] + (real)0.5*dx;

//            const real box_min[] = { -(real)50.0-(real)0.5*dx, -(real)150.0-(real)0.5*dx, (real)0.0 };
//            const real box_max[] = { -(real)50.0+(real)0.5*dx, -(real)150.0+(real)0.5*dx,       dx  };

            const real box_min[] = { source_posi[0]-(real)0.5*dx, source_posi[1]-(real)0.5*dx, source_posi[2]-(real)0.5*dx };
            const real box_max[] = { source_posi[0]+(real)0.5*dx, source_posi[1]+(real)0.5*dx, source_posi[2]+(real)0.5*dx };
//            const real box_min[] = { source_posi[0]-(real)1.0*dx, source_posi[1]-(real)1.0*dx, source_posi[2]-(real)1.0*dx };
//            const real box_max[] = { source_posi[0]+(real)1.0*dx, source_posi[1]+(real)1.0*dx, source_posi[2]+(real)1.0*dx };

//            // tokyo //
//            const real box_min[] = { -(real)2.0*dx, -(real)2.0*dx, (real)10.0 + (real)8.0+(real)0.0*dx };
//            const real box_max[] = { +(real)2.0*dx, +(real)2.0*dx, (real)10.0 + (real)8.0+(real)4.0*dx  };

            const real xyz[] = { x_tmp, y_tmp, z_tmp };
            const real box_lv_obj = FuncObj::length_from_box(xyz, box_min, box_max);

            if ( FuncObj::is_fluid(lv_obj[id]) && FuncObj::is_obj(box_lv_obj) ) {
                sc_scalar[id] = (real)1.0/(dx*dx*dx);
                sc_scalarn[id] = (real)1.0/(dx*dx*dx);
                if ( field.tree().nodes(i)->nodeCalFlags().Cal() ) { num_source++; }
            }
            // kawamura 2020
            // added --else-- for reset pollution source.
            else
            {
                sc_scalar[id] = 0;
                sc_scalarn[id] = 0;
            }

        }

    }

    MPI_Request requ[1];
    int  num_source_g = 0;
    MPI_Iallreduce(&num_source,     &num_source_g,     1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, &requ[0]);
    MPI_Waitall(1, requ, MPI_STATUSES_IGNORE);

//    std::cout << rank_ << " : num_source, num_source_global = " << num_source << ", " << num_source_g << std::endl;
    std::cout << mpi.rank() << " : num_source, num_source_global = " << num_source << ", " << num_source_g << std::endl;

    mpi.barrier();
    if (mpi.rank() == 0) { std::cout << mpi.rank() << " : " << 0 << " :  num_source_global = " << num_source_g << std::endl; }
    if ( num_source_g == 0 ) { std::cout << "error : source not found " << __PRETTY_FUNCTION__ << std::endl;  }


    for (int l=0; l<n_leaf; l++) {
        const int  lv     = field.tree().nodes(l)->level();
        const int  offset = field.tree().nodes(l)->mesh_offset();

        real*  sc_scalar     = &field.meshValue    (lv).valueObjLS().sc_scalar ()  [offset];
        real*  sc_scalarn    = &field.meshValue_new(lv).valueObjLS().sc_scalar ()  [offset];

        for (int j=0; j<nn_cell; j++) {
            sc_scalar[j] = sc_scalar[j]/(float)num_source_g;
            sc_scalarn[j] = sc_scalar[j]/(float)num_source_g;
        }

    }
}

#endif // ifdef USE_PBVR
