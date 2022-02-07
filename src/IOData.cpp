
#include <string>
#include <zlib.h>

#include "IOData.h"
#include "defineFilenames.h"
#include "defineAMR.h"
#include "FuncLoop.h"
#include "FuncObj.h"
#include "FuncMath.h"
#include "FuncAMRMesh.h"
//#include "WriteHDF5.h"
//#include "ReadHDF5.h"
#include "mpi_wrapper.hpp"

// public //
//void  IOData::
//writeHDF5Values (
//    const Grid*         grids,
//    const Tree&         tree,
//    const MeshValue*    meshValues
//    )
//const
//{
//    const int  n_leaf = tree.number_of_nodes();
//    const int  nx_cell[3] = { DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF };
//
//
//    std::vector<stHDF5Values> HDF5values;
//    std::vector<int>          lv_leaf;
//    std::vector<bool>         flag_visualization;
//
//
//    for (int l=0; l<n_leaf; l++) {
//        const int  lv = tree.nodes(l)->level();
//        const int  nx = DefAMR::NX_LEAF;
//        const int  ny = DefAMR::NX_LEAF;
//        const int  nz = DefAMR::NX_LEAF;
//
//        const Array3D<int>  offsets     = tree.nodes(l)->neighbor_mesh_offsets();
//        const Array3D<int>  offsets_lbm = tree.nodes(l)->neighbor_lbm_mesh_offsets();
//
//        // HDF5 //
//        HDF5values.emplace_back( makeHDF5Values(nx, ny, nz, offsets, offsets_lbm, tree.nodes(l), meshValues[lv]) );
//
//
//        // flag //
//        lv_leaf.push_back( tree.nodes(l)->level() );
//        flag_visualization.push_back( tree.nodes(l)->nodeCalFlags().Cal() );
//    }
//
//
//    const std::string  folder        = Foldernames::io_folder;
//    const std::string  filename      = Filenames::hdf5_value_name0          + "-rank" + std::to_string(rank_) + "-step" + std::to_string(step_);
//    const std::string  filename_lbm  = Filenames::hdf5_value_name0 + "_lbm" + "-rank" + std::to_string(rank_) + "-step" + std::to_string(step_);
//
//    WriteHDF5::write_hdf_datasets (
//        folder,
//        filename,
//        filename_lbm,
//        hdf5_grid_name_,
//        HDF5values,
//        lv_leaf,
//        flag_visualization,
//        nx_cell
//        );
//}
//
//
//void  IOData::
//readHDF5Values (
//    const Grid*         grids,
//    const Tree&         tree,
//          MeshValue*    meshValues
//    )
//const
//{
//    const std::string  folder        = Foldernames::io_folder;
//    const std::string  filename      = Filenames::hdf5_value_name0          + "-rank" + std::to_string(rank_) + "-step" + std::to_string(step_);
//    const std::string  filename_lbm  = Filenames::hdf5_value_name0 + "_lbm" + "-rank" + std::to_string(rank_) + "-step" + std::to_string(step_);
//
//    const int  n_leaf = tree.number_of_nodes();
//    const int  nx_cell[3] = { DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF };
//
//    // read //
//    hid_t   file_id     = H5Fopen( (folder + "/" + filename     + ".h5").c_str(),  H5F_ACC_RDONLY, H5P_DEFAULT);
//    hid_t   file_id_lbm = H5Fopen( (folder + "/" + filename_lbm + ".h5").c_str(),  H5F_ACC_RDONLY, H5P_DEFAULT);
//
//    for (int l=0; l<n_leaf; l++) {
//        stHDF5Values  HDF5value = ReadHDF5::read_hdf_datasets(file_id, file_id_lbm, l);
//
//        // update //
//        const int  lv     = tree.nodes(l)->level();
//        const int  offset = tree.nodes(l)->mesh_offset();
//
//        // coordinates xyz //
//        real*  x = &meshValues[lv].coordinates().x()[offset];
//        real*  y = &meshValues[lv].coordinates().y()[offset];
//        real*  z = &meshValues[lv].coordinates().z()[offset];
//
//        const int  nnp = (nx_cell[0]+1)*(nx_cell[1]+1)*(nx_cell[2]+1);  // Node //
//        real*  xp_tmp = new real[nnp];
//        real*  yp_tmp = new real[nnp];
//        real*  zp_tmp = new real[nnp];
//
//        FuncMath::copy_i(xp_tmp, HDF5value.xp.get(), nnp);
//        FuncMath::copy_i(yp_tmp, HDF5value.yp.get(), nnp);
//        FuncMath::copy_i(zp_tmp, HDF5value.zp.get(), nnp);
//
//        meshValues[lv].coordinates().set_x_from_xp(x, xp_tmp, nx_cell[0], nx_cell[1], nx_cell[2]);
//        meshValues[lv].coordinates().set_x_from_xp(y, yp_tmp, nx_cell[0], nx_cell[1], nx_cell[2]);
//        meshValues[lv].coordinates().set_x_from_xp(z, zp_tmp, nx_cell[0], nx_cell[1], nx_cell[2]);
//
////        meshValues[lv].coordinates().set_x_from_xp(x, HDF5value.xp.get(), nx_cell[0], nx_cell[1], nx_cell[2]);
////        meshValues[lv].coordinates().set_x_from_xp(y, HDF5value.yp.get(), nx_cell[0], nx_cell[1], nx_cell[2]);
////        meshValues[lv].coordinates().set_x_from_xp(z, HDF5value.zp.get(), nx_cell[0], nx_cell[1], nx_cell[2]);
//
//        delete [] xp_tmp;
//        delete [] yp_tmp;
//        delete [] zp_tmp;
//
//
//        const int  strides[] = { 1, DefAMR::NX_LEAF, DefAMR::NX_LEAF*DefAMR::NX_LEAF };
//        meshValues[lv].coordinates().reset_dx(strides);
//
//
//        // values //
//        const int  nx = nx_cell[0];
//        const int  ny = nx_cell[1];
//        const int  nz = nx_cell[2];
//        const int  nn = nx*ny*nz;
//
//        FuncLoop::Loop3d  loop3d(0,nx, 0,ny, 0,nz);
//        loop3d.for_each( [&](const int i, const int j, const int k)
//            {
//                const int  id  = FuncLoop::id(i,j,k, nx  ,ny  ,nz  );
//                const int  idp = FuncLoop::id(i,j,k, nx+1,ny+1,nz+1);
//
//                meshValues[lv].valueObjLS().lv_obj() [offset + id] = HDF5value.lv_obj. get()[idp];
//                meshValues[lv].valueObjLS().rho_obj()[offset + id] = HDF5value.rho_obj.get()[idp];
//                meshValues[lv].valueObjLS().u_obj()  [offset + id] = HDF5value.u_obj.  get()[idp];
//                meshValues[lv].valueObjLS().v_obj()  [offset + id] = HDF5value.v_obj.  get()[idp];
//                meshValues[lv].valueObjLS().w_obj()  [offset + id] = HDF5value.w_obj.  get()[idp];
//
//                meshValues[lv].valueNS().u()  [offset + id] = HDF5value.u.  get()[idp];
//                meshValues[lv].valueNS().v()  [offset + id] = HDF5value.v.  get()[idp];
//                meshValues[lv].valueNS().w()  [offset + id] = HDF5value.w.  get()[idp];
//                meshValues[lv].valueNS().rho()[offset + id] = HDF5value.rho.get()[idp];
//
//                meshValues[lv].valueNS().scalar()[offset + id] = HDF5value.scalar.get()[idp];
//            } );
//
//        // lbm values //
//        const int  nQ         = LBM_velocity_model::nQ;
//        const int  offset_lbm = tree.nodes(l)->neighbor_lbm_mesh_offsets().offset0();
//        FuncMath::copy_i( &meshValues[lv].valueLBM().f_lbm()[offset_lbm], HDF5value.f_lbm.get(), nn*nQ );
//
//        // grid attribute values //
//    }
//
//    H5Fclose(file_id);
//    H5Fclose(file_id_lbm);
//}


void IOData::
writeBinaries(
    const Grid*          grids,
    const Tree&          tree,
    const MeshValue*     meshValues0,
    const MeshValue*     meshValues1,
    const ValueTimeAverage* valueTimeAverage
) const
{
    #ifdef NO_IODATA_RESTART
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skiped @ NO_IODATA_RESTART" << std::endl; }
    #else
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": write restarter t=" << step_ << std::endl; }
    const std::string& prefix = Foldernames::io_folder + "/restart/" + std::to_string(step_);

    // timestep
    if(comm_.is_rank0()) {
        auto&& file = std::ofstream(prefix + "/timestep.txt");
        if(!file) { throw STD_RUNTIME_ERROR("failed to open restart/timestep.txt" ); }
        file << step_;
    }

    // meshValues
    const MeshValue* meshValuess[2] = {meshValues0, meshValues1};
    for(int i=0; i<2; i++) {
        if( meshValuess[i] == nullptr ) { continue; }
        const MeshValue* meshValues = meshValuess[i];
        for(int lv=0; lv<DefAMR::LV_MAX; lv++) {
            const std::string& suffix = "_ptr" + std::to_string(i) + "_rank" + std::to_string(comm_.world().rank()) + "_lv" + std::to_string(lv) + ".dat";
            const std::intptr_t& nn_max = tree.number_of_nodes_lv(lv) * DefAMR::NN_LEAF;
            const MeshValue& meshValue = meshValues[lv];
            const int n_scalars = meshValue.n_scalars();
            auto&& write = [&](const std::string& name, const char* p, const size_t sizeof_type, const int n=1) {
                #ifndef NO_IODATA_GZBIN
                auto&& filename = prefix + "/" + name + suffix + ".gz";
                gzFile gz = gzopen(filename.c_str(), "wb1"); // 0: nocompression, 1: fast, 9: max compression
                try {
                    runtime_assert(gz != NULL, "FileIOError");
                    runtime_assert(gzwrite(gz, p, nn_max * n * sizeof_type)  == nn_max * n * sizeof_type, "FileIOError");
                } catch (std::runtime_error e) {
                    std::cerr << e.what() 
                        << " at: " << filename 
                        << " (" << nn_max << " " << sizeof_type << " " << nn_max*sizeof_type << ")"
                        << std::endl;
                }
                gzclose(gz);
                #else
                auto&& filename = prefix + "/" + name + suffix;
                FILE* fp = fopen(filename.c_str(), "wb");
                try {
                    runtime_assert(fp != NULL, "FileIOError");
                    runtime_assert(fwrite(p, sizeof_type, nn_max * n, fp) == nn_max * n, "FileIOError");
                } catch(std::runtime_error e) {
                    std::cerr << e.what()
                        << " at: " << filename 
                        << " (" << nn_max*n << " " << sizeof_type << " " << nn_max*n*sizeof_type << ")"
                        << std::endl;
                }
                fclose(fp);
                #endif
            };
            #define WRITE_BINARY_0()
            #define WRITE_BINARY_1(SET)
            #define WRITE_BINARY_2(SET, VAL) write( #VAL , \
                    reinterpret_cast<const char*>(meshValue. SET (). VAL ()), \
                    sizeof(decltype(*(meshValue. SET (). VAL ()))) \
            )
            #define WRITE_BINARY_3(SET, VAL, N) write( #VAL , \
                    reinterpret_cast<const char*>(meshValue. SET (). VAL ()), \
                    sizeof(decltype(*(meshValue. SET (). VAL ()))), \
                    N \
            )
            #define WRITE_BINARY_(arg1, arg2, arg3, arg4, arg5, ...) arg5
            #define WRITE_BINARY(...) WRITE_BINARY_(,##__VA_ARGS__,\
                                                    WRITE_BINARY_3(__VA_ARGS__),\
                                                    WRITE_BINARY_2(__VA_ARGS__),\
                                                    WRITE_BINARY_1(__VA_ARGS__),\
                                                    WRITE_BINARY_0(__VA_ARGS__)\
                                                   )
            // coordinate
            WRITE_BINARY(coordinates, i);
            WRITE_BINARY(coordinates, j);
            WRITE_BINARY(coordinates, k);
            {
                using real_type = Coordinates::real_type;
                const auto& coordinates = meshValue.coordinates();
                const real_type x0 = coordinates.x0();
                const real_type y0 = coordinates.y0();
                const real_type z0 = coordinates.z0();
                const real_type dx = coordinates.dx();
                const real_type dy = coordinates.dy();
                const real_type dz = coordinates.dz();
                auto&& filename = prefix + "/xyz0_dxdydz" + suffix;
                FILE* fp = fopen(filename.c_str(), "wb");
                fwrite(&x0, sizeof(real_type), 1, fp);
                fwrite(&y0, sizeof(real_type), 1, fp);
                fwrite(&z0, sizeof(real_type), 1, fp);
                fwrite(&dx, sizeof(real_type), 1, fp);
                fwrite(&dy, sizeof(real_type), 1, fp);
                fwrite(&dz, sizeof(real_type), 1, fp);
                fclose(fp);
            }
            // valueNS
            WRITE_BINARY(valueNS, u);
            WRITE_BINARY(valueNS, v);
            WRITE_BINARY(valueNS, w);
            WRITE_BINARY(valueNS, rho);
            WRITE_BINARY(valueNS, scalar, n_scalars);
            WRITE_BINARY(valueNS, sgs_vis);
            WRITE_BINARY(valueNS, T);
            // valueObjLS
            WRITE_BINARY(valueObjLS, lv_obj);
            WRITE_BINARY(valueObjLS, pad_obj);
            WRITE_BINARY(valueObjLS, rho_obj);
            WRITE_BINARY(valueObjLS, scalar_obj, n_scalars);
            WRITE_BINARY(valueObjLS, u_obj);
            WRITE_BINARY(valueObjLS, v_obj);
            WRITE_BINARY(valueObjLS, w_obj);
            WRITE_BINARY(valueObjLS, T_obj);
            WRITE_BINARY(valueObjLS, scalarn_obj, n_scalars);
            WRITE_BINARY(valueObjLS, rhon_obj);
            WRITE_BINARY(valueObjLS, un_obj);
            WRITE_BINARY(valueObjLS, vn_obj);
            WRITE_BINARY(valueObjLS, wn_obj);
            WRITE_BINARY(valueObjLS, Tn_obj);
            WRITE_BINARY(valueObjLS, hflux_z_obj);
            WRITE_BINARY(valueObjLS, hfluxn_z_obj);
            WRITE_BINARY(valueObjLS, bcTypes_f);
            WRITE_BINARY(valueObjLS, dirichlet_weight);
            WRITE_BINARY(valueObjLS, sc_scalar, n_scalars);
            #undef WRITE_BINARY
            // valueLBM
            for(int q=0; q<LBM_velocity_model::nQ; q++) {
                write( "lbm_" + std::to_string(q), 
                        reinterpret_cast<const char*>(meshValue.valueLBM().f_lbm() + q*nn_max),
                        sizeof(real)
                );
            }
            // valueTimeAverage
            write("u_mean", reinterpret_cast<const char*>(valueTimeAverage[lv].u_mean()), sizeof(real));
            write("v_mean", reinterpret_cast<const char*>(valueTimeAverage[lv].v_mean()), sizeof(real));
            write("w_mean", reinterpret_cast<const char*>(valueTimeAverage[lv].w_mean()), sizeof(real));
            write("T_mean", reinterpret_cast<const char*>(valueTimeAverage[lv].T_mean()), sizeof(real));
            write("vel2_fluc", reinterpret_cast<const char*>(valueTimeAverage[lv].vel2_fluc()), sizeof(real));
        }
    }
    #endif // NO_IODATA_RESTART
}


void IOData::
readBinaries(
    const Grid*          grids,
    const Tree&          tree,
          MeshValue*     meshValues0,
          MeshValue*     meshValues1,
          ValueTimeAverage* valueTimeAverage
) const
{
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": read restarter t=" << step_ << std::endl; }
    const std::string& prefix = Foldernames::io_folder + "/restart/" + std::to_string(step_);

    // timestep check
    {
        auto&& file = std::ifstream(prefix + "/timestep.txt");
        if(!file) { throw STD_RUNTIME_ERROR("failed to open restart/timestep.txt" ); }
        int t;
        file >> t;
        if(step_ != t) { throw STD_RUNTIME_ERROR("invalid timestep of restart"); }
    }

    // meshValues
    MeshValue* meshValuess[2] = {meshValues0, meshValues1};
    for(int i=0; i<2; i++) {
        if( meshValuess[i] == nullptr ) { continue; }
        MeshValue* meshValues = meshValuess[i];
        for(int lv=0; lv<DefAMR::LV_MAX; lv++) {
            const std::string& suffix = "_ptr" + std::to_string(i) + "_rank" + std::to_string(comm_.world().rank()) + "_lv" + std::to_string(lv) + ".dat";
            const std::intptr_t& nn_max = tree.number_of_nodes_lv(lv) * DefAMR::NN_LEAF;
            MeshValue& meshValue = meshValues[lv];
            const int n_scalars = meshValue.n_scalars();
            auto&& read = [&](const std::string& name, char* p, const size_t sizeof_type, const int n=1) {
                #ifndef NO_IODATA_GZBIN
                auto&& filename = prefix + "/" + name + suffix + ".gz";
                try {
                    gzFile gz = gzopen(filename.c_str(), "rb");
                    runtime_assert(gz != NULL, "FileIOError");
                    runtime_assert(gzread(gz, p, nn_max * n * sizeof_type)  == nn_max * n * sizeof_type, "FileIOError");
                    gzclose(gz);
                } catch (std::runtime_error e) {
                    std::cerr << "restart read: failed to open " << filename << std::endl;
                    throw;
                }
                #else
                auto&& filename = prefix + "/" + name + suffix;
                try {
                    FILE* fp = fopen(filename.c_str(), "rb");
                    runtime_assert(fp != NULL, "FileIOError");
                    runtime_assert(fread(p, sizeof_type, nn_max * n, fp) == nn_max * n, "FileIOError");
                } catch (std::runtime_error e) {
                    std::cerr << "restart read: failed to open " << filename << std::endl;
                    throw;
                }
                #endif
            };
            #define READ_BINARY_0()
            #define READ_BINARY_1(SET)
            #define READ_BINARY_2(SET, VAL) read( #VAL , \
                    reinterpret_cast<char*>(meshValue. SET (). VAL ()), \
                    sizeof(decltype(*(meshValue. SET (). VAL ()))) \
            )
            #define READ_BINARY_3(SET, VAL, N) read( #VAL , \
                    reinterpret_cast<char*>(meshValue. SET (). VAL ()), \
                    sizeof(decltype(*(meshValue. SET (). VAL ()))), \
                    N \
            )
            #define READ_BINARY_(arg1, arg2, arg3, arg4, arg5, ...) arg5
            #define READ_BINARY(...) READ_BINARY_(,##__VA_ARGS__,\
                                                  READ_BINARY_3(__VA_ARGS__),\
                                                  READ_BINARY_2(__VA_ARGS__),\
                                                  READ_BINARY_1(__VA_ARGS__),\
                                                  READ_BINARY_0(__VA_ARGS__)\
            )
            // coordinates
            READ_BINARY(coordinates, i);
            READ_BINARY(coordinates, j);
            READ_BINARY(coordinates, k);
            {
                auto&& filename = prefix + "/xyz0_dxdydz" + suffix;
                FILE* fp = fopen(filename.c_str(), "rb");
                using real_type = Coordinates::real_type;
                real_type x0, y0, z0, dx, dy, dz;
                fread(&x0, sizeof(real_type), 1, fp);
                fread(&y0, sizeof(real_type), 1, fp);
                fread(&z0, sizeof(real_type), 1, fp);
                fread(&dx, sizeof(real_type), 1, fp);
                fread(&dy, sizeof(real_type), 1, fp);
                fread(&dz, sizeof(real_type), 1, fp);
                meshValue.coordinates().set_x0(x0, y0, z0);
                meshValue.coordinates().set_dx(dx, dy, dz);
                fclose(fp);
            }
            // valueNS
            READ_BINARY(valueNS, u);
            READ_BINARY(valueNS, v);
            READ_BINARY(valueNS, w);
            READ_BINARY(valueNS, rho);
            READ_BINARY(valueNS, scalar, n_scalars);
            READ_BINARY(valueNS, sgs_vis);
            READ_BINARY(valueNS, T);
            // valueObjLS
            READ_BINARY(valueObjLS, lv_obj);
            READ_BINARY(valueObjLS, pad_obj);
            READ_BINARY(valueObjLS, rho_obj);
            READ_BINARY(valueObjLS, scalar_obj, n_scalars);
            READ_BINARY(valueObjLS, u_obj);
            READ_BINARY(valueObjLS, v_obj);
            READ_BINARY(valueObjLS, w_obj);
            READ_BINARY(valueObjLS, T_obj);
            READ_BINARY(valueObjLS, scalarn_obj, n_scalars);
            READ_BINARY(valueObjLS, rhon_obj);
            READ_BINARY(valueObjLS, un_obj);
            READ_BINARY(valueObjLS, vn_obj);
            READ_BINARY(valueObjLS, wn_obj);
            READ_BINARY(valueObjLS, Tn_obj);
            READ_BINARY(valueObjLS, hflux_z_obj);
            READ_BINARY(valueObjLS, hfluxn_z_obj);
            READ_BINARY(valueObjLS, bcTypes_f);
            READ_BINARY(valueObjLS, dirichlet_weight);
            READ_BINARY(valueObjLS, sc_scalar, n_scalars);
            #undef READ_BINARY
            // valueLBM
            for(int q=0; q<LBM_velocity_model::nQ; q++) {
                read( "lbm_" + std::to_string(q), 
                        reinterpret_cast<char*>(meshValue.valueLBM().f_lbm() + q*nn_max),
                        sizeof(real)
                );
            }
            // valueTimeAverage
            read("u_mean", reinterpret_cast<char*>(valueTimeAverage[lv].u_mean()), sizeof(real));
            read("v_mean", reinterpret_cast<char*>(valueTimeAverage[lv].v_mean()), sizeof(real));
            read("w_mean", reinterpret_cast<char*>(valueTimeAverage[lv].w_mean()), sizeof(real));
            read("T_mean", reinterpret_cast<char*>(valueTimeAverage[lv].T_mean()), sizeof(real));
            read("vel2_fluc", reinterpret_cast<char*>(valueTimeAverage[lv].vel2_fluc()), sizeof(real));
        }
    }
}


void  IOData::
writeVTKFile (
    const Grid*          grids,
    const Tree&          tree,
    const Parameters&    parameters,
    const MeshValue*     meshValues,
    const VTKOutputScale vtkOutputScale,
    const int            filter_bits,
    const std::vector<real>& zSlices,
    const std::vector<real>& ySlices,
    const std::vector<real>& xSlices
    )
const
{
#ifdef NO_IODATA_PARAVIEW
    if(comm_.is_rank0()) { std::cout << "skip VTK output (due to NO_IODATA_PARAVIEW defined)" << std::endl; }
#else


#ifdef PARAVIEW_ENS0
    if(comm_.col_id_wo_offset() != 0) { return; } // skip output for ensembles except ens 0
#endif
    const int  ndiv = (vtkOutputScale == VTKOutputScale::Full      ) ? 1 :
                      (vtkOutputScale == VTKOutputScale::Downsize2 ) ? 2 :
                      (vtkOutputScale == VTKOutputScale::Downsize4 ) ? 4 :
                      (vtkOutputScale == VTKOutputScale::Slice_Full) ? 1 :
                                                                       1 ;

    const real c_ref = parameters.c_ref_lbm();

//    std::cout << __PRETTY_FUNCTION__ << std::endl;

    const int  n_leaf = tree.number_of_nodes();
    const int  nx = DefAMR::NX_LEAF;
    const int  ny = DefAMR::NX_LEAF;
    const int  nz = DefAMR::NX_LEAF;

//    const int  nn  = nx*ny*nz;              // Cell //
//    const int  nnp = (nx+1)*(ny+1)*(nz+1);  // Node //
    const int  nn_ndiv  = nx/ndiv*ny/ndiv*nz/ndiv;              // Cell //
    const int  nnp_ndiv = (nx/ndiv+1)*(ny/ndiv+1)*(nz/ndiv+1);  // Node //

    // Points & Cells //
    std::vector<uint64_t> index;

    std::vector<float>    xyz_vtu; // node //
    std::vector<int64_t>  connectivity_vtu; // node //
    std::vector<int64_t>  offsets_vtu; // cell //
    std::vector<uint8_t>  types_vtu; // cell //

    // PointData or CellData (decided by PARAVIEW_CELLDATA) //
    std::vector<float>    vel_vtu;
    std::vector<float>    rho_vtu;

    std::vector<float>    lv_obj_vtu;
    std::vector<float>    T_obj_vtu;
    std::vector<float>    vel_obj_vtu;

    std::vector<float>    T_vtu;
    std::vector<float>    hflux_z_obj_vtu;


    // CellData //
    std::vector<uint8_t>  amr_lv_vtu;
    std::vector<uint16_t> rank_vtu;


    int64_t  leaf_count  = 0;
    int64_t  point_count = 0;
    int64_t  cell_count  = 0;
    const int n_scalars = meshValues[0].n_scalars();
    std::vector< std::vector<float> >   scalars_vtu(n_scalars);

    for (int l=0; l<n_leaf; l++) {
        const int  lv = tree.nodes(l)->level();
        const int  nx = DefAMR::NX_LEAF;
        const int  ny = DefAMR::NX_LEAF;
        const int  nz = DefAMR::NX_LEAF;

        const auto& meshValue = meshValues[lv];
        const auto& coordinates = meshValue.coordinates();
        const int nn_max    = meshValue.nn_max();

        const Array3D<int>  offsets     = tree.nodes(l)->neighbor_mesh_offsets();
        const Array3D<int>  offsets_lbm = tree.nodes(l)->neighbor_lbm_mesh_offsets();

        const int  offset     = offsets.offset0();
        const int  offset_lbm = offsets_lbm.offset0();

        // skip leaf cal_flag //
        if ( !tree.nodes(l)->nodeCalFlags().Cal() ) { continue; }

        if( vtkOutputScale == VTKOutputScale::Slice_Full ) {
            // skip leaf zSlice //
            const real z0 = coordinates.znode(offset);
            const real z1 = coordinates.znode(offset, nz);
            bool is_output_z = false;
            for(const real& z_val : zSlices) { is_output_z |= (z0 <= z_val && z_val <= z1); }
            // skip leaf ySlice //
            bool is_output_y = false;
            const real y0 = coordinates.ynode(offset);
            const real y1 = coordinates.ynode(offset, ny);
            for(const real& y_val : ySlices) { is_output_y |= (y0 <= y_val && y_val <= y1); }
            // skip leaf xSlice //
            bool is_output_x = false;
            const real x0 = coordinates.xnode(offset);
            const real x1 = coordinates.xnode(offset, nx);
            for(const real& x_val : xSlices) { is_output_x |= (x0 <= x_val && x_val <= x1); }
            // merge skip leaf xyz //
            bool is_skip = true;
            if(!zSlices.empty()) { is_skip &= (!is_output_z); }
            if(!ySlices.empty()) { is_skip &= (!is_output_y); }
            if(!xSlices.empty()) { is_skip &= (!is_output_x); }
            if(is_skip) { continue; }
        }

        // debug //
        if ( tree.nodes(l)->node_type() == NodeTypes::InvalidT ) { std::cout << __PRETTY_FUNCTION__ << " : NodeTypes::InvalidT leaf = " << l << std::endl; exit(-1); }

        // xyz_vtu //
        for (int k=0; k<nz+1; k=k+ndiv) {
        for (int j=0; j<ny+1; j=j+ndiv) {
        for (int i=0; i<nx+1; i=i+ndiv) {
            xyz_vtu.push_back( coordinates.xnode(offset, i) );
            xyz_vtu.push_back( coordinates.ynode(offset, j) );
            xyz_vtu.push_back( coordinates.znode(offset, k) );

            point_count++;
        }
        }
        }


        // PointData or CellData (decided by PARAVIEW_CELLDATA) //
        #ifdef PARAVIEW_CELLDATA
        for (int k=0; k<nz; k=k+ndiv) {
        for (int j=0; j<ny; j=j+ndiv) {
        for (int i=0; i<nx; i=i+ndiv) {
        #else
        for (int k=0; k<nz+1; k=k+ndiv) {
        for (int j=0; j<ny+1; j=j+ndiv) {
        for (int i=0; i<nx+1; i=i+ndiv) {
        #endif
            const int  nx_leaf = DefAMR::NX_LEAF;

            auto&& set_data = [&](std::vector<float>& v, const real* p, const real& factor=1) -> float {
                union { float val; std::uint32_t bitwise; } u;
                #ifdef PARAVIEW_CELLDATA
                u.val = factor * FuncAMRMesh::RawData(p, i, j, k, nx_leaf, offsets);
                #else
                u.val = factor * FuncAMRMesh::CellToNode( p, i,j,k, nx_leaf, offsets );
                #endif
                u.bitwise &= (0xffffffffu << filter_bits);
                v.push_back(u.val);
                return u.val;
            };

            // vel_vtu //
            set_data(vel_vtu, meshValue.valueNS().u(), c_ref );
            set_data(vel_vtu, meshValue.valueNS().v(), c_ref );
            set_data(vel_vtu, meshValue.valueNS().w(), c_ref );

            // vel_obj_vtu //
            set_data(vel_obj_vtu, meshValue.valueObjLS().u_obj(), c_ref );
            set_data(vel_obj_vtu, meshValue.valueObjLS().v_obj(), c_ref );
            set_data(vel_obj_vtu, meshValue.valueObjLS().w_obj(), c_ref );

            // rho_vtu //
            set_data(rho_vtu, meshValue.valueNS().rho() );

            // lv_obj_vtu //
            set_data(lv_obj_vtu, meshValue.valueObjLS().lv_obj() );
            set_data( T_obj_vtu, meshValue.valueObjLS(). T_obj() );

            // scalar //
            for(int n=0; n<n_scalars; n++) {
                set_data( scalars_vtu.at(n), meshValue.valueNS().scalar (n) );
            }

            // temperature //
            set_data(T_vtu, meshValue.valueNS().T());
            set_data(hflux_z_obj_vtu, meshValue.valueObjLS(). hflux_z_obj() );

            // stat //

        }
        }
        }


        // Cells //
        // connectivity_vtu, offsets_vtu, types_vtu //
        for (int k=0; k<nz; k=k+ndiv) {
        for (int j=0; j<ny; j=j+ndiv) {
        for (int i=0; i<nx; i=i+ndiv) {
            // connectivity_vtu //
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv  , j/ndiv  , k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv+1, j/ndiv  , k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv  , j/ndiv+1, k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv+1, j/ndiv+1, k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv  , j/ndiv  , k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv+1, j/ndiv  , k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv  , j/ndiv+1, k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*leaf_count + FuncLoop::id(i/ndiv+1, j/ndiv+1, k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );

            // offsets_vtu //
//            offsets_vtu.push_back( (cell_count + 1)*8 );
            offsets_vtu.push_back( (nn_ndiv*leaf_count + FuncLoop::id(i/ndiv,j/ndiv,k/ndiv, nx/ndiv,ny/ndiv,nz/ndiv) + 1)*8 );

            // types_vtu //
            types_vtu.push_back( 11 ); // 11 : cube //

            cell_count++;
        }
        }
        }

        // CellData //
        for (int k=0; k<nz; k=k+ndiv) {
        for (int j=0; j<ny; j=j+ndiv) {
        for (int i=0; i<nx; i=i+ndiv) {
            amr_lv_vtu.push_back( lv );
            rank_vtu.push_back( comm_.col_vector().rank() );
        }
        }
        }

        // indices //
        leaf_count++;
    }


    ParaviewVTU  paraviewVTU(comm_, vtkOutputScale);
    paraviewVTU.set_number_of_points_and_cells(point_count, cell_count);
    paraviewVTU.set_Point_and_Cells( xyz_vtu.data(), connectivity_vtu.data(), offsets_vtu.data(), types_vtu.data() );

#ifdef PARAVIEW_CELLDATA
#define PUSH_BACK_PHYS push_back_CellData
#else
#define PUSH_BACK_PHYS push_back_PointData
#endif

    paraviewVTU.PUSH_BACK_PHYS("vel",      (BYTE*)vel_vtu.data(),       3, VTKDataType::Float32);
    paraviewVTU.PUSH_BACK_PHYS("rho",      (BYTE*)rho_vtu.data(),       1, VTKDataType::Float32);
    paraviewVTU.PUSH_BACK_PHYS("lv_obj",   (BYTE*)lv_obj_vtu.data(),    1, VTKDataType::Float32);
    paraviewVTU.PUSH_BACK_PHYS("T_obj",    (BYTE*)T_obj_vtu.data(),     1, VTKDataType::Float32);
    paraviewVTU.PUSH_BACK_PHYS("vel_obj",  (BYTE*)vel_obj_vtu.data(),   3, VTKDataType::Float32);
    for(int n=0; n<n_scalars; n++) {
        paraviewVTU.PUSH_BACK_PHYS("scalar"+std::to_string(n),   (BYTE*)scalars_vtu.at(n).data (),   1, VTKDataType::Float32);
    }
    paraviewVTU.PUSH_BACK_PHYS("T",           (BYTE*)T_vtu.data(),           1, VTKDataType::Float32);
    paraviewVTU.PUSH_BACK_PHYS("hflux_z_obj", (BYTE*)hflux_z_obj_vtu.data(), 1, VTKDataType::Float32);

    paraviewVTU.push_back_CellData("amr_lv", (BYTE*)amr_lv_vtu.data(), 1,  VTKDataType::UInt8);
    paraviewVTU.push_back_CellData("rank",   (BYTE*)rank_vtu.data(),   1,  VTKDataType::UInt16);

    paraviewVTU.OutputVTUFiles(step_);
    if (comm_.col_vector().rank() == 0) { paraviewVTU.OutputPVTUFiles(step_); }

#endif
}


//// private //
//stHDF5Values  IOData::
//makeHDF5Values (
//    const int           nx,
//    const int           ny,
//    const int           nz,
//    const Array3D<int>& offsets,
//    const Array3D<int>& offsets_lbm,
//    const Node*         node,
//    const MeshValue&    meshValue
//    )
//const
//{
//    const int  nn  = nx*ny*nz;              // Cell //
//    const int  nnp = (nx+1)*(ny+1)*(nz+1);  // Node //
//
//    const int  offset     = offsets.offset0();
//    const int  offset_lbm = offsets_lbm.offset0();
//
//    // coordinates xyz //
//    std::unique_ptr<float[]>  xp(new float[nnp]);
//    std::unique_ptr<float[]>  yp(new float[nnp]);
//    std::unique_ptr<float[]>  zp(new float[nnp]);
//
//    // values //
//    std::unique_ptr<float[]>  lv_obj (new float[nnp]);
//    std::unique_ptr<float[]>  rho_obj(new float[nnp]);
//    std::unique_ptr<float[]>  u_obj  (new float[nnp]);
//    std::unique_ptr<float[]>  v_obj  (new float[nnp]);
//    std::unique_ptr<float[]>  w_obj  (new float[nnp]);
//
//    std::unique_ptr<float[]>  u  (new float[nnp]);
//    std::unique_ptr<float[]>  v  (new float[nnp]);
//    std::unique_ptr<float[]>  w  (new float[nnp]);
//    std::unique_ptr<float[]>  rho(new float[nnp]);
//
//    std::unique_ptr<float[]>  scalar(new float[nnp]);
//
//    // lbm values //
//    const int  nQ  = LBM_velocity_model::nQ;
//    std::unique_ptr<float[]>  f_lbm(new float[nn*nQ]);
//
//    // HDF5 grid information //
//    std::unique_ptr<int[]>  amr_lv(new int[1]);
////    std::unique_ptr<int[]>  amr_posx(new int[1]);
////    std::unique_ptr<int[]>  amr_posy(new int[1]);
////    std::unique_ptr<int[]>  amr_posz(new int[1]);
//    std::unique_ptr<int[]>  mpi_rank(new int[1]);
//
//
//    // copy //
//    const real  dx0 = meshValue.coordinates().x(offset + 1    ) - meshValue.coordinates().x(offset);
//    const real  dy0 = meshValue.coordinates().y(offset + nx   ) - meshValue.coordinates().y(offset);
//    const real  dz0 = meshValue.coordinates().z(offset + nx*ny) - meshValue.coordinates().z(offset);
//
//
//    // Loop3d 0-nx+1 //
//    FuncLoop::Loop3d  loop3dp( 0,nx+1, 0,ny+1, 0,nz+1 );
//
//    // coordinates xyz @ node //
//    loop3dp.for_each( [&, offset, dx0, dy0, dz0](const int i, const int j, const int k)
//        {
//            const int  idp = FuncLoop::id(i,j,k, nx+1,ny+1,nz+1);
//
//            xp[idp] = meshValue.coordinates().x(offset) + dx0*i;
//            yp[idp] = meshValue.coordinates().y(offset) + dy0*j;
//            zp[idp] = meshValue.coordinates().z(offset) + dz0*k;
//        } );
//
//
//    // values //
//    loop3dp.for_each( [&, offsets](const int i, const int j, const int k)
//        {
//            const int  idp = FuncLoop::id(i,j,k, nx+1,ny+1,nz+1);
//
//            const int  nx_leaf = DefAMR::NX_LEAF;
////            const int  ir  = (i < nx) ? i : nx - 1;
////            const int  jr  = (j < ny) ? j : ny - 1;
////            const int  kr  = (k < nz) ? k : nz - 1;
//            const int  ir  = i;
//            const int  jr  = j;
//            const int  kr  = k;
//
//            // NodeToCell //
//            lv_obj [idp] = FuncAMRMesh::RawData( meshValue.valueObjLS().lv_obj(),  ir,jr,kr, nx_leaf, offsets );
//            rho_obj[idp] = FuncAMRMesh::RawData( meshValue.valueObjLS().rho_obj(), ir,jr,kr, nx_leaf, offsets );
//            u_obj  [idp] = FuncAMRMesh::RawData( meshValue.valueObjLS().u_obj(),   ir,jr,kr, nx_leaf, offsets );
//            v_obj  [idp] = FuncAMRMesh::RawData( meshValue.valueObjLS().v_obj(),   ir,jr,kr, nx_leaf, offsets );
//            w_obj  [idp] = FuncAMRMesh::RawData( meshValue.valueObjLS().w_obj(),   ir,jr,kr, nx_leaf, offsets );
//
//            u  [idp]    = FuncAMRMesh::RawData( meshValue.valueNS().u  (), ir,jr,kr, nx_leaf, offsets );
//            v  [idp]    = FuncAMRMesh::RawData( meshValue.valueNS().v  (), ir,jr,kr, nx_leaf, offsets );
//            w  [idp]    = FuncAMRMesh::RawData( meshValue.valueNS().w  (), ir,jr,kr, nx_leaf, offsets );
//            rho[idp]    = FuncAMRMesh::RawData( meshValue.valueNS().rho(), ir,jr,kr, nx_leaf, offsets );
//
//            scalar[idp]    = FuncAMRMesh::RawData( meshValue.valueNS().scalar(), ir,jr,kr, nx_leaf, offsets );
//        } );
//
//    // lbm values //
//    FuncMath::copy_i_omp( f_lbm.get(), &meshValue.valueLBM().f_lbm()[offset_lbm], nn*nQ );
//
//
//    // HDF5 grid information //
//    amr_lv  [0] = node->level();
////    amr_posx[0] = node->x();
////    amr_posy[0] = node->y();
////    amr_posz[0] = node->z();
//    mpi_rank[0] = rank_;
//
//
//    return  stHDF5Values {
//                // coordinates xyz //
//                std::move(xp),
//                std::move(yp),
//                std::move(zp),
//                // values //
//                std::move(lv_obj),
//                std::move(rho_obj),
//                std::move(u_obj),
//                std::move(v_obj),
//                std::move(w_obj),
//                std::move(u),
//                std::move(v),
//                std::move(w),
//                std::move(rho),
//                std::move(scalar),
//                // lbm //
//                std::move(f_lbm),
//                // HDF5 grid information //
//                std::move(amr_lv),
////                std::move(amr_posx),
////                std::move(amr_posy),
////                std::move(amr_posz),
//                std::move(mpi_rank)
//                };
//}
