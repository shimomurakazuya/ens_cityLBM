#include "IOData.h"
#include "defineFilenames.h"
#include "defineAMR.h"
#include "FuncLoop.h"
#include "FuncMath.h"
#include "FuncAMRMesh.h"
//#include "WriteHDF5.h"
//#include "ReadHDF5.h"


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


void  IOData::
writeVTKFile (
    const Grid*          grids,
    const Tree&          tree,
    const Parameters&    parameters,
    const MeshValue*     meshValues,
    const VTKOutputScale vtkOutputScale
    )
const
{
    const int  ndiv = (vtkOutputScale == VTKOutputScale::Full     ) ? 1 :
                      (vtkOutputScale == VTKOutputScale::Downsize2) ? 2 :
                      (vtkOutputScale == VTKOutputScale::Downsize4) ? 4 :
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

    // PointData //
    std::vector<float>    vel_vtu;
    std::vector<float>    rho_vtu;

    std::vector<float>    lv_obj_vtu;
    std::vector<float>    T_obj_vtu;
    std::vector<float>    vel_obj_vtu;

    std::vector<float>    scalar_vtu;
    std::vector<float>    T_vtu;

    // CellData //
    std::vector<uint8_t>  amr_lv_vtu;
    std::vector<uint16_t> rank_vtu;


    int64_t  leaf_count  = 0;
    int64_t  point_count = 0;
    int64_t  cell_count  = 0;

    for (int l=0; l<n_leaf; l++) {
        const int  lv = tree.nodes(l)->level();
        const int  nx = DefAMR::NX_LEAF;
        const int  ny = DefAMR::NX_LEAF;
        const int  nz = DefAMR::NX_LEAF;

        const MeshValue& meshValue = meshValues[lv];

        const Array3D<int>  offsets     = tree.nodes(l)->neighbor_mesh_offsets();
        const Array3D<int>  offsets_lbm = tree.nodes(l)->neighbor_lbm_mesh_offsets();

        const int  offset     = offsets.offset0();
        const int  offset_lbm = offsets_lbm.offset0();

        // coordinates //
        const real  dx0 = (meshValue.coordinates().x(offset + 1    ) - meshValue.coordinates().x(offset));
        const real  dy0 = (meshValue.coordinates().y(offset + nx   ) - meshValue.coordinates().y(offset));
        const real  dz0 = (meshValue.coordinates().z(offset + nx*ny) - meshValue.coordinates().z(offset));


        // debug //
        if ( tree.nodes(l)->node_type() == NodeTypes::InvalidT ) { std::cout << __PRETTY_FUNCTION__ << " : NodeTypes::InvalidT leaf = " << l << std::endl; exit(-1); }

        // xyz_vtu //
        for (int k=0; k<nz+1; k=k+ndiv) {
        for (int j=0; j<ny+1; j=j+ndiv) {
        for (int i=0; i<nx+1; i=i+ndiv) {
            xyz_vtu.push_back( meshValue.coordinates().x(offset) + dx0*i );
            xyz_vtu.push_back( meshValue.coordinates().y(offset) + dy0*j );
            xyz_vtu.push_back( meshValue.coordinates().z(offset) + dz0*k );

            point_count++;
        }
        }
        }


        // PointData //
        for (int k=0; k<nz+1; k=k+ndiv) {
        for (int j=0; j<ny+1; j=j+ndiv) {
        for (int i=0; i<nx+1; i=i+ndiv) {
            const int  nx_leaf = DefAMR::NX_LEAF;

#if 0
            // vel_vtu //
            vel_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueNS().u(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueNS().v(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueNS().w(), i,j,k, nx_leaf, offsets ) * c_ref );

            // rho_vtu //
            rho_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueNS().rho(), i,j,k, nx_leaf, offsets ) );

            // lv_obj_vtu //
            lv_obj_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueObjLS().lv_obj(), i,j,k, nx_leaf, offsets ) );

            // vel_obj_vtu //
            vel_obj_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueObjLS().u_obj(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_obj_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueObjLS().v_obj(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_obj_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueObjLS().w_obj(), i,j,k, nx_leaf, offsets ) * c_ref );

            // scalar //
            scalar_vtu.push_back( FuncAMRMesh::RawData( meshValue.valueNS().scalar(), i,j,k, nx_leaf, offsets ) );
#else
            // vel_vtu //
            vel_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueNS().u(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueNS().v(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueNS().w(), i,j,k, nx_leaf, offsets ) * c_ref );

            // rho_vtu //
            rho_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueNS().rho(), i,j,k, nx_leaf, offsets ) );

            // lv_obj_vtu //
            lv_obj_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueObjLS().lv_obj(), i,j,k, nx_leaf, offsets ) );
            T_obj_vtu. push_back( FuncAMRMesh::CellToNode( meshValue.valueObjLS().T_obj(),  i,j,k, nx_leaf, offsets ) );

            // vel_obj_vtu //
            vel_obj_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueObjLS().u_obj(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_obj_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueObjLS().v_obj(), i,j,k, nx_leaf, offsets ) * c_ref );
            vel_obj_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueObjLS().w_obj(), i,j,k, nx_leaf, offsets ) * c_ref );

            // scalar //
            scalar_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueNS().scalar(), i,j,k, nx_leaf, offsets ) );

            T_vtu.push_back( FuncAMRMesh::CellToNode( meshValue.valueNS().T(), i,j,k, nx_leaf, offsets ) );
#endif
        }
        }
        }

        if ( !tree.nodes(l)->nodeCalFlags().Cal() ) { continue; }

        // Cells //
        // connectivity_vtu, offsets_vtu, types_vtu //
        for (int k=0; k<nz; k=k+ndiv) {
        for (int j=0; j<ny; j=j+ndiv) {
        for (int i=0; i<nx; i=i+ndiv) {
            // connectivity_vtu //
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv  , j/ndiv  , k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv+1, j/ndiv  , k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv  , j/ndiv+1, k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv+1, j/ndiv+1, k/ndiv  , nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv  , j/ndiv  , k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv+1, j/ndiv  , k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv  , j/ndiv+1, k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );
            connectivity_vtu.push_back( nnp_ndiv*l + FuncLoop::id(i/ndiv+1, j/ndiv+1, k/ndiv+1, nx/ndiv+1,ny/ndiv+1,nz/ndiv+1) );

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
            rank_vtu.push_back( rank_ );
        }
        }
        }


        // indices //
        leaf_count++;
    }


    ParaviewVTU  paraviewVTU(vtkOutputScale);
    paraviewVTU.set_number_of_points_and_cells(point_count, cell_count);
    paraviewVTU.set_Point_and_Cells( xyz_vtu.data(), connectivity_vtu.data(), offsets_vtu.data(), types_vtu.data() );

    paraviewVTU.push_back_PointData("vel",      (BYTE*)vel_vtu.data(),       3, VTKDataType::Float32);
    paraviewVTU.push_back_PointData("rho",      (BYTE*)rho_vtu.data(),       1, VTKDataType::Float32);
    paraviewVTU.push_back_PointData("lv_obj",   (BYTE*)lv_obj_vtu.data(),    1, VTKDataType::Float32);
    paraviewVTU.push_back_PointData("T_obj",    (BYTE*)T_obj_vtu.data(),     1, VTKDataType::Float32);
    paraviewVTU.push_back_PointData("vel_obj",  (BYTE*)vel_obj_vtu.data(),   3, VTKDataType::Float32);
    paraviewVTU.push_back_PointData("scalar",   (BYTE*)scalar_vtu.data(),    1, VTKDataType::Float32);
    paraviewVTU.push_back_PointData("T",        (BYTE*)T_vtu.data(),         1, VTKDataType::Float32);

    paraviewVTU.push_back_CellData("amr_lv", (BYTE*)amr_lv_vtu.data(), 1,  VTKDataType::UInt8);
    paraviewVTU.push_back_CellData("rank",   (BYTE*)rank_vtu.data(),   1,  VTKDataType::UInt16);

    paraviewVTU.OutputVTUFiles(step_);
    if (rank_ == 0) { paraviewVTU.OutputPVTUFiles(step_); }

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
