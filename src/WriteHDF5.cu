#if 0

#include "WriteHDF5.h"
#include "WriteXDMF.h"
#include "defineLBM.h"
#include "defineAMR.h"
#include <sstream>


namespace  WriteHDF5 {


void
write_hdf_datasets (
    const std::string                   folder,
    const std::string                   filename,
    const std::string                   filename_lbm,
    const std::string                   grid_name0,
    const std::vector<stHDF5Values>&    HDF5Values,
    const std::vector<int>&             lv_leaf,
    const std::vector<bool>&            flag_visualization,
    const int                           nx_cell[]
    )
{
    const unsigned long long  nx = nx_cell[0];
    const unsigned long long  ny = nx_cell[1];
    const unsigned long long  nz = nx_cell[2];

//    std::cout << "nx, ny, nz = " << nx << ", " << ny << ", " << nz << std::endl;

    std::vector<WriteFormat>  ioGeometryDataType;
    std::vector<WriteFormat>  ioAttributeDataType;
    std::vector<WriteFormat>  ioLBMAttributeDataType;
    std::vector<WriteFormat>  ioGridAttributeDataType;

    const int  number_of_leaves = HDF5Values.size();
    for (int i=0; i<number_of_leaves; i++) {
        const std::string grid_name = grid_name0 + std::to_string(i);

        ioGeometryDataType.emplace_back ( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "X", (std::to_string(i)).c_str(), "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].xp.get() );
        ioGeometryDataType.emplace_back ( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "Y", (std::to_string(i)).c_str(), "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].yp.get() );
        ioGeometryDataType.emplace_back ( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "Z", (std::to_string(i)).c_str(), "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].zp.get() );


        // values //
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "lv_obj",  (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].lv_obj.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "rho_obj", (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].rho_obj.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "u_obj",   (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].u_obj.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "v_obj",   (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].v_obj.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "w_obj",   (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].w_obj.get() );

        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "u",   (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].u.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "v",   (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].v.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "w",   (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].w.get() );
        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "rho", (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].rho.get() );

        ioAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "scalar", (std::to_string(i)).c_str(),  "Node"}, std::vector<hsize_t>{nz+1, ny+1, nx+1}, DataTypes::TypeFloat, HDF5Values[i].scalar.get() );



        // lbm values //
        const int  nQ = LBM_velocity_model::nQ;
        const int  nn = nx*ny*nz;
        for (int iv=0; iv<nQ; iv++) {
            const std::string  string_f_lbm("f" + std::to_string(iv) + "_lbm");

            ioLBMAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename_lbm.c_str(), string_f_lbm.c_str(), (std::to_string(i)).c_str(),  "Cell"}, std::vector<hsize_t>{nz, ny, nx}, DataTypes::TypeFloat, &HDF5Values[i].f_lbm.get()[nn*iv] );
        }

        // HDF5 grid information //
        ioGridAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "amr_lv",   (std::to_string(i)).c_str(),  "Grid"}, std::vector<hsize_t>{1, 1, 1}, DataTypes::TypeInt, HDF5Values[i].amr_lv.get() );
//        ioGridAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "amr_posx",   (std::to_string(i)).c_str(),  "Grid"}, std::vector<hsize_t>{1, 1, 1}, DataTypes::TypeInt, HDF5Values[i].amr_posx.get() );
//        ioGridAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "amr_posy",   (std::to_string(i)).c_str(),  "Grid"}, std::vector<hsize_t>{1, 1, 1}, DataTypes::TypeInt, HDF5Values[i].amr_posy.get() );
//        ioGridAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "amr_posz",   (std::to_string(i)).c_str(),  "Grid"}, std::vector<hsize_t>{1, 1, 1}, DataTypes::TypeInt, HDF5Values[i].amr_posz.get() );
        ioGridAttributeDataType.emplace_back( NameList{grid_name.c_str(), folder.c_str(), filename.c_str(), "mpi_rank", (std::to_string(i)).c_str(),  "Grid"}, std::vector<hsize_t>{1, 1, 1}, DataTypes::TypeInt, HDF5Values[i].mpi_rank.get() );

    }


    // Write the data file.
    write_hdf5_and_xdmf_datasets (folder, filename, filename_lbm, ioGeometryDataType, ioAttributeDataType, ioLBMAttributeDataType, ioGridAttributeDataType, lv_leaf, flag_visualization);
}


void
write_hdf_dataset_impl (
    const hid_t&  file_id,
    const std::vector<WriteFormat>& vec
    )
{
//    for (auto elem : vec) {
//        H5LTmake_dataset(file_id,  ("/" + elem.nameList.dataset_name()).c_str(), elem.dims.size(), elem.dims.data(), elem.type_id(), elem.val);
//    }

    for (int i=0; i<vec.size(); i++) {
        auto elem = vec[i];
        H5LTmake_dataset(file_id,  ("/" + elem.nameList.dataset_name()).c_str(), elem.dims.size(), elem.dims.data(), elem.type_id(), elem.val);
    }
}


void
write_hdf5_and_xdmf_datasets (
    const std::string               folder,
    const std::string               filename,
    const std::string               filename_lbm,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioLBMAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType,
    const std::vector<int>&         lv_leaf,
    const std::vector<bool>&        flag_visualization
    )
{
    /// xdmf ///
    // full //
    write_xdmf_dataset(folder, filename, ioGeometryDataType, ioAttributeDataType, ioLBMAttributeDataType, ioGridAttributeDataType);

    // lv_leaf //
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        std::string  lv_filename = "lv" + std::to_string(lv) + "_" + filename;
        write_xdmf_dataset_for_lv(folder, lv_filename, ioGeometryDataType, ioAttributeDataType, ioGridAttributeDataType, lv_leaf, lv);
    }

    // visulalization //
    std::string  vis_filename = "vis_" + filename;
    write_xdmf_dataset_for_visualization(folder, vis_filename, ioGeometryDataType, ioAttributeDataType, ioGridAttributeDataType, flag_visualization);


    /// hdf5 ///
    write_hdf5_dataset(folder, filename, filename_lbm, ioGeometryDataType, ioAttributeDataType, ioLBMAttributeDataType, ioGridAttributeDataType);
}


void
write_hdf5_dataset (
    const std::string               folder,
    const std::string               filename,
    const std::string               filename_lbm,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioLBMAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType
    )
{
    write_hdf5_dataset_for_visualization (
        folder,
        filename,
        ioGeometryDataType,
        ioAttributeDataType,
        ioGridAttributeDataType
        );


    write_hdf5_dataset_lbm (
        folder,
        filename_lbm,
        ioLBMAttributeDataType
        );

//    write_hdf5_dataset_lbm (
//        folder,
//        filename_lbm + "_",
//        ioLBMAttributeDataType
//        );
}


void
write_hdf5_dataset_lbm (
    const std::string               folder,
    const std::string               filename_lbm,
    const std::vector<WriteFormat>& ioLBMAttributeDataType
    )
{
    // hdf5 lbm //
    static hid_t     file_id_lbm;
    file_id_lbm = H5Fcreate((folder + "/" + filename_lbm + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    write_hdf_dataset_impl (file_id_lbm, ioLBMAttributeDataType);

    H5Fclose(file_id_lbm);
}


void
write_hdf5_dataset_for_visualization (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType
    )
{
    // hdf5 //
    static hid_t     file_id;
    file_id = H5Fcreate((folder + "/" + filename + ".h5").c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    write_hdf_dataset_impl (file_id, ioGeometryDataType);
    write_hdf_dataset_impl (file_id, ioAttributeDataType);
    write_hdf_dataset_impl (file_id, ioGridAttributeDataType);

    H5Fclose(file_id);
}


void
write_xdmf_dataset (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioLBMAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType
    )
{
    // xdmf //
    static FILE *xmf = 0;
    xmf = fopen((folder + "/" + filename + ".xmf").c_str(), "w");

    WriteXDMF::write_xdmf_header (xmf);
    WriteXDMF::write_xdmf_datas  (xmf, ioGeometryDataType, ioAttributeDataType, ioLBMAttributeDataType, ioGridAttributeDataType);
    WriteXDMF::write_xdmf_footer (xmf);

    fclose(xmf);
}


void
write_xdmf_dataset_for_visualization (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType,
    const std::vector<bool>&        flag_visualization
    )
{
    // xdmf //
    static FILE *xmf = 0;
    xmf = fopen((folder + "/" + filename + ".xmf").c_str(), "w");

    WriteXDMF::write_xdmf_header (xmf);
    WriteXDMF::write_xdmf_datas_for_visualization (xmf, ioGeometryDataType, ioAttributeDataType, ioGridAttributeDataType, flag_visualization);
    WriteXDMF::write_xdmf_footer (xmf);

    fclose(xmf);
}


void
write_xdmf_dataset_for_lv (
    const std::string               folder,
    const std::string               filename,
    const std::vector<WriteFormat>& ioGeometryDataType,
    const std::vector<WriteFormat>& ioAttributeDataType,
    const std::vector<WriteFormat>& ioGridAttributeDataType,
    const std::vector<int>&         lv_leaf,
    const int                       lv
    )
{
    // xdmf //
    static FILE *xmf = 0;
    xmf = fopen((folder + "/" + filename + ".xmf").c_str(), "w");

    WriteXDMF::write_xdmf_header (xmf);
    WriteXDMF::write_xdmf_datas_for_lv (xmf, ioGeometryDataType, ioAttributeDataType, ioGridAttributeDataType, lv_leaf, lv);
    WriteXDMF::write_xdmf_footer (xmf);

    fclose(xmf);
}


};

#endif
