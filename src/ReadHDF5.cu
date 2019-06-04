#if 0

#include "ReadHDF5.h"


namespace  ReadHDF5 {


stHDF5Values
read_hdf_datasets(
    const hid_t  file_id,
    const hid_t  file_id_lbm,
    const int    i
    )
{
    constexpr int ndim  = 3;


    // read //
    std::unique_ptr<float[]>  xp( read_hdf5_dataset<float>(ndim, file_id, ("X" + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  yp( read_hdf5_dataset<float>(ndim, file_id, ("Y" + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  zp( read_hdf5_dataset<float>(ndim, file_id, ("Z" + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );

    // values //
    std::unique_ptr<float[]>  lv_obj ( read_hdf5_dataset<float>(ndim, file_id, ("lv_obj"  + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  rho_obj( read_hdf5_dataset<float>(ndim, file_id, ("rho_obj" + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  u_obj  ( read_hdf5_dataset<float>(ndim, file_id, ("u_obj"   + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  v_obj  ( read_hdf5_dataset<float>(ndim, file_id, ("v_obj"   + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  w_obj  ( read_hdf5_dataset<float>(ndim, file_id, ("w_obj"   + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );

    std::unique_ptr<float[]>  u  ( read_hdf5_dataset<float>(ndim, file_id, ("u"   + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  v  ( read_hdf5_dataset<float>(ndim, file_id, ("v"   + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  w  ( read_hdf5_dataset<float>(ndim, file_id, ("w"   + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );
    std::unique_ptr<float[]>  rho( read_hdf5_dataset<float>(ndim, file_id, ("rho" + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );

    std::unique_ptr<float[]>  scalar( read_hdf5_dataset<float>(ndim, file_id, ("scalar" + std::to_string(i)).c_str(), H5T_NATIVE_FLOAT) );

    // lbm values //
    const int  nQ = LBM_velocity_model::nQ;
    std::vector <std::string>  strings_f_lbm;
    for (int iv=0; iv<nQ; iv++) {
        const std::string  string_f_lbm("f" + std::to_string(iv) + "_lbm");
        strings_f_lbm.push_back(string_f_lbm + std::to_string(i));
    }
    std::unique_ptr<float[]>  f_lbm( read_hdf5_dataset_series<float>(ndim, file_id_lbm, strings_f_lbm, H5T_NATIVE_FLOAT) );


    // move //
    stHDF5Values  HDF5value {
                    move(xp),
                    move(yp),
                    move(zp),
                    move(lv_obj),
                    move(rho_obj),
                    move(u_obj),
                    move(v_obj),
                    move(w_obj),
                    move(u),
                    move(v),
                    move(w),
                    move(rho),
                    move(scalar),
                    move(f_lbm)
                    };


    return  HDF5value;
}


};

#endif
