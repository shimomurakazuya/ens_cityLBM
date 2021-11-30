#pragma once
#ifndef PARAVIEWVTU_H_
#define PARAVIEWVTU_H_


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>
#include "definePrecision.h"
#include "defineFilenames.h"
#include "MPICommEnsemble.h"
#include "mpi_wrapper.hpp"

enum VTKDataType {
    UInt64,
    UInt32,
    UInt16,
    UInt8,
    Int64,
    Int32,
    Int16,
    Float32,
};


enum VTKOutputScale {
    Full,
    Downsize2,
    Downsize4,
    Slice_Full
};


class  ParaviewVTU {
private:
    struct DataArray {
        std::string plot_name;
        BYTE*       val;
        int         num_components;
        VTKDataType dataType;


        DataArray(std::string _plot_name, BYTE* _val, int _num_components, VTKDataType _dataType):
            plot_name(_plot_name),
            val(_val),
            num_components(_num_components),
            dataType(_dataType){}
    };

private:
    int  number_of_points_ = 0;
    int  number_of_cells_  = 0;

    std::vector<BYTE>  data_;
    int64_t     offset_ = 0;

    float*      xyz_;
    int64_t*    connectivity_;
    int64_t*    offsets_;
    uint8_t*    types_;

    std::vector<DataArray>  pointData_;
    std::vector<DataArray>  cellData_;

    const VTKDataType  header_type_ = UInt64;
    const std::string  format_ = "appended";

    const MPICommEnsemble comm_;
    const util::mpi::int_type rank_;
    const VTKOutputScale  vtkOutputScale_;

public:
    ParaviewVTU() = delete;

    ParaviewVTU(
    const MPICommEnsemble comm, 
    VTKOutputScale vtkOutputScale = VTKOutputScale::Full
    ) : 
    comm_(comm), 
    rank_(comm.world().rank()),
    vtkOutputScale_(vtkOutputScale) 
    { }

    ~ParaviewVTU() {};

public:
    void set_number_of_points_and_cells(int number_of_points, int number_of_cells)
    {
        number_of_points_ = number_of_points;
        number_of_cells_  = number_of_cells;
    }

    void set_Point_and_Cells(
        float*    xyz,
        int64_t*  connectivity,
        int64_t*  offsets,
        uint8_t*  types
        )
    {
        xyz_ = xyz;

        connectivity_ = connectivity;
        offsets_      = offsets;
        types_        = types;
    }

    void push_back_PointData(std::string plot_name, BYTE* val, int num_components, VTKDataType datatype)
    {
        DataArray tmp(plot_name, val, num_components, datatype);
        pointData_.emplace_back( tmp );
    }

    void push_back_CellData(std::string plot_name, BYTE* val, int num_components, VTKDataType datatype)
    {
        DataArray tmp(plot_name, val, num_components, datatype);
        cellData_.emplace_back( tmp );
    }

    void OutputVTUFiles(int step);
    void OutputPVTUFiles(int step);

private:
    // vtu //
    void SetDataArray(std::ofstream& fout, std::string name, VTKDataType type, int number_of_components, BYTE* val, uint64_t total_elem);

    void HeaderVTUFile(std::ofstream& fout);
    void FooterVTUFile(std::ofstream& fout);

    void HeaderPieceFile(std::ofstream& fout);
    void FooterPieceFile(std::ofstream& fout);

    void HeaderAppendedData(std::ofstream& fout);
    void FooterAppendedData(std::ofstream& fout);

    void HeaderData(std::ofstream& fout, std::string str);
    void FooterData(std::ofstream& fout, std::string str);

    void WriteAppendedData(std::ofstream& fout);

    // pvtu //
    void SetPDataArray(std::ofstream& fout, std::string name, VTKDataType type, int number_of_components);

    void HeaderPVTUFile(std::ofstream& fout);
    void SetPieceSource(std::ofstream& fout, int step, int nprocs, int iproc0); // iproc0: head of comm.col_vector()

private:
    std::string cout_type(VTKDataType datatype)
    {
        return  ( datatype == VTKDataType::UInt64  ) ? "UInt64" :
                ( datatype == VTKDataType::UInt32  ) ? "UInt32" :
                ( datatype == VTKDataType::UInt16  ) ? "UInt16" :
                ( datatype == VTKDataType::UInt8   ) ? "UInt8" :
                ( datatype == VTKDataType::Int64   ) ? "Int64"  :
                ( datatype == VTKDataType::Int32   ) ? "Int32"  :
                ( datatype == VTKDataType::Int16   ) ? "Int16"  :
                ( datatype == VTKDataType::Float32 ) ? "Float32"  :
                                                    " ";
    }


    uint64_t get_sizeof(VTKDataType datatype)
    {
        return  ( datatype == VTKDataType::UInt64  ) ? sizeof(uint64_t) :
                ( datatype == VTKDataType::UInt32  ) ? sizeof(uint32_t) :
                ( datatype == VTKDataType::UInt16  ) ? sizeof(uint16_t) :
                ( datatype == VTKDataType::UInt8   ) ? sizeof(uint8_t) :
                ( datatype == VTKDataType::Int64   ) ? sizeof(int64_t)  :
                ( datatype == VTKDataType::Int32   ) ? sizeof(int32_t)  :
                ( datatype == VTKDataType::Int16   ) ? sizeof(int16_t)  :
                ( datatype == VTKDataType::Float32 ) ? sizeof(float)    :
                                                       sizeof(BYTE);
    }

    template<typename T>
    void
    add_appended_data(std::vector<BYTE>& vec, const T a)
    {
        const int m = sizeof(T)/sizeof(BYTE);
        const BYTE* tmp = (BYTE*) &a;

        vec.insert( vec.end(), tmp, tmp + m );
    }


    template<typename T>
    void
    add_appended_data(std::vector<BYTE>& vec, const T* a, const int n)
    {
        const int m = sizeof(T)/sizeof(BYTE) * n;
        const BYTE* tmp = (BYTE*)a;

        vec.insert( vec.end(), tmp, tmp + m );
    }


    std::string io_folder() { return Foldernames::io_folder; }
    std::string to_str_zf(int i, int digit=4) {  return std::string(digit - std::to_string(i).length(), '0') + std::to_string(i); }

    std::string fname_vtu(int step, int rank)
    {
        const auto str_rank = to_str_zf(rank);
        const auto str_step = to_str_zf(step);
        return (vtkOutputScale_ == VTKOutputScale::Full      ) ? Filenames::vtu_value_name0                + "-rank" + str_rank + "-step" + str_step + ".vtu":
               (vtkOutputScale_ == VTKOutputScale::Downsize2 ) ? Filenames::vtu_value_name0 + "-downsize2" + "-rank" + str_rank + "-step" + str_step + ".vtu":
               (vtkOutputScale_ == VTKOutputScale::Downsize4 ) ? Filenames::vtu_value_name0 + "-downsize4" + "-rank" + str_rank + "-step" + str_step + ".vtu":
               (vtkOutputScale_ == VTKOutputScale::Slice_Full) ? Filenames::vtu_value_name0 + "-slice"     + "-rank" + str_rank + "-step" + str_step + ".vtu":
                                                                 Filenames::vtu_value_name0 + "-unknown"   + "-rank" + str_rank + "-step" + str_step + ".vtu";
    }


    std::string fname_pvtu(int step)
    {
        const auto str_ens = to_str_zf(comm_.col_id());
        const auto str_step = to_str_zf(step);
        return (vtkOutputScale_ == VTKOutputScale::Full     )  ? Filenames::pvtu_value_name0                + "-ens" + str_ens + "-step" + str_step + ".pvtu":
               (vtkOutputScale_ == VTKOutputScale::Downsize2)  ? Filenames::pvtu_value_name0 + "-downsize2" + "-ens" + str_ens + "-step" + str_step + ".pvtu":
               (vtkOutputScale_ == VTKOutputScale::Downsize4)  ? Filenames::pvtu_value_name0 + "-downsize4" + "-ens" + str_ens + "-step" + str_step + ".pvtu":
               (vtkOutputScale_ == VTKOutputScale::Slice_Full) ? Filenames::pvtu_value_name0 + "-slice"     + "-ens" + str_ens + "-step" + str_step + ".pvtu":
                                                                 Filenames::pvtu_value_name0 + "-unknown"   + "-ens" + str_ens + "-step" + str_step + ".pvtu";
    }

};


#endif
