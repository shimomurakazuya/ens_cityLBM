#pragma once
#ifndef HDFTYPES_H_
#define HDFTYPES_H_


#include <iostream>
#include <memory>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <vector>


namespace  HDFTypes {


enum  DataTypes  {
    TypeFloat,
    TypeInt
};


struct  NameList {
    std::string     grid_name;
    std::string     folder_name;
    std::string     file_name;
    std::string     attribute_name;
    std::string     index_name;
    std::string     pos_name;


    NameList(
        const std::string   _grid_name,
        const std::string   _folder_name,
        const std::string   _file_name,
        const std::string   _attribute_name,
        const std::string   _index_name,
        const std::string   _pos_name
        ) :
        grid_name(_grid_name),
        folder_name(_folder_name),
        file_name(_file_name),
        attribute_name(_attribute_name),
        index_name(_index_name),
        pos_name(_pos_name)
    {}

public:
    const std::string  dataset_name()  const { return  attribute_name + index_name; }

};


struct  WriteFormat {
public:
    NameList                nameList;
    std::vector<hsize_t>    dims;
    DataTypes               dataType;
    void*                   val;


    WriteFormat(
        NameList                _nameList,
        std::vector<hsize_t>    _dims,
        DataTypes               _dataType,
        void*                   _val
        ) :
        nameList(_nameList),
        dims(_dims),
        dataType(_dataType),
        val(_val)
    {}


public:
    hid_t  type_id() {
        if      (dataType == DataTypes::TypeFloat) { return  H5T_NATIVE_FLOAT; }
        else if (dataType == DataTypes::TypeInt  ) { return  H5T_NATIVE_INT; }
        else                                       { return  (hid_t)NULL; }
        }

    std::string  numberType() {
        if      (dataType == DataTypes::TypeFloat) { return  "Float"; }
        else if (dataType == DataTypes::TypeInt  ) { return  "Int"; }
        else                                       { return  NULL; }
        }

    int  precision() {
        if      (dataType == DataTypes::TypeFloat) { return  4; }
        else if (dataType == DataTypes::TypeInt  ) { return  4; }
        else                                       { return  0; }
        }

};


};


#endif
