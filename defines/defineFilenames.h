#pragma once
#ifndef DEFINEFILENAMES_H_
#define DEFINEFILENAMES_H_


namespace  Foldernames {
    const std::string  input_folder = "../io/input";
    const std::string  output_folder = "../io/output";
    const std::string  io_folder = "../io/iofiles";
};


namespace  Filenames {
    const std::string  info_name0 = "info";

    const std::string  parameter_name0 = "parameter";
    const std::string  grid_name0      = "grid";
    const std::string  tree_name0      = "tree";

    const std::string  taskID_name0    = "taskID";

    const std::string  hdf5_value_name0 = "hdf5_value";

    const std::string  vtu_value_name0   = "vtu_value";
    const std::string  pvtu_value_name0  = "pvtu_value";

    const std::string  ElapsedTimeInfo_name0 = "elapsedTimeInfo";

//    const std::string  Map_name0 = "map/output10x10.ssv";
//    const std::string  Map_name0 = "map/output300x300.ssv";
//    const std::string  Map_name0 = "map/output800x800.ssv";
//    const std::string  Map_name0 = "map/output1400x1400.ssv";
//    const std::string  Map_name0 = "map/output1600x1600.ssv";
    const std::string  Map_name0 = "map/output2950x2950.ssv";

    const std::string ChannelFlow_statics = "channel_flow_statics";
    const std::string NaturalConvection2d_statics = "natural_convection2d_statics";
    const std::string NaturalConvection3d_statics = "natural_convection3d_statics";
    const std::string TaylorGreen_statics = "taylor_green_statics";
};


#endif
