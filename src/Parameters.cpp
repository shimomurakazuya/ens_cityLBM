#include "Parameters.h"
#include "defineAMR.h"
#include "defineFilenames.h"
#include <fstream>
#include <cmath>

#include "bytes_type.hpp"


void  Parameters::
init(const OptionParser& optionParser)
{
    initGPUDevice           (optionParser);

    initCoefFluidPropertyLBM(optionParser);
    initCoefTime            (optionParser);
    initCoefCFROutput       (optionParser);

    initCoefDomain          (optionParser);
}

void Parameters::
coutParameters() const
{
    auto to_hex = [](real d) { return *reinterpret_cast<util::bytes_type<real>*>(&d); };
    std::cout << std::hex << std::showbase;
    std::cout << "test struct Paramters" << std::endl;
    auto&& dom = coefDomain();
    std::cout << " real3 global_domain_min = " 
        << dom.x_global_domain_min << " " 
        << dom.y_global_domain_min << " " 
        << dom.z_global_domain_min << std::endl;
    std::cout << " real3-hex global_domain_min = " 
        << to_hex(dom.x_global_domain_min) << " " 
        << to_hex(dom.y_global_domain_min) << " " 
        << to_hex(dom.z_global_domain_min) << std::endl;
    std::cout << " real3 global_domain_max = " 
        << dom.x_global_domain_max << " " 
        << dom.y_global_domain_max << " " 
        << dom.z_global_domain_max << std::endl;
    std::cout << " real3-hex global_domain_max = " 
        << to_hex(dom.x_global_domain_max) << " " 
        << to_hex(dom.y_global_domain_max) << " " 
        << to_hex(dom.z_global_domain_max) << std::endl;
    std::cout << " real3 global_domain_center() = " 
        << dom.x_global_domain_center() << " " 
        << dom.y_global_domain_center() << " " 
        << dom.z_global_domain_center() << std::endl;
    std::cout << " real3-hex global_domain_center() = " 
        << to_hex(dom.x_global_domain_center()) << " " 
        << to_hex(dom.y_global_domain_center()) << " " 
        << to_hex(dom.z_global_domain_center()) << std::endl;
    std::cout << std::dec;
}


void  Parameters::
copy(const Parameters& parameters)
{
    copyGPUDevice           (parameters.gpuDevice());
    copyCoefFluidPropertyLBM(parameters.coefFluidPropertyLBM());
    copyCoefTime            (parameters.coefTime());
    copyCoefCFROutput       (parameters.coefCFROutput());
    copyCoefDomain          (parameters.coefDomain());
}


void  Parameters::
readParameters(const OptionParser&  optionParser, const std::string  filename)
{
    init(optionParser);

    stIOStep  stIOstep;

    // read //
    std::ifstream  fin;
    fin.open(filename, std::ios::binary);
    fin.read( ( char * ) &stIOstep, sizeof( stIOStep ) );
    fin.close();


    coefTime_.step_now  = stIOstep.step;
    coefTime_.step_init = stIOstep.step;
    coefTime_.time_now  = stIOstep.time;
    coefTime_.time_init = stIOstep.time;

//    std::cout << coefTime_.step_now  << std::endl;
//    std::cout << coefTime_.step_init << std::endl;
//    std::cout << coefTime_.time_now  << std::endl;
//    std::cout << coefTime_.time_init << std::endl;
}


void  Parameters::
writeParameters(const std::string  filename, int fstep)
const
{
    const int   step = coefTime_.step_now;
    const real  time = coefTime_.time_now;
    std::cout << __PRETTY_FUNCTION__
              << " : time, step, fstep = "
              << time << ", " << step << ", " << fstep << std::endl;

    stIOStep  stIOstep{ step, time, fstep };

    // write //
    std::ofstream  fout;
    fout.open(filename, std::ios::binary);
    fout.write( ( char * ) &stIOstep, sizeof( stIOStep ) );
    fout.close();
}


std::string  Parameters::
filename(const int step)
const
{
    return    Foldernames::io_folder + "/"
            + Filenames::parameter_name0
            + "-step" + std::to_string(step)
            + ".dat";
}


// private //
void  Parameters::
initGPUDevice(const OptionParser& optionParser)
{
    const auto gpu_per_node = optionParser.gpu_per_node();
    const auto rank = comm_.world().rank();
    gpuDevice_.gpu_per_node = gpu_per_node;
    gpuDevice_.id_gpu       = rank % gpu_per_node;

#ifdef USE_NVCC
    cudaSetDevice(gpuDevice_.id_gpu);
#endif

    std::cout << "rank, device = " << rank << ", " << gpuDevice_.id_gpu << std::endl;
}


void  Parameters::
initCoefFluidPropertyLBM(const OptionParser& optionParser)
{
    const real  vel_ref_lbm = optionParser.velocity_lbm(0);
    const real  cfl_ref_lbm = optionParser.velocity_lbm(1);

    const real  c_ref = vel_ref_lbm / cfl_ref_lbm;

//    coefFluidPropertyLBM_.vel_ref_lbm  = vel_ref_lbm;
//    coefFluidPropertyLBM_.cfl_ref_lbm  = cfl_ref_lbm;
    coefFluidPropertyLBM_.c_ref_lbm    = c_ref;
}


void  Parameters::
copyGPUDevice(const GPUDevice& gpuDevice)
{
    gpuDevice_.id_gpu       = gpuDevice.id_gpu;
    gpuDevice_.gpu_per_node = gpuDevice.gpu_per_node;
}


void  Parameters::
copyCoefFluidPropertyLBM(const CoefFluidPropertyLBM& coefFluidPropertyLBM)
{
    coefFluidPropertyLBM_.c_ref_lbm = coefFluidPropertyLBM.c_ref_lbm;
}


void  Parameters::
initCoefTime(const OptionParser& optionParser)
{
    const real  time_init = 0.0;
    const real  _time_end = optionParser.time_end();


    // dx, dt //
    const int  nx_tmp  = optionParser.number_of_grid_point(0);
    const int  dim_dir = optionParser.number_of_grid_point(1);
    const real domain_length = static_cast<real>( optionParser.domain_length(dim_dir) );
    const real dx_leaf_lv0 = domain_length / nx_tmp;
    const real dx_lv0 = dx_leaf_lv0 / DefAMR::NX_LEAF;

    const real dt0 = dx_lv0 / c_ref_lbm();
    const int  step = static_cast<int>( ceil( _time_end/dt0 ) );
    const real time_end = step*dt0;


    coefTime_.time_now  = time_init;
    coefTime_.time_init = time_init;
    coefTime_.time_end  = time_end;
    coefTime_.dt0       = dt0;

    coefTime_.step_now  = time_init/dt0;
    coefTime_.step_init = time_init/dt0;
    coefTime_.step_end  = time_end/dt0;
}


void  Parameters::
copyCoefTime(const CoefTime& coefTime)
{
    coefTime_.time_now  = coefTime.time_now;
    coefTime_.time_init = coefTime.time_init;
    coefTime_.time_end  = coefTime.time_end;
    coefTime_.dt0       = coefTime.dt0;

    coefTime_.step_now  = coefTime.step_now ;
    coefTime_.step_init = coefTime.step_init;
    coefTime_.step_end  = coefTime.step_end;
}


void  Parameters::
initCoefCFROutput(const OptionParser& optionParser)
{
    coefCFROutput_.cout_step        = optionParser.cfrm_steps(0);
    coefCFROutput_.fout_step        = optionParser.cfrm_steps(1);
    coefCFROutput_.rout_step        = optionParser.cfrm_steps(2);
    coefCFROutput_.mout_step        = optionParser.cfrm_steps(3);

    coefCFROutput_.cout_flag        = optionParser.cfrm_flags(0);
    coefCFROutput_.fout_flag        = optionParser.cfrm_flags(1);
    coefCFROutput_.rout_flag        = optionParser.cfrm_flags(2);
    coefCFROutput_.mout_flag        = optionParser.cfrm_flags(3);
}


void  Parameters::
copyCoefCFROutput(const CoefCFROutput& coefCFROutput)
{
    coefCFROutput_.cout_step        = coefCFROutput.cout_step;
    coefCFROutput_.fout_step        = coefCFROutput.fout_step;
    coefCFROutput_.rout_step        = coefCFROutput.rout_step;
    coefCFROutput_.mout_step        = coefCFROutput.mout_step;

    coefCFROutput_.cout_flag        = coefCFROutput.cout_flag;
    coefCFROutput_.fout_flag        = coefCFROutput.fout_flag;
    coefCFROutput_.rout_flag        = coefCFROutput.rout_flag;
    coefCFROutput_.mout_flag        = coefCFROutput.mout_flag;
}


void  Parameters::
initCoefDomain(const OptionParser& optionParser)
{
    coefDomain_.x_global_domain_min = optionParser.domain_min(0);
    coefDomain_.y_global_domain_min = optionParser.domain_min(1);
    coefDomain_.z_global_domain_min = optionParser.domain_min(2);

    coefDomain_.x_global_domain_max = optionParser.domain_min(0) + optionParser.domain_length(0);
    coefDomain_.y_global_domain_max = optionParser.domain_min(1) + optionParser.domain_length(1);
    coefDomain_.z_global_domain_max = optionParser.domain_min(2) + optionParser.domain_length(2);
}


void  Parameters::
copyCoefDomain(const CoefDomain& coefDomain)
{
    coefDomain_.x_global_domain_min = coefDomain.x_global_domain_min;
    coefDomain_.y_global_domain_min = coefDomain.y_global_domain_min;
    coefDomain_.z_global_domain_min = coefDomain.z_global_domain_min;

    coefDomain_.x_global_domain_max = coefDomain.x_global_domain_max;
    coefDomain_.y_global_domain_max = coefDomain.y_global_domain_max;
    coefDomain_.z_global_domain_max = coefDomain.z_global_domain_max;
}
