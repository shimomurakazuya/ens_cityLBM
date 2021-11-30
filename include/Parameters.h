#pragma once
#ifndef PARAMETERS_H_
#define PARAMETERS_H_


#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include "defineCUDA.h"


#include "definePrecision.h"
#include "Parser.h"

#include "MPICommEnsemble.h"


struct  GPUDevice {
    int  id_gpu;
    int  gpu_per_node;
};


struct  CoefFluidPropertyLBM {
    real  c_ref_lbm;
};


struct  CoefTime {
    // second //
    real  time_now;
    real  time_init;
    real  time_end;
    real  dt0;

    int   step_now;
    int   step_init;
    int   step_end;
};


struct  CoefCFROutput {
    int  cout_step = 1;
    int  fout_step = 1;

    bool  cout_flag;
    bool  fout_flag;


    // func //
    int   t_cout(int t) const { return  t/cout_step; }
    int   t_fout(int t) const { return  t/fout_step; }

    bool  is_cout_step       (int t) const { return  (t%cout_step        == 0); }
    bool  is_fout_step       (int t) const { return  (t%fout_step        == 0); }

    bool  is_cout_step(int t, int step_end) const { return ( (t%cout_step == 0) || (t == step_end) ); }
    bool  is_fout_step(int t, int step_end) const { return ( (t%fout_step == 0) || (t == step_end) ); }

    bool  is_just_before_cout_step       (int t) const { return  (t%cout_step        == cout_step       -1); }
    bool  is_just_before_fout_step       (int t) const { return  (t%fout_step        == fout_step       -1); }
};


struct  CoefDomain {
    real  x_global_domain_min;
    real  y_global_domain_min;
    real  z_global_domain_min;

    real  x_global_domain_max;
    real  y_global_domain_max;
    real  z_global_domain_max;

    real x_global_domain_center() const { return (x_global_domain_min + x_global_domain_max)*0.5; }
    real y_global_domain_center() const { return (y_global_domain_min + y_global_domain_max)*0.5; }
    real z_global_domain_center() const { return (z_global_domain_min + z_global_domain_max)*0.5; }
};


struct  stIOStep {
    int  step;
    real time;
    int  fstep;
};


class  Parameters {
private:
    const MPICommEnsemble comm_;
    GPUDevice             gpuDevice_;

    CoefFluidPropertyLBM  coefFluidPropertyLBM_;
    CoefTime              coefTime_;
    CoefCFROutput         coefCFROutput_;
    CoefDomain            coefDomain_;

    int                   gpu_per_node_;

public:
    Parameters () = delete;
    Parameters (const MPICommEnsemble comm):  comm_(comm) {}
    ~Parameters () {}

public:
    // CoefFluidPropertyLBM //
    real  c_ref_lbm() const { return coefFluidPropertyLBM_.c_ref_lbm; }

    // CoefTime //
    real  time_init() const { return coefTime_.time_init; }
    real  time_end () const { return coefTime_.time_end; }
    real  time_now()  const { return coefTime_.time_now; }
    real  dt0()       const { return coefTime_.dt0; }

    int   step_now()  const { return coefTime_.step_now; }
    int   step_init() const { return coefTime_.step_init; }
    int   step_end()  const { return coefTime_.step_end; }

    void  time_step_update(int step, real time) {
        coefTime_.step_now = step;
        coefTime_.time_now = time;
    }

    // struct //
    const GPUDevice&            gpuDevice()             const { return  gpuDevice_; }
    const CoefFluidPropertyLBM& coefFluidPropertyLBM()  const { return  coefFluidPropertyLBM_; }
    const CoefTime&             coefTime()              const { return  coefTime_; }
    const CoefCFROutput&        coefCFROutput()         const { return  coefCFROutput_; }
    const CoefDomain&           coefDomain()            const { return  coefDomain_; }

public:
    bool  is_step_end(int t) const { return (t == step_end() - 1); }

public:
    void  init(const OptionParser& optionParser);
    void  copy(const Parameters& parameters);

    void  readParameters (const OptionParser&  optionParser, const std::string  filename);
    void  writeParameters(const std::string  filename, int step) const;
    void  coutParameters() const;

    std::string  filename(const int step) const;

private:
    void  initGPUDevice(const OptionParser& optionParser);
    void  copyGPUDevice(const GPUDevice& gpuDevice);

    void  initCoefFluidPropertyLBM(const OptionParser& optionParser);
    void  copyCoefFluidPropertyLBM(const CoefFluidPropertyLBM& coefFluidPropertyLBM);

    void  initCoefTime(const OptionParser& optionParser);
    void  copyCoefTime(const CoefTime& coefTime);

    void  initCoefCFROutput(const OptionParser& optionParser);
    void  copyCoefCFROutput(const CoefCFROutput& coefCFROutput);

    void  initCoefDomain(const OptionParser& optionParser);
    void  copyCoefDomain(const CoefDomain& coefDomain);
};


#endif
