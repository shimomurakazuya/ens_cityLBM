#pragma once
#ifndef NUDGINGCOEFADAPTATION_H_
#define NUDGINGCOEFADAPTATION_H_


#include <iostream>
#include <string>
#include "definePrecision.h"
#include "defineLBM.h"
#include "Field.h"
#include "MPICommEnsemble.h"
#include "Date.h"
#include "FastPostprocessDA.h"
#include "ParticleFilter.h"


class  NudgingCoefAdaptation {
private:
    const MPICommEnsemble comm_;

    ParticleFilter particleFilter_;

    /// hyper parameters ///
    static constexpr double sigma_x_ = 0.25;
    
    #ifdef NUDSTAT_VELSDEV // ( u' u' )^0.5
    static constexpr double sigma_y_ = 0.20;
    #else
    static constexpr double sigma_y_ = 0.020;
    #endif

    static constexpr float coef_nudging_init_{1.0};
    static constexpr float coef_nudging_min_{0.10};
    static constexpr float coef_nudging_max_{10.0};

public:
    NudgingCoefAdaptation() = delete;
    NudgingCoefAdaptation(const MPICommEnsemble comm): comm_(comm), particleFilter_(comm)
    {
        particleFilter_.init( comm_.row_vector().size(), coef_nudging_init_, sigma_x_, sigma_y_, coef_nudging_min_, coef_nudging_max_ );
    }
    ~NudgingCoefAdaptation(){}

private:
    ParticleFilter& particleFilter() { return particleFilter_;}

public:
    float coef_nudging    () const { return particleFilter_.x()[ comm_.row_vector().rank() ]; }
    float coef_ave_nudging() const { return particleFilter_.x().sum() / particleFilter_.x().size(); }
    float coef_min_nudging() const { return particleFilter_.x().min(); }
    float coef_max_nudging() const { return particleFilter_.x().max(); }

public:
    void UpdateCoef(
        const int           
            #ifndef USE_WRF0
            time_num
            #endif
            ,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
        const MeshValue*    meshValues,
        const ValueTimeAverage* valueTimeAverage,
        const bool          use_pf
        )
    {
        if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__; }

        if (use_pf == false) { return; }
        if (comm_.is_rank0()) { std::cout << std::endl << __PRETTY_FUNCTION__ << " : time = " << time_num << " minutes" << std::endl; }
    
        std::vector<double> _errors_ens(comm_.row_vector().size()); for (int i=0; i<comm_.row_vector().size(); i++) { _errors_ens[i] = 0.0; }
        std::vector<double>  errors_ens(comm_.row_vector().size()); for (int i=0; i<comm_.row_vector().size(); i++) {  errors_ens[i] = 0.0; }

        const double error_ens = MeanErrorObservedData_st1(time_num, grids, tree, parameters, meshValues, valueTimeAverage);
        _errors_ens[ comm_.row_vector().rank() ] = error_ens;

        if (comm_.is_rank0()) { std::cout << "MPI_Allreduce" << std::endl; }
        MPI_Allreduce(_errors_ens.data(), errors_ens.data(), comm_.row_vector().size(), MPI_DOUBLE, MPI_MAX, comm_.row_vector().comm());

        if (comm_.is_rank0()) { std::cout << "particleFilter_.update" << std::endl; }
        particleFilter_.update( errors_ens.data(), errors_ens.size() );
    }


    void load_NudgingCoefAdaptation(const int t)
    {
        particleFilter_.load_ParticleFiter(t);
    }

    void save_NudgingCoefAdaptation(const int t) const
    {
        particleFilter_.save_ParticleFiter(t);
    }

private:
    float MeanErrorObservedData_st1(
        const int           
            #ifndef USE_WRF0
            time_num
            #endif
            ,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
        const MeshValue*    meshValues,
        const ValueTimeAverage* valueTimeAverage
        )
    {
        if (comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__; }

        Date date = timenum2utc(time_num);
        const std::string fname = Foldernames::input_folder + "/observed/st1/2003/07/16/" + date.hhmmss_str() + "_utc.dat";
    
        std::vector<float> x, y, z;
        struct val_st1 {
            float umean, vmean, wmean;
            float usdev, vsdev, wsdev;
        };
    
        std::vector<val_st1> val;
        std::string str_buf1, str_buf2;
    
        std::ifstream fin(fname.c_str(), std::ios::in);
        while (std::getline(fin, str_buf1)) {    
            std::istringstream i_stream(str_buf1);
    
            std::vector<float> vt;
            while (getline(i_stream, str_buf2, ' ')) {
                vt.push_back(stof(str_buf2));
            }
    
            x.emplace_back( vt[0] );
            y.emplace_back( vt[1] );
            z.emplace_back( vt[2] );
            val.emplace_back( val_st1{ vt[3], vt[4], vt[5], vt[6], vt[7], vt[8] } );
        }
        fin.close();
    
        FastPostprocessDA postDA(comm_);
        postDA.setupdata(x,y,z, tree, meshValues);
        postDA.InterpolateMonitorData(tree, parameters, meshValues, valueTimeAverage);
    
        float error = 0.0;
        for (int i=0; i<val.size(); i++) {
            if (postDA.calcflg()[i] != 1) { continue; }
    
            #ifdef NUDSTAT_VELSDEV // u' u'

            const float vel2_sdev_obs = sqrt( pow(val[i].usdev, 2) + pow(val[i].vsdev, 2) + pow(val[i].wsdev, 2) );
            const float vel2_sdev_cal = sqrt( postDA.vel2_fluc()[i] );

            const float eval_obs = vel2_sdev_obs;
            const float eval_cal = vel2_sdev_cal;

            #else // intensity of turbulence [-]

            constexpr float ep = 1.0e-8;
            const float vel_sdev_obs = sqrt( pow(val[i].usdev, 2) + pow(val[i].vsdev, 2) + pow(val[i].wsdev, 2) ) + ep;
            const float vel_mean_obs = sqrt( pow(val[i].umean, 2) + pow(val[i].vmean, 2) + pow(val[i].wmean, 2) ) + ep;

            const float vel_sdev_cal = sqrt( postDA.vel2_fluc()[i] ) + ep;
            const float vel_mean_cal = sqrt( pow(postDA.u_mean()[i], 2) + pow(postDA.v_mean()[i], 2) + pow(postDA.w_mean()[i], 2) ) + ep;

            const float eval_obs = vel_sdev_obs / vel_mean_obs;
            const float eval_cal = vel_sdev_cal / vel_mean_cal;

            #endif
            const float _error_abs = fabs( eval_obs - eval_cal );

            const float error_lim = 2.5*sigma_y_;
            const float error_abs = (_error_abs < error_lim) ? _error_abs : error_lim;
    
            error += error_abs;
        }
    
        // evaluate error
        const auto comm_target = comm_.col_vector().comm();
    
        float error_ens = 0.0;
        MPI_Allreduce(&error, &error_ens, 1, MPI_FLOAT, MPI_SUM, comm_target);

        return error_ens/val.size();
    }

public:

private:
    Date timenum2utc(int timenum){
        constexpr int timenum_base = 360; // minutes
        constexpr int hh_utc_base  = 12;

        const int timenum_offset = (timenum - timenum_base) + 60*hh_utc_base; // timenum from 00:00:00 utc

        const int ss_offset = timenum_offset * 60;
        const int mm_offset = timenum_offset;
        const int hh_offset = timenum_offset / 60;

        const int ss = ss_offset%60;
        const int mm = mm_offset%60;
        const int hh = hh_offset%24;

        return Date(hh, mm, ss);
    }
};


#endif

