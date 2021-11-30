#pragma once
#ifndef PARTICLEFILTERST_H_
#define PARTICLEFILTERST_H_


#include <iostream>
#include <string>
#include "definePrecision.h"
#include "defineLBM.h"
#include "Field.h"
#include "MPICommEnsemble.h"
#include "Date.h"
#include "FastPostprocessDA.h"


class  ParticleFilterSt {
private:
    const MPICommEnsemble comm_;

public:
    ParticleFilterSt() = delete;
    ParticleFilterSt(const MPICommEnsemble comm): comm_(comm) {}
    ~ParticleFilterSt(){}

private:
//    float coef_mean_nudging_{0.5};
    float coef_mean_nudging_{1.0};
//    float coef_mean_nudging_{2.0};
//    float coef_mean_nudging_{4.0};

    float coef_wid_nudging_{0.5}; // < 1.0
//    float coef_wid_nudging_{0.8}; // < 1.0

public:
    float coef_mean_nudging() const { return coef_mean_nudging_; }
    float coef_wid_nudging()  const { return coef_wid_nudging_; }
    float coef_nudging()      const { return f_with_ens_spread(coef_mean_nudging_, coef_wid_nudging_); }

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
              MeshValue*    meshValues,
        const ValueTimeAverage* valueTimeAverage,
        const bool          use_pf = true
        )
    {
        if (use_pf == false) { return; }
        if (comm_.is_rank0()) { std::cout << std::endl << __PRETTY_FUNCTION__ << " : time = " << time_num << " minutes" << std::endl; }
        const auto comm_target = comm_.col_vector().comm();
    
        const float error_ens = ErrorObservedData_st1(time_num, grids, tree, parameters, meshValues, valueTimeAverage);
    
        constexpr float ep = 1.0e-9;
        float _eval_ens = 1.0/(error_ens + ep);
        float _eval_g = 0.0;
        MPI_Allreduce(&_eval_ens, &_eval_g, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        _eval_g /=  comm_.col_vector().size();
    
        const float eval_ens = _eval_ens / _eval_g;

        const float coef_nudging_ens = f_with_ens_spread(coef_mean_nudging_, coef_wid_nudging_);

        float coef_nud_tmp = eval_ens * coef_nudging_ens;
        float coef_nud_tmp_g = 0.0;
        MPI_Allreduce(&coef_nud_tmp, &coef_nud_tmp_g, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

        coef_nud_tmp_g /= comm_.col_vector().size();

        coef_mean_nudging_ = coef_nud_tmp_g;
    }

private:
    float ErrorObservedData_st1(
        const int           
            #ifndef USE_WRF0
            time_num
            #endif
            ,
        const Grid*         grids,
        const Tree&         tree,
        const Parameters&   parameters,
              MeshValue*    meshValues,
        const ValueTimeAverage* valueTimeAverage
        )
    {
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
    
            #if 0 // velocity fluctuation [m^2/s^2]

            const float vel2_fluc_obs = pow(val[i].usdev, 2) + pow(val[i].vsdev, 2) + pow(val[i].wsdev, 2);
            const float vel2_fluc_cal = postDA.vel2_fluc()[i];

            const float eval_obs = vel2_fluc_obs;
            const float eval_cal = vel2_fluc_cal;

            #else // intensity of turbulence [-]

            constexpr float ep = 1.0e-8;
            const float vel_sdev_obs = sqrt( pow(val[i].usdev, 2) + pow(val[i].vsdev, 2) + pow(val[i].wsdev, 2) ) + ep;
            const float vel_mean_obs = sqrt( pow(val[i].umean, 2) + pow(val[i].vmean, 2) + pow(val[i].wmean, 2) ) + ep;

            const float vel_sdev_cal = sqrt( postDA.vel2_fluc()[i] ) + ep;
            const float vel_mean_cal = sqrt( pow(postDA.u_mean()[i], 2) + pow(postDA.v_mean()[i], 2) + pow(postDA.w_mean()[i], 2) ) + ep;

            const float eval_obs = vel_sdev_obs / vel_mean_obs;
            const float eval_cal = vel_sdev_cal / vel_mean_cal;

            #endif
    
            error += fabs( eval_obs - eval_cal );
        }
    
        // evaluate error
        const auto comm_target = comm_.col_vector().comm();
    
        float error_ens = 0.0;
        MPI_Allreduce(&error, &error_ens, 1, MPI_FLOAT, MPI_SUM, comm_target);

        return error_ens;
    }


public:
    float f_with_ens_spread(const float mean, const float wid) const
    {
        return comm_.n_cols() == 1
           ? mean
           : mean * ( (1.0-wid) + comm_.ensemble_id() * 2.0*wid/(comm_.ensemble_size() - 1) );
    }

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

