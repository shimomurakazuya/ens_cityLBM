#pragma once
#ifndef FUNCHEATFLUXMODEL_H__
#define FUNCHEATFLUXMODEL_H__


#include <iostream>
#include "definePrecision.h"
#include "defineFluidProperty.h"


namespace  FuncHeatFluxModel {


template<typename T>
__HOST__ __DEVICE__
T q3w_HeatFlux(const T u_chi, const T T_chi, const T Tb, const T dzs, const T Tref)
{
    constexpr int iter_max = 5;

    constexpr double ep    = 5.0e-5;
//    constexpr double ep    = 1.0e-5;
//    constexpr double ep    = 1.0e-6; // nan
    //    desert
//    constexpr double z0m   = 0.00267;
//    constexpr double z0h   = 0.000049;
    //    Gobi
    constexpr double z0m   = 0.00280 + 0.00030;
    constexpr double z0h   = 0.000011;

    //    LOHDIM-LES
//    constexpr double z0m   = 0.3;
//    constexpr double z0h   = 0.0001;

    constexpr double kappa = 0.4;
    constexpr double g     = 9.80665;
//    constexpr double g     = -9.80665;

    double  L = (Tb - T_chi) >= 0.0 ?
                            (15.0*dzs + ep) * 100.0 
                         : -(15.0*dzs + ep) * 100.0;
    double  psi_m = 0.0;
    double  psi_h = 0.0;

    double  u_s = 0.0;
//    T  u_s = sqrt(kappa * u_chi/dzs);

    double  q3w = 0.0;

    auto lim_zero_div = [](const T a){
        auto sgn_a = a >= 0.0 ? 1.0 : -1.0;
        return sgn_a * (fabs(a) + ep);
    };

//    auto init_x_L = [&](T& _x, T& _L, const T _Tb, const T _T_chi){
//        if ( (_Tb - _T_chi) > 0.0 ) { _L =  150.0; }
//        else                        { _L = -150.0; }
//
//        _x = pow( (1.0 - 15.0*dzs/_L), 0.25 );
//    };

    auto calc_x = [&](const double _L){
        auto _L_tmp = lim_zero_div(_L);
        auto _tmp =  (1.0 - 15.0*dzs/_L_tmp);
        auto tmp = (_tmp < ep) ? ep : _tmp;

        return pow(tmp, 0.25);
    };

    auto update_L = [&](double& _L, const double _u_s, const double _q3w){
        auto _q3w_tmp = lim_zero_div(_q3w);
        _L = - (pow(fabs(_u_s), 3) * Tref) / (kappa * g * _q3w_tmp);
    };

    auto update_psi_m = [&](double& _psi_m, const double _L){
        auto _x = calc_x(_L);
        auto _L_tmp = lim_zero_div(_L);
        auto tmp_st = -4.7 * dzs/_L_tmp;
        auto tmp_un = 2.0*log(0.5*(1.0+_x)) + log(0.5*(1.0+_x*_x)) - 2.0*atan(_x) + M_PI*0.5;

        _psi_m = (_L >= 0.0) ? tmp_st : tmp_un;
    };

    auto update_psi_h = [&](double& _psi_h, const double _L){
        auto _x = calc_x(_L);
        auto _L_tmp = lim_zero_div(_L);
        auto tmp_st = -7.8 * dzs/_L_tmp;
        auto tmp_un = 2.0*log(0.5*(1.0+_x*_x));

        _psi_h = (_L >= 0.0) ? tmp_st : tmp_un;
    };

    auto update_u_s = [&](double& _u_s, const double _psi_m){
        auto fdiv = fabs( lim_zero_div(log(dzs/z0m) - _psi_m) );

        _u_s = (kappa * u_chi) / fdiv;
    };

    auto update_q3w = [&](double& _q3w, const double _psi_h){
        auto fdiv = fabs( lim_zero_div(log(dzs/z0h) - _psi_h) );

        _q3w = (kappa * u_chi * (Tb - T_chi)) / fdiv;
    };


    for (int _i=0; _i<iter_max; _i++) {
        update_psi_m(psi_m, L);
        update_psi_h(psi_h, L);
        update_u_s  (u_s, psi_m);
        update_q3w  (q3w, psi_h);
        update_L    (L, u_s, q3w);
    }

    return (T)q3w;
}


};


#endif
