#pragma once
#ifndef FUNCDIFFERENTIALOPERATOR_H_
#define FUNCDIFFERENTIALOPERATOR_H_


#include "FuncGradientOperator.h"
#include "FuncObj.h"


namespace FuncDifferentialOperator {


// euler
template<typename T>
__HOST__ __DEVICE__
T flux_euler_1st(T f0, T f1, T vel, double dx, double dt)
{
    const T fs = FuncGradientOperator::fs_upwind_c2s(f0, f1, vel);

    return fs*vel*dt;
}


template<typename T>
__HOST__ __DEVICE__
T flux_euler_2nd(T f0, T f1, T vel, double dx, double dt)
{
//    const T fs = fs_upwind_c2s(f0, f1, vel);
    const T fs = FuncGradientOperator::fs_central_v2s(f0, f1)    ;

    return fs*vel*dt;
}


template<typename T>
__HOST__ __DEVICE__
T flux_euler_4th(T f0, T f1, T f2, T f3, T vel, double dx, double dt)
{
    const T fs  = FuncGradientOperator::fs_central_v2s(f0, f1, f2, f3)    ;

    return fs*vel*dt;
}


template<typename T>
__HOST__ __DEVICE__
T flux_euler_6th(T f0, T f1, T f2, T f3, T f4, T f5, T vel, double dx, double dt)
{
    using namespace FuncGradientOperator;
    const T fs  = fs_central_v2s(f0, f1, f2, f3, f4, f5);

    return fs*vel*dt;
}


template<typename T>
__HOST__ __DEVICE__
T flux_euler_weno(T f0, T f1, T f2, T f3, T f4, T f5, T vel, double dx, double dt)
{
    using namespace FuncGradientOperator;
    const T fs  = weno_v2s(f0, f1, f2, f3, f4, f5, vel);

    return fs*vel*dt;
}


// semilag
template<typename T>
__HOST__ __DEVICE__
T flux_semilag_2nd(T f0, T f1, T vel, double dx, double dt)
{
    using namespace FuncGradientOperator;
//    const T fs = fs_upwind_c2s(f0, f1, vel);
    const T fs = fs_central_v2s(f0, f1)    ;
    const T fx = fx_central_v2s(f0, f1, dx);

    const T a0 = fs;
    const T ax = fx;

    const T mveldt = -vel*dt;
    return - ( a0*mveldt + ax*mveldt*mveldt*(T)0.5 );
}


template<typename T>
__HOST__ __DEVICE__
T flux_semilag_4th(T f0, T f1, T f2, T f3, T vel, double dx, double dt)
{
    using namespace FuncGradientOperator;

    const T fs  = fs_central_v2s (f0, f1, f2, f3)    ;
    const T fx  = fx_central_v2s (f0, f1, f2, f3, dx);
    const T fxx = fxx_central_v2s(f0, f1, f2, f3, dx);
    const T fx3 = fx3_central_v2s(f0, f1, f2, f3, dx);

    const T a0  = fs  ;
    const T ax  = fx  ;
    const T axx = fxx / ((T)2.0);
    const T ax3 = fx3 / ((T)6.0);

    const T mveldt = -vel*dt;
    return - ( a0*mveldt + ax*pow(mveldt, 2)*(T)0.5 + (T)axx*pow(mveldt, 3)/(T)3.0 + (T)ax3*pow(mveldt, 4)/(T)4.0 );
}


template<typename T>
__HOST__ __DEVICE__
T flux_semilag_6th(T f0, T f1, T f2, T f3, T f4, T f5, T vel, double dx, double dt)
{
    using namespace FuncGradientOperator;

    const T fs  = fs_central_v2s (f0, f1, f2, f3, f4, f5);
    const T fx  = fx_central_v2s (f0, f1, f2, f3, f4, f5, dx);
    const T fxx = fxx_central_v2s(f0, f1, f2, f3, f4, f5, dx);
    const T fx3 = fx3_central_v2s(f0, f1, f2, f3, f4, f5, dx);
    const T fx4 = fx4_central_v2s(f0, f1, f2, f3, f4, f5, dx);
    const T fx5 = fx5_central_v2s(f0, f1, f2, f3, f4, f5, dx);

    const T a0  = fs ;
    const T ax  = fx ;
    const T axx = fxx / ((T)2.0);
    const T ax3 = fx3 / ((T)6.0);
    const T ax4 = fx4 / ((T)24.0);
    const T ax5 = fx5 / ((T)120.0);

    const T mveldt = -vel*dt;
    return - ( a0*mveldt + ax*pow(mveldt, 2)*(T)0.5 + (T)axx*pow(mveldt, 3)/(T)3.0 + (T)ax3*pow(mveldt, 4)/(T)4.0 + (T)ax4*pow(mveldt, 5)/(T)5.0 + (T)ax5*pow(mveldt, 6)/(T)6.0 );
}


template<typename T>
__HOST__ __DEVICE__
T flux_semilag_weno(T f0, T f1, T f2, T f3, T f4, T f5, T vel, double dx, double dt)
{
    using namespace FuncGradientOperator;
    const T fs  = weno_v2s       (f0, f1, f2, f3, f4, f5, vel);
    const T fx  = fx_central_v2s (f0, f1, f2, f3, f4, f5, dx);
    const T fxx = fxx_central_v2s(f0, f1, f2, f3, f4, f5, dx);
    const T fx3 = fx3_central_v2s(f0, f1, f2, f3, f4, f5, dx);
    const T fx4 = fx4_central_v2s(f0, f1, f2, f3, f4, f5, dx);
    const T fx5 = fx5_central_v2s(f0, f1, f2, f3, f4, f5, dx);

    const T a0  = fs ;
    const T ax  = fx ;
    const T axx = fxx / ((T)2.0);
    const T ax3 = fx3 / ((T)6.0);
    const T ax4 = fx4 / ((T)24.0);
    const T ax5 = fx5 / ((T)120.0);

    const T mveldt = -vel*dt;
    return - ( a0*mveldt + ax*pow(mveldt, 2)*(T)0.5 + (T)axx*pow(mveldt, 3)/(T)3.0 + (T)ax3*pow(mveldt, 4)/(T)4.0 + (T)ax4*pow(mveldt, 5)/(T)5.0 + (T)ax5*pow(mveldt, 6)/(T)6.0 );
}


template<typename T>
__HOST__ __DEVICE__
T flux_euler(T f0, T f1, T f2, T f3, T f4, T f5, T vel, double dx, double dt)
{
    return flux_euler_weno(f0, f1, f2, f3, f4, f5, vel, dx, dt);

//    constexpr T weight = 0.99;
//    return          weight  * flux_euler_weno(f0, f1, f2, f3, f4, f5, vel, dx, dt)
//           + (1.0 - weight) * flux_euler_1st(f2, f3, vel, dx, dt);
}


template<typename T>
__HOST__ __DEVICE__
T flux_semilag(T f0, T f1, T f2, T f3, T f4, T f5, T vel, double dx, double dt)
{
//    return flux_semilag_2nd(f2, f3, vel, dx, dt);
//    return flux_semilag_4th(f1, f2, f3, f4, vel, dx, dt);
//    return flux_semilag_6th(f0, f1, f2, f3, f4, f5, vel, dx, dt);
    return flux_semilag_weno(f0, f1, f2, f3, f4, f5, vel, dx, dt);
}


template<typename T>
__HOST__ __DEVICE__
T convection1d_fvm_euler_st7(
    T f0, T f1, T f2, T f3, T f4, T f5, T f6,
    T vel0, T vel1,
    double dx,   
    double dt
    )
{
    const T flux0 =  flux_euler(f0, f1, f2, f3, f4, f5, vel0, dx, dt);
    const T flux1 =  flux_euler(f1, f2, f3, f4, f5, f6, vel1, dx, dt);

    return - (flux1 - flux0)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T convection1d_fvm_euler_w_bc_st7(
    T f0, T f1, T f2, T f3, T f4, T f5, T f6,
    T vel0, T vel1,
    T lv0,  T lv1,
    double dx,   
    double dt
    )
{
    const T flux0 = FuncObj::is_fluid(lv0) ?  flux_euler(f0, f1, f2, f3, f4, f5, vel0, dx, dt) : (T)0.0;
    const T flux1 = FuncObj::is_fluid(lv1) ?  flux_euler(f1, f2, f3, f4, f5, f6, vel1, dx, dt) : (T)0.0;

    return - (flux1 - flux0)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T convection1d_fvm_semilag_w_bc(
    T f0, T f1, T f2, T f3, T f4, T f5, T f6,
    T vel0, T vel1,
    T lv0,  T lv1,
    double dx,   
    double dt
    )
{
    const T flux0 = FuncObj::is_fluid(lv0) ?  flux_semilag(f0, f1, f2, f3, f4, f5, vel0, dx, dt) : (T)0.0;
    const T flux1 = FuncObj::is_fluid(lv1) ?  flux_semilag(f1, f2, f3, f4, f5, f6, vel1, dx, dt) : (T)0.0;

    return - (flux1 - flux0)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T advection1d_fdm_euler_st3(
    T f0, T f1, T f2,
    T u0, T u1, T u2,
    double dx,
    double dt
    )
{
    const T d0 = f1 - f0;
    const T d1 = f2 - f1;

    const T vel = u1;
    const T fxu = FuncGradientOperator::fs_upwind_c2s (d0, d1, u1);
    const T fxc = FuncGradientOperator::fs_central_c2s(d0, d1);

    const T w  = sqrt( fabs( d0 - d1 ) / ( fabs(d0) + fabs(d1) + 1.0e-7 ) );
    const T fx = w*fxu + ((T)1.0 - w)*fxc;

    return - (vel*fx*dt)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T advection1d_fdm_euler_neumann_st3(
    T f0, T f1, T f2,
    T u0, T u1, T u2,
    bool obj0, bool obj1, bool obj2,
    double dx,
    double dt
    )
{
    const T d0 = (obj0 || obj1) ? 0.0 : f1 - f0;
    const T d1 = (obj1 || obj2) ? 0.0 : f2 - f1;

    const T us0 = obj0 ? 0.0 : (T)0.5 * (u0 + u1);
    const T us1 = obj2 ? 0.0 : (T)0.5 * (u2 + u1);

    const T vel = u1;
    const T fxu = FuncGradientOperator::fs_upwind_c2s (d0, d1, u1);
    const T fxc = FuncGradientOperator::fs_central_c2s(d0, d1);

    const T w  = sqrt( fabs( d0 - d1 ) / ( fabs(d0) + fabs(d1) + 1.0e-7 ) );
    const T fx = w*fxu + ((T)1.0 - w)*fxc;

    return - (vel*fx*dt)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T advection1d_fdm_euler_st7(
    T f0, T f1, T f2, T f3, T f4, T f5, T f6,
    T vel,
    double dx,
    double dt
    )
{
    const T d0 = f1 - f0;
    const T d1 = f2 - f1;
    const T d2 = f3 - f2;
    const T d3 = f4 - f3;
    const T d4 = f5 - f4;
    const T d5 = f6 - f5;

    const T fx = FuncGradientOperator::weno_v2s(d0, d1, d2, d3, d4, d5, vel);

    return - (vel*fx*dt)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T advection1d_fdm_euler_neumann_st7(
    T f0, T f1, T f2, T f3, T f4, T f5, T f6,
    T vel,
    bool obj0, bool obj1, bool obj2, bool obj3, bool obj4, bool obj5, bool obj6,
    double dx,
    double dt
    )
{
    const T d0 = (obj0 || obj1 || obj2) ? 0.0 : f1 - f0;
    const T d1 = (obj1 || obj2        ) ? 0.0 : f2 - f1;
    const T d2 = (obj2 || obj3        ) ? 0.0 : f3 - f2;
    const T d3 = (obj3 || obj4        ) ? 0.0 : f4 - f3;
    const T d4 = (obj4 || obj5        ) ? 0.0 : f5 - f4;
    const T d5 = (obj4 || obj5 || obj6) ? 0.0 : f6 - f5;

    const T fx = FuncGradientOperator::weno_v2s(d0, d1, d2, d3, d4, d5, vel);

    return - (vel*fx*dt)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T advection1d_fdm_euler_upwind(
    T f0, T f1, T f2,
    T vel,
    double dx,
    double dt
    )
{
    const T d0 = f1 - f0;
    const T d1 = f2 - f1;

    const T fx = FuncGradientOperator::fs_upwind_c2s(d0, d1, vel);

    return - (vel*fx*dt)/dx;
}


} // namespace //


#endif
