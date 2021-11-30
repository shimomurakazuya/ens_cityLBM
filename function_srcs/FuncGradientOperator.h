#pragma once
#ifndef FUNCGRADIENTOPERATOR_H_
#define FUNCGRADIENTOPERATOR_H_


namespace FuncGradientOperator {


// c : center
// s : surface
// v : volume


// 2 stencils //
template<typename T>
__HOST__ __DEVICE__
T fs_upwind_c2s(T f0, T f1, T vel)
{
    return  vel > (T)0.0 ? f0 : f1;
}

template<typename T>
__HOST__ __DEVICE__
T fs_central_c2s(T f0, T f1)
{
    return  (f0 + f1) * (T)0.5;
}

template<typename T>
__HOST__ __DEVICE__
T fx_central_c2s(T f0, T f1, double dx)
{
    return  (f1 - f0)/dx;
}

// 2 stencils (fvm) //
template<typename T>
__HOST__ __DEVICE__
T fs_upwind_v2s(T f0, T f1, T vel)
{
    return  vel > (T)0.0 ? f0 : f1;
}

template<typename T>
__HOST__ __DEVICE__
T fs_central_v2s(T f0, T f1)
{
    return  (f0 + f1) * (T)0.5;
}

template<typename T>
__HOST__ __DEVICE__
T fx_central_v2s(T f0, T f1, double dx)
{
    return  (f1 - f0)/dx;
}

// 3 stelcils //
template<typename T>
__HOST__ __DEVICE__
T fx_upwind_c2c(T f0, T f1, T f2, double dx, T vel)
{
    return  vel > (T)0.0 ? (f1 - f0)/dx : (f2 - f1)/dx;
}


template<typename T>
__HOST__ __DEVICE__
T fx_central_c2c(T f0, T f1, T f2, double dx)
{
    return  (f2 - f0)/((T)2.0*dx);
}


template<typename T>
__HOST__ __DEVICE__
T fxx_central_c2c(T f0, T f1, T f2, double dx)
{
    return  (f2 - (T)2.0*f1 + f0)/(dx*dx);
}

// 4 stencils //
template<typename T>
__HOST__ __DEVICE__
T fs_central_c2s(T f0, T f1, T f2, T f3)
{
    return  ( -f0 + (T)9.0*f1 + (T)9.0*f2 - f3 ) / ((T)16.0);
}

template<typename T>
__HOST__ __DEVICE__
T fx_central_c2s(T f0, T f1, T f2, T f3, double dx)
{
    return  (f0 - 27.0*f1 + 27.0*f2 - f3) / ( (T)24.0 * dx );
}

template<typename T>
__HOST__ __DEVICE__
T fxx_central_c2s(T f0, T f1, T f2, T f3, double dx)
{
    return  (f0 - f1 - f2 + f3) / ( (T)2.0 * dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx3_central_c2s(T f0, T f1, T f2, T f3, double dx)
{
    return  ( -f0 + (T)3.0*f1 - (T)3.0*f2 + f3) / ( dx*dx*dx );
}

// 4 stencils (fvm) //
template<typename T>
__HOST__ __DEVICE__
T fs_central_v2s(T f0, T f1, T f2, T f3)
{
    return  ( -f0 + (T)7.0*f1 + (T)7.0*f2 - f3 ) / ((T)12.0);
}

template<typename T>
__HOST__ __DEVICE__
T fs_upwind_v2s(T f0, T f1, T f2, T f3, T vel)
{
    return  vel > (T)0.0 ?
          (  -f0 + (T)5.0*f1 - (T)2.0*f2 ) / ((T)6.0)
        : (  -f3 + (T)5.0*f2 - (T)2.0*f1 ) / ((T)6.0) ;
}

template<typename T>
__HOST__ __DEVICE__
T fx_central_v2s(T f0, T f1, T f2, T f3, double dx)
{
    return  (f0 - 15.0*f1 + 15.0*f2 - f3) / ( (T)12.0 * dx );
}

template<typename T>
__HOST__ __DEVICE__
T fxx_central_v2s(T f0, T f1, T f2, T f3, double dx)
{
    return  (f0 - f1 - f2 + f3) / ( (T)2.0 * dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx3_central_v2s(T f0, T f1, T f2, T f3, double dx)
{
    return  ( -f0 + (T)3.0*f1 - (T)3.0*f2 + f3) / ( dx*dx*dx );
}

// 5 stencils
template<typename T>
__HOST__ __DEVICE__
T fx_central_c2c(T f0, T f1, T f2, T f3, T f4, double dx)
{
    return (  f0 - (T)8.0*f1 + (T)8.0*f3 - f4 ) / ( (T)12.0*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fxx_central_c2c(T f0, T f1, T f2, T f3, T f4, double dx)
{
    return ( -f0 + (T)16.0*f1 - (T)30.0*f2 + (T)16.0*f3 - f4 ) / ( (T)12.0*dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx3_central_c2c(T f0, T f1, T f2, T f3, T f4, double dx)
{
    return ( -f0 + (T)2.0*f1 - (T)2.0*f3 + f4 ) / ( (T)2.0*dx*dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx4_central_c2c(T f0, T f1, T f2, T f3, T f4, double dx)
{
    return (  f0 - (T)4.0*f1 + (T)6.0*f2 - (T)4.0*f3 + f4 ) / ( dx*dx*dx*dx );
}

// 6 stencils
template<typename T>
__HOST__ __DEVICE__
T fs_central_c2s(T f0, T f1, T f2, T f3, T f4, T f5)
{
    return  (  (T)3.0*f0 - (T)25.0*f1 + (T)150.0*f2 + (T)150.0*f3 - (T)25.0*f4 + (T)3.0*f5 ) / ((T)256.0);
}

template<typename T>
__HOST__ __DEVICE__
T fx_central_c2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( (T)-9.0*f0 + (T)125.0*f1 - (T)2250.0*f2 + (T)2250.0*f3 - (T)125.0*f4 + (T)9.0*f5) / ( (T)1920.0 * dx );
}

template<typename T>
__HOST__ __DEVICE__
T fxx_central_c2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ((T)-5.0*f0 + (T)39.0*f1 - (T)34.0*f2 - (T)34.0*f3 + (T)39.0*f4 - (T)5.0*f5) / ( (T)48.0 * dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx3_central_c2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( f0 - (T)13.0*f1 + (T)34.0*f2 - (T)34.0*f3 + (T)13.0*f4 - f5 ) / ( (T)8.0*dx*dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx4_central_c2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( f0 - (T)3.0*f1 + (T)2.0*f2 + (T)2.0*f3 - (T)3.0*f4 + f5 ) / ( 2.0*dx*dx*dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx5_central_c2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( -f0 + (T)5.0*f1 - (T)10.0*f2 + (T)10.0*f3 - (T)5.0*f4 + f5 ) / ( dx*dx*dx*dx*dx );
}

// 6 stencils (fvm)
template<typename T>
__HOST__ __DEVICE__
T fs_upwind_v2s(T f0, T f1, T f2, T f3, T f4, T f5, T vel)
{
    return  vel > (T)0.0 ?
          (  (T)2.0*f0 - (T)13.0*f1 + (T)47.0*f2 + (T)27.0*f3 - (T)3.0*f4 ) / ((T)60.0)
        : (  (T)2.0*f5 - (T)13.0*f4 + (T)47.0*f3 + (T)27.0*f2 - (T)3.0*f1 ) / ((T)60.0) ;
}

template<typename T>
__HOST__ __DEVICE__
T fs_central_v2s(T f0, T f1, T f2, T f3, T f4, T f5)
{
    return  (  (T)f0 - (T)8.0*f1 + (T)37.0*f2 + (T)37.0*f3 - (T)8.0*f4 + f5 ) / ((T)60.0);
}

template<typename T>
__HOST__ __DEVICE__
T fx_central_v2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( (T)-2.0*f0 + (T)25.0*f1 - (T)245.0*f2 + (T)245.0*f3 - (T)25.0*f4 + (T)2.0*f5) / ( (T)180.0 * dx );
}

template<typename T>
__HOST__ __DEVICE__
T fxx_central_v2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( -f0 + (T)7.0*f1 - (T)6.0*f2 - (T)6.0*f3 + (T)7.0*f4 - f5) / ( (T)8.0 * dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx3_central_v2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  (  f0 - (T)11.0*f1 + (T)28.0*f2 - (T)28.0*f3 + (T)11.0*f4 - f5) / ( (T)6.0 * dx*dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx4_central_v2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  (  f0 - (T)3.0*f1 + (T)2.0*f2 + (T)2.0*f3 - (T)3.0*f4 + f5) / ( (T)2.0 * dx*dx*dx*dx );
}

template<typename T>
__HOST__ __DEVICE__
T fx5_central_v2s(T f0, T f1, T f2, T f3, T f4, T f5, double dx)
{
    return  ( -f0 + (T)5.0*f1 - (T)10.0*f2 + (T)10.0*f3 - (T)5.0*f4 + f5) / ( dx*dx*dx*dx*dx );
}

// 5 stencils
template<typename T>
__HOST__ __DEVICE__
T weno_c2c(T f0, T f1, T f2, T f3, T f4, T vel, double dx)
{
    const double ddx = 1.0/dx;

    const T fxpm = (vel >= 0.0) ?
                     (-f3 + 3.0*f2 - 3.0*f1 + f0) * (0.5*ddx)
                  : -(-f1 + 3.0*f2 - 3.0*f3 + f4) * (0.5*ddx);

    constexpr T ep = 1.0e-10;
    const T r = (vel >= 0.0) ?
                    ( pow(f2 - 2.0*f1 + f0, 2) + ep ) / ( pow(f3 - 2.0*f2 + f1, 2) + ep )
                  : ( pow(f4 - 2.0*f3 + f2, 2) + ep ) / ( pow(f3 - 2.0*f2 + f1, 2) + ep );

    const T fx_weno =  (f3 - f1)*(0.5*ddx) + fxpm/(1.0 + 2.0*r*r);

    return fx_weno;
}

// 6 stencils (fvm, weno)
template<typename T>
__HOST__ __DEVICE__
T weno_v2s(T f0, T f1, T f2, T f3, T f4, T f5, T vel)
{
    constexpr double C1312 = 13.0/12.0;
    constexpr double C16   = 1.0/6.0;
    constexpr double ep    = 1.0e-17;

    double v0, v1, v2, v3, v4;
    if (vel > (T)0.0) {
        v0 = f0;
        v1 = f1;
        v2 = f2;
        v3 = f3;
        v4 = f4;
    }
    else {
        v0 = f5;
        v1 = f4;
        v2 = f3;
        v3 = f2;
        v4 = f1;
    }

    double tmp0, tmp1;
    // IS //
    tmp0 = v0 - 2.0*v1 + v2;
    tmp1 = v0 - 4.0*v1 + 3.0*v2;
    const double is0 = C1312*tmp0*tmp0 + 0.25*tmp1*tmp1 + ep;

    tmp0 = v1 - 2.0*v2 + v3;
    tmp1 = v1 - v3;
    const double is1 = C1312*tmp0*tmp0 + 0.25*tmp1*tmp1 + ep;

    tmp0 = v2 - 2.0*v3 + v4;
    tmp1 = 3.0*v2 - 4.0*v3 + v4;
    const double is2 = C1312*tmp0*tmp0 + 0.25*tmp1*tmp1 + ep;

    const double al0 = 0.1/(is0*is0);
    const double al1 = 0.6/(is1*is1);
    const double al2 = 0.3/(is2*is2);
    const double sal = 1.0/(al0 + al1 + al2);

    const double om0 = al0*sal;
    const double om1 = al1*sal;
    const double om2 = al2*sal;

    const double ps0 = ( 2.0*v0 - 7.0*v1 + 11.0*v2)*C16;
    const double ps1 = (-1.0*v1 + 5.0*v2 +  2.0*v3)*C16;
    const double ps2 = ( 2.0*v2 + 5.0*v3 -  1.0*v4)*C16;

    return  om0*ps0 + om1*ps1 + om2*ps2;
}


template<typename T>
__HOST__ __DEVICE__
T weno_v2s_wiki(T f0, T f1, T f2, T f3, T f4, T f5, T vel)
{
    constexpr double ep    = 1.0e-18;

    double v0, v1, v2, v3, v4;
    if (vel > (T)0.0) {
        v0 = f0;
        v1 = f1;
        v2 = f2;
        v3 = f3;
        v4 = f4;
    }
    else {
        v0 = f5;
        v1 = f4;
        v2 = f3;
        v3 = f2;
        v4 = f1;
    }

    const double u0 =  3.0/8.0*v0 - 5.0/4.0*v1 + 15.0/8.0*v2;
    const double u1 = -1.0/8.0*v1 + 3.0/4.0*v2 +  3.0/8.0*v3;
    const double u2 =  3.0/8.0*v2 + 3.0/4.0*v3 -  1.0/8.0*v4;

    const double beta0 = (  4.0*v0*v0 - 19.0*v0*v1 + 25.0*v1*v1 + 11.0*v0*v2 - 31.0*v1*v2 + 10.0*v2*v2 )/3.0 + ep;
    const double beta1 = (  4.0*v1*v1 - 13.0*v1*v2 + 13.0*v2*v2 +  5.0*v1*v3 - 13.0*v2*v3 +  4.0*v3*v3 )/3.0 + ep;
    const double beta2 = ( 10.0*v2*v2 - 31.0*v2*v3 + 25.0*v3*v3 + 11.0*v2*v4 - 19.0*v3*v4 +  4.0*v4*v4 )/3.0 + ep;
    
    const double _w0 =  1.0/16.0 / (beta0*beta0);
    const double _w1 = 10.0/16.0 / (beta1*beta1);
    const double _w2 =  5.0/16.0 / (beta2*beta2);

    const double _wall  = _w0 + _w1 + _w2;
    const double _dwall = 1.0/_wall;

    const double w0 = _w0 * _dwall;
    const double w1 = _w1 * _dwall;
    const double w2 = _w2 * _dwall;

    return  w0*u0 + w1*u1 + w2*u2;
}


template<typename T>
__HOST__ __DEVICE__
T weno_v2s_yms(T f0, T f1, T f2, T f3, T f4, T f5, T vel)
{
    constexpr double C1312 = 13.0/12.0;
    constexpr double C16   = 1.0/6.0;
    constexpr double ep    = 1.0e-17;

    double v0, v1, v2, v3, v4;
    if (vel > (T)0.0) {
        v0 = f0;
        v1 = f1;
        v2 = f2;
        v3 = f3;
        v4 = f4;
    }
    else {
        v0 = f5;
        v1 = f4;
        v2 = f3;
        v3 = f2;
        v4 = f1;
    }

    double tmp0, tmp1;
    // IS //
    tmp0 = v0 - 2.0*v1 + v2;
    tmp1 = v0 - 4.0*v1 + 3.0*v2;
    const double is0 = C1312*tmp0*tmp0 + 0.25*tmp1*tmp1 + ep;

    tmp0 = v1 - 2.0*v2 + v3;
    tmp1 = v1 - v3;
    const double is1 = C1312*tmp0*tmp0 + 0.25*tmp1*tmp1 + ep;

    tmp0 = v2 - 2.0*v3 + v4;
    tmp1 = 3.0*v2 - 4.0*v3 + v4;
    const double is2 = C1312*tmp0*tmp0 + 0.25*tmp1*tmp1 + ep;

    const double al0 = 0.1/(is0*is0);
    const double al1 = 0.6/(is1*is1);
    const double al2 = 0.3/(is2*is2);
    const double sal = 1.0/(al0 + al1 + al2);

    const double om0 = al0*sal;
    const double om1 = al1*sal;
    const double om2 = al2*sal;

    const double ps0 = (2.0*v0 - 7.0*v1 + 11.0*v2)*C16;
    const double ps1 = (-v1 + 5.0*v2 + 2.0*v3)*C16;
    const double ps2 = (2.0*v2 + 5.0*v3 - v4)*C16;

    return  om0*ps0 + om1*ps1 + om2*ps2;
}


//template<typename T>
//__HOST__ __DEVICE__
//T fs_3rd_upwind(T f0, T f1, T f2, T f3, T vel)
//{
//    return  vel > (T)0.0 ?
//        ((T)2.0*f2 + (T)5.0*f1 - f0) / (T)6.0 :
//        ((T)2.0*f1 + (T)5.0*f2 - f3) / (T)6.0;
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T fs_3rd_weno(T f0, T f1, T f2, T f3, T vel)
//{
//    const int   p = 1;
//    const T ep = 1.0e-16;
//
//    const T sgn_u = vel > (T)0.0 ? (T)1.0 : -(T)1.0;
//
//    const T fi = vel > (T)0.0 ? f1 : f2; // 1st-upwind
////  const T fi = vel > 0.0 ? (1.5*f1 - 0.5*f0) : (1.5*f2 - 0.5*f3); // 2nd-upwind
//
//    const T weno_fi[3] = {
//                    (T)1.5*f1 - (T)0.5*f0,
//                    (T)1.5*f2 - (T)0.5*f3,
//                    (T)0.5*f1 + (T)0.5*f2  } ;
//
//    const T weno_fx[3] = {
//                    f1 - f0,
//                    f3 - f2,
//                    f2 - f1  };
//
//    // weno weight //
//    const T ideal_ww[3] = {
//                    ((T)1.0 + sgn_u) / (T)6.0,
//                    ((T)1.0 - sgn_u) / (T)6.0,
//                    (T)2.0/(T)3.0 };
//
//    // smooth indicator //
//    T   ww[3], aa[3];
//    for (int ii=0; ii<3; ii++) {
//        const T IS =
//              pow(weno_fi[ii] - fi, 2)
//            + pow(weno_fx[ii]     , 2) * 4.0/3.0
//            - sgn_u * (weno_fi[ii] - fi) * (weno_fx[ii]);
//
//        aa[ii] = ideal_ww[ii] / (pow(IS, p) + ep);
//    }
//    const T aa_all = aa[0] + aa[1] + aa[2] + ep;
//
//    // weight
//    for (int ii=0; ii<3; ii++) {
//        ww[ii] = aa[ii] / aa_all;
//    }
//
//    return  (ww[0]*weno_fi[0] + ww[1]*weno_fi[1] + ww[2]*weno_fi[2]);
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T fs_4th_central(T f0, T f1, T f2, T f3)
//{
//    return  -(f3 - (T)7.0*f2 - (T)7.0*f1 + f0) / (T)12.0;
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T fs_5th_upwind(T f0, T f1, T f2, T f3, T f4, T f5, T vel)
//{
//    return  vel > (T)0.0 ?
//            ((T)2.0*f0 - (T)13.0*f1 + (T)47.0*f2 + (T)27.0*f3 - (T)3.0*f4) / (T)60.0 :
//            ((T)2.0*f5 - (T)13.0*f4 + (T)47.0*f3 + (T)27.0*f2 - (T)3.0*f1) / (T)60.0 ;
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T fs_5th_weno(T f0, T f1, T f2, T f3, T f4, T f5, T vel)
//{
//    const int   p  = 2;
//    const double ep = 1.0e-32;
//
//    const double sgn_u = vel > 0.0 ? 1.0 : -1.0;
//
//    const double fi = vel > 0.0 ?
//                        (-f1 + 5.0*f2 + 2.0*f3)/6.0 :
//                        (-f4 + 5.0*f3 + 2.0*f2)/6.0; // 3rd-upwind
//
//    // weno //
//    const double weno_fi[4] = {
//                        (  2.0*f0 - 7.0*f1 + 11.0*f2)/6.0,
//                        (-        f1 + 5.0*f2 +  2.0*f3)/6.0,
//                        (  2.0*f2 + 5.0*f3 -         f4)/6.0,
//                        ( 11.0*f3 - 7.0*f4 +  2.0*f5)/6.0 };
//
//    const double weno_fx[4] = {
//                        (      f0 - 3.0*f1 +  2.0*f2),
//                        ( -       f2 +        f3),
//                        ( -       f2 +        f3),
//                        ( -2.0*f3 + 3.0*f4 -         f5) };
//
//    const double weno_fxx[4] = {
//                        ( f0 - 2.0*f1 + f2),
//                        ( f1 - 2.0*f2 + f3),
//                        ( f2 - 2.0*f3 + f4),
//                        ( f3 - 2.0*f4 + f5) };
//
//    // weight //
//    const double ideal_ww[4] = {
//                        (1.0 +        sgn_u) /20.0,
//                        (9.0 + 3.0*sgn_u) /20.0,
//                        (9.0 - 3.0*sgn_u) /20.0,
//                        (1.0 -        sgn_u) /20.0 };
//
//    double   ff     = 0.0;
//    double   aa_all = ep;
//
//    // smooth indicator //
//    for (int ii=0; ii<4; ii++) {
//        const double IS =
//              ( pow(weno_fi[ii]-fi, 2) + pow(weno_fx[ii], 2) + pow(weno_fxx[ii], 2) )
//            - ( 2.0*(weno_fi[ii]-fi)*weno_fx[ii] + 2.0*weno_fx[ii]*weno_fxx[ii] ) /2.0 * sgn_u
//            + ( pow(weno_fx[ii], 2) + weno_fx[ii]*weno_fxx[ii] + pow(weno_fxx[ii], 2) ) /3.0
//            - ( weno_fx[ii]*weno_fxx[ii] ) /4.0 * sgn_u
//            + ( 0.25*pow(weno_fxx[ii], 2) ) /5.0;
//
//        const double tmp = ideal_ww[ii] / (pow(IS, p) + ep);
//
//        ff     += tmp*weno_fi[ii];
//        aa_all += tmp;
//    }
//
//    return  ff / aa_all;
//}


//template<typename T>
//__HOST__ __DEVICE__
//void poly_5th_weno(
//    T& fs, T& fx, T& fxx,
//    T f0, T f1, T f2, T f3, T f4, T f5,
//    T vel, T dx
//    )
//{
//    const int   p  = 2;
//    const T ep = 1.0e-16;
//
//    const T sgn_u = vel > (T)0.0 ? (T)1.0 : -(T)1.0;
//
//    const T fi = vel > (T)0.0 ?
//                        (-f1 + (T)5.0*f2 + (T)2.0*f3)/(T)6.0 :
//                        (-f4 + (T)5.0*f3 + (T)2.0*f2)/(T)6.0; // 3rd-upwind
//
//    // weno //
//    const T weno_fi[4] = {
//                        ( (T) 2.0*f0 - (T)7.0*f1 + (T)11.0*f2)/(T)6.0,
//                        (-        f1 + (T)5.0*f2 + (T) 2.0*f3)/(T)6.0,
//                        ( (T) 2.0*f2 + (T)5.0*f3 -         f4)/(T)6.0,
//                        ( (T)11.0*f3 - (T)7.0*f4 + (T) 2.0*f5)/(T)6.0 };
//
//    const T weno_fx[4] = {
//                        (  (T)    f0 - (T)3.0*f1 +  (T)2.0*f2),
//                        ( -       f2 +        f3),
//                        ( -       f2 +        f3),
//                        ( -(T)2.0*f3 + (T)3.0*f4 -         f5) };
//
//    const T weno_fxx[4] = {
//                        ( f0 - (T)2.0*f1 + f2),
//                        ( f1 - (T)2.0*f2 + f3),
//                        ( f2 - (T)2.0*f3 + f4),
//                        ( f3 - (T)2.0*f4 + f5) };
//
//    // weight //
//    const T ideal_ww[4] = {
//                        ((T)1.0 +        sgn_u) /(T)20.0,
//                        ((T)9.0 + (T)3.0*sgn_u) /(T)20.0,
//                        ((T)9.0 - (T)3.0*sgn_u) /(T)20.0,
//                        ((T)1.0 -        sgn_u) /(T)20.0 };
//
//    T   ff     = 0.0;
//    T   ffx    = 0.0;
//    T   ffxx   = 0.0;
//    T   aa_all = ep;
//
//    // smooth indicator //
//    for (int ii=0; ii<4; ii++) {
//        const T IS =
//              ( pow(weno_fi[ii]-fi, 2) + pow(weno_fx[ii], 2) + pow(weno_fxx[ii], 2) )
//            - ( (T)2.0*(weno_fi[ii]-fi)*weno_fx[ii] + (T)2.0*weno_fx[ii]*weno_fxx[ii] ) /(T)2.0 * sgn_u
//            + ( pow(weno_fx[ii], 2) + weno_fx[ii]*weno_fxx[ii] + pow(weno_fxx[ii], 2) ) /(T)3.0
//            - ( weno_fx[ii]*weno_fxx[ii] ) /(T)4.0 * sgn_u
//            + ( (T)0.25*pow(weno_fxx[ii], 2) ) /(T)5.0;
//
//        const T tmp = ideal_ww[ii] / (pow(IS, p) + ep);
//
//        ff     += tmp*weno_fi [ii];
//        ffx    += tmp*weno_fx [ii];
//        ffxx   += tmp*weno_fxx[ii];
//        aa_all += tmp;
//    }
//
//    const T daa_all = (T)1.0/aa_all;
//    const T ddx = (T)1.0/dx;
//
//    // update
//    fs  = ff  *daa_all;
//    fx  = ffx *daa_all*ddx;
//    fxx = ffxx*daa_all*ddx*ddx;
//}
//
//
//// average //
//template<typename T>
//__HOST__ __DEVICE__
//T average_3stencil(const T *f, int m)
//{
//    return  fs_2nd_central (f[-m], f[0], f[+m]);
//}
//
//
//// fx //
//template<typename T>
//__HOST__ __DEVICE__
//T fx_1st_central(T f0, T f1, T dx)
//{
//    return  (f1 - f0) / dx;
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T fx_1st_upwind(T f0, T f1, T f2, T dx, T vel)
//{
//    return  vel > (T)0.0 ?
//            (f1 - f0) / dx :
//            (f2 - f1) / dx;
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T fx_2nd_central(T f0, T f1, T f2, T dx)
//{
//    return  (f2 - f0) / ((T)2.0*dx);
//}
//
//
//// u fx //
//template<typename T>
//__HOST__ __DEVICE__
//T ufx_1st_central(T u0, T u1, T f0, T f1, T dx)
//{
//    return  fs_1st_central(u0, u1) * fx_1st_central(f0, f1, dx);
//}
//
//
//template<typename T>
//__HOST__ __DEVICE__
//T usfx_1st_central(T u0, T u1, T u2, T f0, T f1, T dx)
//{
//    return  fs_2nd_central(u0, u1, u2) * fx_1st_central(f0, f1, dx);
//}
//
//
//// semi lagrangian //
//template<typename T>
//__HOST__ __DEVICE__
//T fs_3rd_upwind_poly(
//    T& fs,
//    T& fx,
//    T& fxx,
//    T  f0,
//    T  f1,
//    T  f2,
//    T  f3,
//    T  vel,
//    T  dx
//    )
//{
//    if (vel > (T)0.0) {
//        fs  = ((T)2.0*f2 + (T)5.0*f1 - f0) / (T)6.0;
//        fx  = (f2 - f1) / dx;
//        fxx = (f2 - (T)2.0*f1 + f0) / (dx*dx);
//    }
//    else {
//        fs  = ((T)2.0*f1 + (T)5.0*f2 - f3) / (T)6.0;
//        fx  = (f2 - f1) / dx;
//        fxx = (f3 - (T)2.0*f2 + f1) / (dx*dx);
//    }
//}



// taylor expansion //
template<typename T>
__HOST__ __DEVICE__
T taylor_expansion_dt(
    T ft,
    T ftt,
    double dt
    )
{
    return  ft + ftt*dt*(T)0.5;
}


template<typename T>
__HOST__ __DEVICE__
T taylor_expansion_dt(
    T ft,
    T ftt,
    T ft3,
    double dt
    )
{
    return  taylor_expansion_dt (ft, ftt, dt) + ft3*pow(dt, 2)/(T)6.0;
}


template<typename T>
__HOST__ __DEVICE__
T taylor_expansion_dt(
    T ft,
    T ftt,
    T ft3,
    T ft4,
    double dt
    )
{
    return  taylor_expansion_dt(ft, ftt, ft3, dt) + ft4*pow(dt, 3)/(T)24.0;
}


//template<typename T>
//__HOST__ __DEVICE__
//T f_average_w_obj(
//    const int ids[],
//    const T*  val,
//    const T*  lv_obj
//    )
//{
//    constexpr T cF = (T)1000.0;
//    constexpr T cO = (T)1.0;
//
//    T   f      = (real)0.0;
//    T   weight = (real)0.0;
//    T   tmp_weight;
//    for (int kk=-1; kk<=1; kk++) {
//    for (int jj=-1; jj<=1; jj++) {
//    for (int ii=-1; ii<=1; ii++) {
//        const int id_tmp = Index::id(ids, ii, jj, kk);
//        tmp_weight = ( FuncObj::is_fluid( lv_obj[ id_tmp ] ) ) ? cF : cO;
//        f += val[ id_tmp ] * tmp_weight;  weight += tmp_weight;
//    }
//    }
//    }
//
//    return  ( f/weight );
//}


} // namespace //


#endif
