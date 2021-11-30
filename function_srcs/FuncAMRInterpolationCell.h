#pragma once
#ifndef FUNCAMRINTERPOLATIONCELL_H_
#define FUNCAMRINTERPOLATIONCELL_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include "defineAMR.h"
#include "defineLBM.h"
//#include "Array3D.h"
#include "FuncLBMAMR.h"
#include "FuncLBM.h"
#include "FuncCumulantLBM.h"

#include "Index.h"


namespace  FuncAMRInterpolationCell {


template <typename T>
__HOST__ __DEVICE__
T  average8(
    const T*  val,
    const int indices[]
    )
{
    constexpr T inv8 = 1.0/8.0;
    return
    (   val[ indices[0] ] +  val[ indices[1] ] +  val[ indices[2] ] +  val[ indices[3] ]
     +  val[ indices[4] ] +  val[ indices[5] ] +  val[ indices[6] ] +  val[ indices[7] ]  ) * inv8;
}


template <typename T>
__HOST__ __DEVICE__
T  average8A(
    const T*  val0,
    const T*  val1,
    const int indices[]
    )
{
    constexpr T inv16 = 1.0/16.0;
    return
    (   val0[ indices[0] ] +  val0[ indices[1] ] +  val0[ indices[2] ] +  val0[ indices[3] ]
     +  val0[ indices[4] ] +  val0[ indices[5] ] +  val0[ indices[6] ] +  val0[ indices[7] ]
     +  val1[ indices[0] ] +  val1[ indices[1] ] +  val1[ indices[2] ] +  val1[ indices[3] ]
     +  val1[ indices[4] ] +  val1[ indices[5] ] +  val1[ indices[6] ] +  val1[ indices[7] ]  ) * inv16;
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_flat(
    const T*    val,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm  // fdm: 0, lbm: nn_leaf*idv //
    )
{
    const T val0  =   val[ Index::id( ix, jy, kz, offsets ) + stride_lbm ];
    return val0;
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_linear_fvm(
    const T*    val,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    const T val0  =   val[ Index::id( ix, jy, kz, offsets ) + stride_lbm ];
    const T val_x = ( val[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_y = ( val[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_z = ( val[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5;


    const auto val_integral = [&](const T& X, const T& Y, const T& Z){
                    constexpr T D2  = (T)0.5;

                    const T XYZ = X*Y*Z;

                    return        val0 *XYZ
                            + ( X*val_x + Y*val_y + Z*val_z ) *XYZ*D2;
        };

//    constexpr T DX  = (T)0.5;
    constexpr T DX3 = (T)0.125;

    const T X0 = X - (T)0.25;
    const T X1 = X + (T)0.25;

    const T Y0 = Y - (T)0.25;
    const T Y1 = Y + (T)0.25;

    const T Z0 = Z - (T)0.25;
    const T Z1 = Z + (T)0.25;

    return ( - val_integral((T)X0, (T)Y0, (T)Z0) + val_integral((T)X1, (T)Y0, (T)Z0)
             + val_integral((T)X0, (T)Y1, (T)Z0) - val_integral((T)X1, (T)Y1, (T)Z0)
             + val_integral((T)X0, (T)Y0, (T)Z1) - val_integral((T)X1, (T)Y0, (T)Z1)
             - val_integral((T)X0, (T)Y1, (T)Z1) + val_integral((T)X1, (T)Y1, (T)Z1)  ) / (DX3);
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_quad_fvm(
    const T*    val,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    const T val0  =   val[ Index::id( ix, jy, kz, offsets ) + stride_lbm ];
    const T val_x = ( val[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_y = ( val[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_z = ( val[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5;

    const T val_xy = ( ( val[ Index::id( ix+1,  jy+1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy+1,  kz  , offsets ) + stride_lbm ] )
                     - ( val[ Index::id( ix+1,  jy-1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.25;
    const T val_yz = ( ( val[ Index::id( ix  ,  jy+1,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz+1, offsets ) + stride_lbm ] )
                     - ( val[ Index::id( ix  ,  jy+1,  kz-1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25;
    const T val_zx = ( ( val[ Index::id( ix+1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix+1,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                     - ( val[ Index::id( ix-1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25;

    const T val_xyz = ( ( ( val[ Index::id( ix+1,  jy+1,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy+1,  kz+1, offsets ) + stride_lbm ] )
                        - ( val[ Index::id( ix+1,  jy-1,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy-1,  kz+1, offsets ) + stride_lbm ] ) )
                      - ( ( val[ Index::id( ix+1,  jy+1,  kz-1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy+1,  kz-1, offsets ) + stride_lbm ] )
                        - ( val[ Index::id( ix+1,  jy-1,  kz-1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) ) * (T)0.125;


    const auto val_integral = [&](const T& X, const T& Y, const T& Z){
                    constexpr T D2  = (T)0.5;
                    constexpr T D4  = (T)0.25;
                    constexpr T D8  = (T)0.125;

                    const T XYZ = X*Y*Z;

                    return        val0 *XYZ
                            + ( X*val_x + Y*val_y + Z*val_z ) *XYZ*D2
                            + ( X*Y*val_xy + Y*Z*val_yz + Z*X*val_zx ) *XYZ*D4
                            + ( X*Y*Z*val_xyz ) *XYZ*D8;
        };

//    constexpr T DX  = (T)0.5;
//    constexpr T DX3 = (T)0.125;
    constexpr T D_DX3 = (T)8.0;

    const T X0 = X - (T)0.25;
    const T X1 = X + (T)0.25;

    const T Y0 = Y - (T)0.25;
    const T Y1 = Y + (T)0.25;

    const T Z0 = Z - (T)0.25;
    const T Z1 = Z + (T)0.25;

    return ( - val_integral((T)X0, (T)Y0, (T)Z0) + val_integral((T)X1, (T)Y0, (T)Z0)
             + val_integral((T)X0, (T)Y1, (T)Z0) - val_integral((T)X1, (T)Y1, (T)Z0)
             + val_integral((T)X0, (T)Y0, (T)Z1) - val_integral((T)X1, (T)Y0, (T)Z1)
             - val_integral((T)X0, (T)Y1, (T)Z1) + val_integral((T)X1, (T)Y1, (T)Z1)  ) * (D_DX3);
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_quad_fvm_FA(
    const T*    val,
    const T*    valn,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    const T val0  = ( val [ Index::id( ix, jy, kz, offsets ) + stride_lbm ]
                    + valn[ Index::id( ix, jy, kz, offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_x = ( ( val [ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5
                    + ( valn[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5 ) * (T)0.5;
    const T val_y = ( ( val [ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5
                    + ( valn[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5 ) * (T)0.5;
    const T val_z = ( ( val [ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5
                    + ( valn[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5 ) * (T)0.5;

    const T val_xy = ( ( ( val [ Index::id( ix+1,  jy+1,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy+1,  kz  , offsets ) + stride_lbm ] )
                       - ( val [ Index::id( ix+1,  jy-1,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.25
                     + ( ( valn[ Index::id( ix+1,  jy+1,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy+1,  kz  , offsets ) + stride_lbm ] )
                       - ( valn[ Index::id( ix+1,  jy-1,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.25 ) * (T)0.5;
    const T val_yz = ( ( ( val [ Index::id( ix  ,  jy+1,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy-1,  kz+1, offsets ) + stride_lbm ] )
                       - ( val [ Index::id( ix  ,  jy+1,  kz-1, offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25
                     + ( ( valn[ Index::id( ix  ,  jy+1,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy-1,  kz+1, offsets ) + stride_lbm ] )
                       - ( valn[ Index::id( ix  ,  jy+1,  kz-1, offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25 ) * (T)0.5;
    const T val_zx = ( ( ( val [ Index::id( ix+1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix+1,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                       - ( val [ Index::id( ix-1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25
                     + ( ( valn[ Index::id( ix+1,  jy  ,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix+1,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                       - ( valn[ Index::id( ix-1,  jy  ,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25 ) * (T)0.5;

    const T val_xyz = ( ( ( ( val [ Index::id( ix+1,  jy+1,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy+1,  kz+1, offsets ) + stride_lbm ] )
                          - ( val [ Index::id( ix+1,  jy-1,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy-1,  kz+1, offsets ) + stride_lbm ] ) )
                        - ( ( val [ Index::id( ix+1,  jy+1,  kz-1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy+1,  kz-1, offsets ) + stride_lbm ] )
                          - ( val [ Index::id( ix+1,  jy-1,  kz-1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) ) * (T)0.125
                      + ( ( ( valn[ Index::id( ix+1,  jy+1,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy+1,  kz+1, offsets ) + stride_lbm ] )
                          - ( valn[ Index::id( ix+1,  jy-1,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy-1,  kz+1, offsets ) + stride_lbm ] ) )
                        - ( ( valn[ Index::id( ix+1,  jy+1,  kz-1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy+1,  kz-1, offsets ) + stride_lbm ] )
                          - ( valn[ Index::id( ix+1,  jy-1,  kz-1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) ) * (T)0.125 ) * (T)0.5;


    const auto val_integral = [&](const T& X, const T& Y, const T& Z){
                    constexpr T D2  = (T)0.5;
                    constexpr T D4  = (T)0.25;
                    constexpr T D8  = (T)0.125;

                    const T XYZ = X*Y*Z;

                    return        val0 *XYZ
                            + ( X*val_x + Y*val_y + Z*val_z ) *XYZ*D2
                            + ( X*Y*val_xy + Y*Z*val_yz + Z*X*val_zx ) *XYZ*D4
                            + ( X*Y*Z*val_xyz ) *XYZ*D8;
        };

//    constexpr T DX  = (T)0.5;
//    constexpr T DX3 = (T)0.125;
    constexpr T D_DX3 = (T)8.0;

    const T X0 = X - (T)0.25;
    const T X1 = X + (T)0.25;

    const T Y0 = Y - (T)0.25;
    const T Y1 = Y + (T)0.25;

    const T Z0 = Z - (T)0.25;
    const T Z1 = Z + (T)0.25;

    return ( - val_integral((T)X0, (T)Y0, (T)Z0) + val_integral((T)X1, (T)Y0, (T)Z0)
             + val_integral((T)X0, (T)Y1, (T)Z0) - val_integral((T)X1, (T)Y1, (T)Z0)
             + val_integral((T)X0, (T)Y0, (T)Z1) - val_integral((T)X1, (T)Y0, (T)Z1)
             - val_integral((T)X0, (T)Y1, (T)Z1) + val_integral((T)X1, (T)Y1, (T)Z1)  ) * (D_DX3);
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_cubic_fvm(
    const T*    val,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    const T val0  =   val[ Index::id( ix, jy, kz, offsets ) + stride_lbm ];
    const T val_x = ( val[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_y = ( val[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_z = ( val[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5;

    const T val_xy = ( ( val[ Index::id( ix+1,  jy+1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy+1,  kz  , offsets ) + stride_lbm ] )
                     - ( val[ Index::id( ix+1,  jy-1,  kz  , offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.25;
    const T val_yz = ( ( val[ Index::id( ix  ,  jy+1,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz+1, offsets ) + stride_lbm ] )
                     - ( val[ Index::id( ix  ,  jy+1,  kz-1, offsets ) + stride_lbm ] - val[ Index::id( ix  ,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25;
    const T val_zx = ( ( val[ Index::id( ix+1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix+1,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                     - ( val[ Index::id( ix-1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25;

    const T val_xyz = ( ( ( val[ Index::id( ix+1,  jy+1,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy+1,  kz+1, offsets ) + stride_lbm ] )
                        - ( val[ Index::id( ix+1,  jy-1,  kz+1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy-1,  kz+1, offsets ) + stride_lbm ] ) )
                      - ( ( val[ Index::id( ix+1,  jy+1,  kz-1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy+1,  kz-1, offsets ) + stride_lbm ] )
                        - ( val[ Index::id( ix+1,  jy-1,  kz-1, offsets ) + stride_lbm ] - val[ Index::id( ix-1,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) ) * (T)0.125;

    const T val_xx = ( val[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - (T)2.0 * val[ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + val[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] );
    const T val_yy = ( val[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - (T)2.0 * val[ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + val[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] );
    const T val_zz = ( val[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - (T)2.0 * val[ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + val[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] );


    const auto val_integral = [&](const T& X, const T& Y, const T& Z){
                    constexpr T D2  = (T)0.5;
                    constexpr T D4  = (T)0.25;
                    constexpr T D8  = (T)0.125;
                    constexpr T D3  = (T)0.333333333333333333;

                    const T XYZ = X*Y*Z;

                    return        val0 *XYZ
                            + ( X*val_x + Y*val_y + Z*val_z ) *XYZ*D2
                            + ( X*Y*val_xy + Y*Z*val_yz + Z*X*val_zx ) *XYZ*D4
                            + ( X*Y*Z*val_xyz ) *XYZ*D8
                            + ( X*X*val_xx + Y*Y*val_yy + Z*Z*val_zz )*(T)0.5 *XYZ*D3;
        };

//    constexpr T DX  = (T)0.5;
    constexpr T DX3 = (T)0.125;

    const T X0 = X - (T)0.25;
    const T X1 = X + (T)0.25;

    const T Y0 = Y - (T)0.25;
    const T Y1 = Y + (T)0.25;

    const T Z0 = Z - (T)0.25;
    const T Z1 = Z + (T)0.25;

    return ( - val_integral((T)X0, (T)Y0, (T)Z0) + val_integral((T)X1, (T)Y0, (T)Z0)
             + val_integral((T)X0, (T)Y1, (T)Z0) - val_integral((T)X1, (T)Y1, (T)Z0)
             + val_integral((T)X0, (T)Y0, (T)Z1) - val_integral((T)X1, (T)Y0, (T)Z1)
             - val_integral((T)X0, (T)Y1, (T)Z1) + val_integral((T)X1, (T)Y1, (T)Z1)  ) / (DX3);
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_cubic_fvm_FA(
    const T*    val,
    const T*    valn,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    const T val0  = ( val [ Index::id( ix, jy, kz, offsets ) + stride_lbm ]
                    + valn[ Index::id( ix, jy, kz, offsets ) + stride_lbm ] ) * (T)0.5;
    const T val_x = ( ( val [ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5
                    + ( valn[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) * (T)0.5 ) * (T)0.5;
    const T val_y = ( ( val [ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5
                    + ( valn[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) * (T)0.5 ) * (T)0.5;
    const T val_z = ( ( val [ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5
                    + ( valn[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) * (T)0.5 ) * (T)0.5;

    const T val_xy = ( ( ( val [ Index::id( ix+1,  jy+1,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy+1,  kz  , offsets ) + stride_lbm ] )
                       - ( val [ Index::id( ix+1,  jy-1,  kz  , offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.25
                     + ( ( valn[ Index::id( ix+1,  jy+1,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy+1,  kz  , offsets ) + stride_lbm ] )
                       - ( valn[ Index::id( ix+1,  jy-1,  kz  , offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.25 ) * (T)0.5;
    const T val_yz = ( ( ( val [ Index::id( ix  ,  jy+1,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy-1,  kz+1, offsets ) + stride_lbm ] )
                       - ( val [ Index::id( ix  ,  jy+1,  kz-1, offsets ) + stride_lbm ] - val [ Index::id( ix  ,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25
                     + ( ( valn[ Index::id( ix  ,  jy+1,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy-1,  kz+1, offsets ) + stride_lbm ] )
                       - ( valn[ Index::id( ix  ,  jy+1,  kz-1, offsets ) + stride_lbm ] - valn[ Index::id( ix  ,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25 ) * (T)0.5;
    const T val_zx = ( ( ( val [ Index::id( ix+1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix+1,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                       - ( val [ Index::id( ix-1,  jy  ,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25
                     + ( ( valn[ Index::id( ix+1,  jy  ,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix+1,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                       - ( valn[ Index::id( ix-1,  jy  ,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.25 ) * (T)0.5;

    const T val_xyz = ( ( ( ( val [ Index::id( ix+1,  jy+1,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy+1,  kz+1, offsets ) + stride_lbm ] )
                          - ( val [ Index::id( ix+1,  jy-1,  kz+1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy-1,  kz+1, offsets ) + stride_lbm ] ) )
                        - ( ( val [ Index::id( ix+1,  jy+1,  kz-1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy+1,  kz-1, offsets ) + stride_lbm ] )
                          - ( val [ Index::id( ix+1,  jy-1,  kz-1, offsets ) + stride_lbm ] - val [ Index::id( ix-1,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) ) * (T)0.125
                      + ( ( ( valn[ Index::id( ix+1,  jy+1,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy+1,  kz+1, offsets ) + stride_lbm ] )
                          - ( valn[ Index::id( ix+1,  jy-1,  kz+1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy-1,  kz+1, offsets ) + stride_lbm ] ) )
                        - ( ( valn[ Index::id( ix+1,  jy+1,  kz-1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy+1,  kz-1, offsets ) + stride_lbm ] )
                          - ( valn[ Index::id( ix+1,  jy-1,  kz-1, offsets ) + stride_lbm ] - valn[ Index::id( ix-1,  jy-1,  kz-1, offsets ) + stride_lbm ] ) ) ) * (T)0.125 ) * (T)0.5;


    const T val_xx = (  ( val [ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - (T)2.0 * val [ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + val [ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] )
                      + ( valn[ Index::id( ix+1,  jy  ,  kz  , offsets ) + stride_lbm ] - (T)2.0 * valn[ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + valn[ Index::id( ix-1,  jy  ,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.5;
    const T val_yy = (  ( val [ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - (T)2.0 * val [ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + val [ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] )
                      + ( valn[ Index::id( ix  ,  jy+1,  kz  , offsets ) + stride_lbm ] - (T)2.0 * valn[ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + valn[ Index::id( ix  ,  jy-1,  kz  , offsets ) + stride_lbm ] ) ) * (T)0.5;
    const T val_zz = (  ( val [ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - (T)2.0 * val [ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + val [ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] )
                      + ( valn[ Index::id( ix  ,  jy  ,  kz+1, offsets ) + stride_lbm ] - (T)2.0 * valn[ Index::id( ix,  jy,  kz, offsets ) + stride_lbm ] + valn[ Index::id( ix  ,  jy  ,  kz-1, offsets ) + stride_lbm ] ) ) * (T)0.5;


    const auto val_integral = [&](const T& X, const T& Y, const T& Z){
                    constexpr T D2  = (T)0.5;
                    constexpr T D4  = (T)0.25;
                    constexpr T D8  = (T)0.125;
                    constexpr T D3  = (T)0.333333333333333333;

                    const T XYZ = X*Y*Z;

                    return        val0 *XYZ
                            + ( X*val_x + Y*val_y + Z*val_z ) *XYZ*D2
                            + ( X*Y*val_xy + Y*Z*val_yz + Z*X*val_zx ) *XYZ*D4
                            + ( X*Y*Z*val_xyz ) *XYZ*D8
                            + ( X*X*val_xx + Y*Y*val_yy + Z*Z*val_zz )*(T)0.5 *XYZ*D3;
        };

//    constexpr T DX  = (T)0.5;
//    constexpr T DX3 = (T)0.125;
    constexpr T D_DX3 = (T)8.0;

    const T X0 = X - (T)0.25;
    const T X1 = X + (T)0.25;

    const T Y0 = Y - (T)0.25;
    const T Y1 = Y + (T)0.25;

    const T Z0 = Z - (T)0.25;
    const T Z1 = Z + (T)0.25;

    return ( - val_integral((T)X0, (T)Y0, (T)Z0) + val_integral((T)X1, (T)Y0, (T)Z0)
             + val_integral((T)X0, (T)Y1, (T)Z0) - val_integral((T)X1, (T)Y1, (T)Z0)
             + val_integral((T)X0, (T)Y0, (T)Z1) - val_integral((T)X1, (T)Y0, (T)Z1)
             - val_integral((T)X0, (T)Y1, (T)Z1) + val_integral((T)X1, (T)Y1, (T)Z1)  ) * (D_DX3);
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell(
    const T*    val,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    // fvm //
//    const T val_int = Val_interpolation_cell_linear_fvm(val, ix,jy,kz, offsets, stride_lbm, X,Y,Z);
    const T val_int = Val_interpolation_cell_quad_fvm(val, ix,jy,kz, offsets, stride_lbm, X,Y,Z);
//    const T val_int = Val_interpolation_cell_cubic_fvm(val, ix,jy,kz, offsets, stride_lbm, X,Y,Z);

    return val_int;
}


template <typename T>
__HOST__ __DEVICE__
T  Val_interpolation_cell_FA(
    const T*    val,
    const T*    valn,
    const int   ix,
    const int   jy,
    const int   kz,
    const int   offsets[],
    const int   stride_lbm, // fdm: 0, lbm: nn_leaf*idv //
    const T     X,
    const T     Y,
    const T     Z
    )
{
    // fvm //
    const T val_int = Val_interpolation_cell_quad_fvm_FA(val, valn, ix,jy,kz, offsets, stride_lbm, X,Y,Z);
//    const T val_int = Val_interpolation_cell_cubic_fvm_FA(val, valn, ix,jy,kz, offsets, stride_lbm, X,Y,Z);

//    const T val_int = ( Val_interpolation_cell_linear_fvm(val , ix,jy,kz, offsets, stride_lbm, X,Y,Z)
//                         + Val_interpolation_cell_linear_fvm(valn, ix,jy,kz, offsets, stride_lbm, X,Y,Z) ) * (T)0.5;

//    const T val_int = ( Val_interpolation_cell_quad_fvm(val , ix,jy,kz, offsets, stride_lbm, X,Y,Z)
//                         + Val_interpolation_cell_quad_fvm(valn, ix,jy,kz, offsets, stride_lbm, X,Y,Z) ) * (T)0.5;

    return val_int;
}


template <typename T>
__HOST__ __DEVICE__
void  L2C_kernel(
    const int   i,
    const int   j,
    const int   k,
          T*    valC,
    const T*    val,
    const int   offsetC,
    const int   offsets[],
    const int   lv_offset[], // [3] : 0 or 1 (left or right side in a coarse mesh) //
    const int   lvL
    )
{
    constexpr int  nx_leaf   = DefAMR::NX_LEAF;
    constexpr int  nx_leafC  = nx_leaf/2;


    // i,j,k    : 0,nx_leafC-1 //
    // coarse //
    const int  ix_stC[3] = { nx_leafC*lv_offset[0], nx_leafC*lv_offset[1], nx_leafC*lv_offset[2] };
    const int  ixC = i + ix_stC[0];
    const int  iyC = j + ix_stC[1];
    const int  izC = k + ix_stC[2];

    const int  idC = Index::id(ixC,iyC,izC) +  offsetC;

    // this //
    const int  ix = i*2;
    const int  jy = j*2;
    const int  kz = k*2;

    const int  id8[8] = {  Index::id(  ix  ,  jy  ,  kz  , offsets  ),
                           Index::id(  ix+1,  jy  ,  kz  , offsets  ),
                           Index::id(  ix  ,  jy+1,  kz  , offsets  ),
                           Index::id(  ix+1,  jy+1,  kz  , offsets  ),
                           Index::id(  ix  ,  jy  ,  kz+1, offsets  ),
                           Index::id(  ix+1,  jy  ,  kz+1, offsets  ),
                           Index::id(  ix  ,  jy+1,  kz+1, offsets  ),
                           Index::id(  ix+1,  jy+1,  kz+1, offsets  )  };


    // update //
    valC[idC] = average8(val, id8);
}


template <typename T0, typename T1>
__HOST__ __DEVICE__
void  LBM_L2C_kernel(
    const int  i,
    const int  j,
    const int  k,
          T0*  valC,
    const T0*  val,
    const int  offsetC,
    const int  offsets[],
    const int  lv_offset[], // [3] : 0 or 1 (left or right side in a coarse mesh) //
    const int  lvL,
    const T1   vis,
    const T1   c_ref,
    const T1   dt0
    )
{
    constexpr int  nx_leaf  = DefAMR::NX_LEAF;
    constexpr int  nn_leaf  = nx_leaf*nx_leaf*nx_leaf;
    constexpr int  nx_leafC = nx_leaf/2;
    constexpr int  nQ = LBM_velocity_model::nQ;


    // i,j,k    : 0,nx_leafC-1 //
    // coarse //
    const int  ix_stC[3] = { nx_leafC*lv_offset[0], nx_leafC*lv_offset[1], nx_leafC*lv_offset[2] };
    const int  ixC = i + ix_stC[0];
    const int  iyC = j + ix_stC[1];
    const int  izC = k + ix_stC[2];

    const int  idC = Index::id(ixC,iyC,izC) +  offsetC;

    // this //
    const int  ix = i*2;
    const int  jy = j*2;
    const int  kz = k*2;

    const int  id8[8] = {  Index::id(  ix  ,  jy  ,  kz  , offsets  ),
                           Index::id(  ix+1,  jy  ,  kz  , offsets  ),
                           Index::id(  ix  ,  jy+1,  kz  , offsets  ),
                           Index::id(  ix+1,  jy+1,  kz  , offsets  ),
                           Index::id(  ix  ,  jy  ,  kz+1, offsets  ),
                           Index::id(  ix+1,  jy  ,  kz+1, offsets  ),
                           Index::id(  ix  ,  jy+1,  kz+1, offsets  ),
                           Index::id(  ix+1,  jy+1,  kz+1, offsets  )  };

    T0  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        fs[idv] = (   val[ id8[0] + nn_leaf*idv ] +  val[ id8[1] + nn_leaf*idv ] +  val[ id8[2] + nn_leaf*idv ] +  val[ id8[3] + nn_leaf*idv ]
                   +  val[ id8[4] + nn_leaf*idv ] +  val[ id8[5] + nn_leaf*idv ] +  val[ id8[6] + nn_leaf*idv ] +  val[ id8[7] + nn_leaf*idv ]
                  ) * (T0)0.125;
    }


    // cal //
    constexpr T1 dtX_dtL = (T1)2.0;
    const T1 dtL = dt0 * powf((T1)0.5, (T1)lvL);
    const T1 dtX = dtL*(T1)2.0;
    const T1 tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const T1 tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );

#if 1 // srt //
    const T1 coefAMR = dtX_dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 coefAMR = dtX/dtL * (tauX)/(tauL);

    T0  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );
#else // cumulant //
    const T1 omega = (T1)1.0 - dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 omega = (T1)1.0 - dtX/dtL * (tauX)/(tauL);

    T0  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    T0  fsAMR[27];
    FuncCumulantLBM::fs_cumulant_lbm(fsAMR, fs, omega, rhos, us,vs,ws);
#endif


    // update //
    for (int idv=0; idv<nQ; idv++) {
        valC[idC + nn_leaf*idv] = fsAMR[idv];
    }

}


template <typename T0, typename T1>
__HOST__ __DEVICE__
void  LBM_L2CA_kernel(
    const int  i,
    const int  j,
    const int  k,
          T0*  valC,
          T0*  valCn,
    const T0*  val,
    const int  offsetC,
    const int  offsets[],
    const int  lv_offset[], // [3] : 0 or 1 (left or right side in a coarse mesh) //
    const int  lvL,
    const T1   vis,
    const T1   c_ref,
    const T1   dt0
    )
{
    constexpr int  nx_leaf  = DefAMR::NX_LEAF;
    constexpr int  nn_leaf  = nx_leaf*nx_leaf*nx_leaf;
    constexpr int  nx_leafC = nx_leaf/2;
    constexpr int  nQ = LBM_velocity_model::nQ;


    // i,j,k    : 0,nx_leafC-1 //
    // coarse //
    const int  ix_stC[3] = { nx_leafC*lv_offset[0], nx_leafC*lv_offset[1], nx_leafC*lv_offset[2] };
    const int  ixC = i + ix_stC[0];
    const int  iyC = j + ix_stC[1];
    const int  izC = k + ix_stC[2];

    const int  idC = Index::id(ixC,iyC,izC) +  offsetC;

    // this //
    const int  ix = i*2;
    const int  jy = j*2;
    const int  kz = k*2;

    const int  id8[8] = {  Index::id(  ix  ,  jy  ,  kz  , offsets  ),
                           Index::id(  ix+1,  jy  ,  kz  , offsets  ),
                           Index::id(  ix  ,  jy+1,  kz  , offsets  ),
                           Index::id(  ix+1,  jy+1,  kz  , offsets  ),
                           Index::id(  ix  ,  jy  ,  kz+1, offsets  ),
                           Index::id(  ix+1,  jy  ,  kz+1, offsets  ),
                           Index::id(  ix  ,  jy+1,  kz+1, offsets  ),
                           Index::id(  ix+1,  jy+1,  kz+1, offsets  )  };

    T0  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        fs[idv] = (   val[ id8[0] + nn_leaf*idv ] +  val[ id8[1] + nn_leaf*idv ] +  val[ id8[2] + nn_leaf*idv ] +  val[ id8[3] + nn_leaf*idv ]
                   +  val[ id8[4] + nn_leaf*idv ] +  val[ id8[5] + nn_leaf*idv ] +  val[ id8[6] + nn_leaf*idv ] +  val[ id8[7] + nn_leaf*idv ]
                  ) * (T0)0.125;
    }


    // cal //
//    const T1 dtL = dt0 * pow(0.5, lvL);
    const T1 dtL = dt0 * powf((T1)0.5, (T1)lvL);
    const T1 dtX = dtL*(T1)2.0;
    const T1 tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const T1 tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );

#if 1 // srt //
    const T1 coefAMR = dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 coefAMR = dtX/dtL * (tauX)/(tauL);

    T0  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );
#else // cumulant //
    const T1 omega = 1.0 - dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 omega = 1.0 - dtX/dtL * (tauX)/(tauL);

    T0  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    T0  fsAMR[27];
    FuncCumulantLBM::fs_cumulant_lbm(fsAMR, fs, omega, rhos, us,vs,ws);
#endif


    // update //
    for (int idv=0; idv<nQ; idv++) {
        valC [idC + nn_leaf*idv] = fsAMR[idv];
        valCn[idC + nn_leaf*idv] = fsAMR[idv];
    }

}


template <typename T>
__HOST__ __DEVICE__
void  L2F_kernel(
    const int   i,
    const int   j,
    const int   k,
    const int   ii, // ii,jj,kk : 0,1 //
    const int   jj,
    const int   kk,
          T*    valF,
    const T*    val,
    const int   offsetsF[],
    const int   offsets[],
    const int   lvL
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);

    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix = ii*nx_leafC + (int)(i/2);
    const int  jy = jj*nx_leafC + (int)(j/2);
    const int  kz = kk*nx_leafC + (int)(k/2);

    const T  xx = (i%2 ==0) ? -0.25 : 0.25;
    const T  yy = (j%2 ==0) ? -0.25 : 0.25;
    const T  zz = (k%2 ==0) ? -0.25 : 0.25;

    const T val_int = Val_interpolation_cell(val, ix,jy,kz, offsets, 0, xx,yy,zz);

    // update //
    valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] ] = val_int;
}


template <typename T>
__HOST__ __DEVICE__
void  L2FA_kernel(
    const int   i,
    const int   j,
    const int   k,
    const int   ii,
    const int   jj,
    const int   kk,
          T*    valF,
    const T*    val,
    const T*    valn,
    const int   offsetsF[],
    const int   offsets[],
    const int   lvL
    )
{
    constexpr int  nx_leaf  = DefAMR::NX_LEAF;
    constexpr int  nx_leafC = (int)(nx_leaf/2);

    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix = ii*nx_leafC + (int)(i/2);
    const int  jy = jj*nx_leafC + (int)(j/2);
    const int  kz = kk*nx_leafC + (int)(k/2);

    const T  xx = (i%2 ==0) ? -0.25 : 0.25;
    const T  yy = (j%2 ==0) ? -0.25 : 0.25;
    const T  zz = (k%2 ==0) ? -0.25 : 0.25;

    // linear interpolation //
    const T  val_int = Val_interpolation_cell_FA(val,valn, ix,jy,kz, offsets, 0, xx,yy,zz);

    // update //
    valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] ] = val_int;
}


template <typename T0, typename T1>
__HOST__ __DEVICE__
void  LBM_L2F_kernel(
    const int  i,
    const int  j,
    const int  k,
    const int  ii,
    const int  jj,
    const int  kk,
          T0*  valF,
    const T0*  val,
    const int  offsetsF[],
    const int  offsets[],
    const int  lvL,
    const T1   vis,
    const T1   c_ref,
    const T1   dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf = DefAMR::NN_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);
    constexpr int  nQ = LBM_velocity_model::nQ;

    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix = ii*nx_leafC + (int)(i/2);
    const int  jy = jj*nx_leafC + (int)(j/2);
    const int  kz = kk*nx_leafC + (int)(k/2);

    const T0 xx = (i%2 ==0) ? -0.25 : 0.25;
    const T0 yy = (j%2 ==0) ? -0.25 : 0.25;
    const T0 zz = (k%2 ==0) ? -0.25 : 0.25;

    T0  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        const T0 val_int = Val_interpolation_cell(val, ix,jy,kz, offsets, nn_leaf*idv, xx,yy,zz);

        fs[idv] = val_int;;
    }


    // cal //
    const T1 dtL = dt0 * powf((T1)0.5, (T1)lvL);
    const T1 dtX = dtL*(T1)0.5;
    const T1 tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const T1 tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );

#if 1 // srt //
    const T1 coefAMR = dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 coefAMR = dtX/dtL * (tauX)/(tauL);

    T0  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );
#else // cumulant //
    const T1 omega = 1.0 - dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 omega = 1.0 - dtX/dtL * (tauX)/(tauL);

    T0  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    T0  fsAMR[27];
    FuncCumulantLBM::fs_cumulant_lbm(fsAMR, fs, omega, rhos, us,vs,ws);
#endif

    // update //
    for (int idv=0; idv<nQ; idv++) {
        valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] + nn_leaf*idv ] = fsAMR[idv];
    }
}


template <typename T0, typename T1>
__HOST__ __DEVICE__
void  LBM_L2FA_kernel(
    const int  i,
    const int  j,
    const int  k,
    const int  ii,
    const int  jj,
    const int  kk,
          T0*  valF,
    const T0*  val,
    const T0*  valn,
    const int  offsetsF[],
    const int  offsets[],
    const int  lvL,
    const T1   vis,
    const T1   c_ref,
    const T1   dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf = DefAMR::NN_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);
    constexpr int  nQ = LBM_velocity_model::nQ;

    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix = ii*nx_leafC + (int)(i/2);
    const int  jy = jj*nx_leafC + (int)(j/2);
    const int  kz = kk*nx_leafC + (int)(k/2);

    const T0 xx = (i%2 ==0) ? -0.25 : 0.25;
    const T0 yy = (j%2 ==0) ? -0.25 : 0.25;
    const T0 zz = (k%2 ==0) ? -0.25 : 0.25;

    T0  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        const T0 val_int = Val_interpolation_cell_FA(val,valn, ix,jy,kz, offsets, nn_leaf*idv, xx,yy,zz);

        fs[idv] = val_int;;
    }


    // cal //
    const T1 dtL = dt0 * powf((T1)0.5, (T1)lvL);
    const T1 dtX = dtL*(T1)0.5;
    const T1 tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const T1 tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );

#if 1 // srt //
    const T1 coefAMR = dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 coefAMR = dtX/dtL * (tauX)/(tauL);

    T0  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );
#else // cumulant //
    const T1 omega = (T1)1.0 - dtX/dtL * (tauX - (T1)1.0)/(tauL - (T1)1.0);
//    const T1 omega = (T1)1.0 - dtX/dtL * (tauX)/(tauL);

    T0  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    T0  fsAMR[27];
    FuncCumulantLBM::fs_cumulant_lbm(fsAMR, fs, omega, rhos, us,vs,ws);
#endif

    // update //
    for (int idv=0; idv<nQ; idv++) {
        valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] + nn_leaf*idv ] = fsAMR[idv];
    }

}


};


#endif
