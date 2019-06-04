#pragma once
#ifndef FUNCAMRINTERPOLATIONNODE_H_
#define FUNCAMRINTERPOLATIONNODE_H_


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


namespace  FuncAMRInterpolationNode {


template <typename T>
__HOST__ __DEVICE__
T  average8(
    const T*  val,
    const int indices[]
    )
{
    return
    (   val[ indices[0] ] +  val[ indices[1] ] +  val[ indices[2] ] +  val[ indices[3] ]
     +  val[ indices[4] ] +  val[ indices[5] ] +  val[ indices[6] ] +  val[ indices[7] ]  ) / ((T)8);
}


template <typename T>
__HOST__ __DEVICE__
T  average8A(
    const T*  val0,
    const T*  val1,
    const int indices[]
    )
{
    return
    (   val0[ indices[0] ] +  val0[ indices[1] ] +  val0[ indices[2] ] +  val0[ indices[3] ]
     +  val0[ indices[4] ] +  val0[ indices[5] ] +  val0[ indices[6] ] +  val0[ indices[7] ]
     +  val1[ indices[0] ] +  val1[ indices[1] ] +  val1[ indices[2] ] +  val1[ indices[3] ]
     +  val1[ indices[4] ] +  val1[ indices[5] ] +  val1[ indices[6] ] +  val1[ indices[7] ]  ) / ((T)16);
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
    const int   offset,
    const int   lv_offset[], // [3] : 0 or 1 (left or right side in a coarse mesh) //
    const int   lvL
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
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
    const int  iy = j*2;
    const int  iz = k*2;

    const int  id = Index::id(ix,iy,iz) +  offset;


    // update //
    valC[idC] = val[id];

}


template <typename T>
__HOST__ __DEVICE__
void  LBM_L2C_kernel(
    const int   i,
    const int   j,
    const int   k,
          T*    valC,
    const T*    val,
    const int   offsetC,
    const int   offset,
    const int   lv_offset[], // [3] : 0 or 1 (left or right side in a coarse mesh) //
    const int   lvL,
    const real  vis,
    const real  c_ref,
    const real  dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf   = nx_leaf*nx_leaf*nx_leaf;
    constexpr int  nx_leafC  = nx_leaf/2;
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
    const int  iy = j*2;
    const int  iz = k*2;

    const int  id = Index::id(ix,iy,iz) +  offset;

    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        fs[idv] = val[id + nn_leaf*idv];
    }


    // cal //
//    const real dtL = dt0 * pow(0.5, lvL);
    const real dtL = dt0 * powf((real)0.5, (real)lvL);
    const real dtX = dtL*(real)2.0;
    const real tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const real tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );
//    const real coefAMR = dtX/dtL * (tauX - (real)1.0)/(tauL - (real)1.0);
    const real coefAMR = dtX/dtL * (tauX)/(tauL);

    real  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );


    // update //
    for (int idv=0; idv<nQ; idv++) {
        valC[idC + nn_leaf*idv] = fsAMR[idv];
    }

}


template <typename T>
__HOST__ __DEVICE__
void  LBM_L2C_cumulant_kernel(
    const int   i,
    const int   j,
    const int   k,
          T*    valC,
    const T*    val,
    const int   offsetC,
    const int   offset,
    const int   lv_offset[], // [3] : 0 or 1 (left or right side in a coarse mesh) //
    const int   lvL,
    const real  vis,
    const real  c_ref,
    const real  dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf   = nx_leaf*nx_leaf*nx_leaf;
    constexpr int  nx_leafC  = nx_leaf/2;
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
    const int  iy = j*2;
    const int  iz = k*2;

    const int  id = Index::id(ix,iy,iz) +  offset;

    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        fs[idv] = val[id + nn_leaf*idv];
    }


    // cal //
//    const real dtL = dt0 * pow(0.5, lvL);
    const real dtL = dt0 * powf((real)0.5, (real)lvL);
    const real dtX = dtL*(real)2.0;
    const real tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const real tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );
//    const real coefAMR = dtX/dtL * (tauX - (real)1.0)/(tauL - (real)1.0);
//    const real coefAMR = dtX/dtL * (tauX)/(tauL);

    const real omega = 1.0 - dtX/dtL * (tauX)/(tauL); 

    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    real  fIs[27];
    FuncCumulantLBM::fs_cumulant_lbm(fIs, fs, omega, rhos, us,vs,ws);


    // update //
    for (int idv=0; idv<nQ; idv++) {
        valC[idC + nn_leaf*idv] = fIs[idv];
    }

}


template <typename T>
__HOST__ __DEVICE__
void  L2F_kernel(
    const int   i,
    const int   j,
    const int   k,
    const int   ii,
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
    const int  ix_C[8] = {  Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ) };

    // update //
    valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] ] = average8( val, ix_C );
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
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);

    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix_C[8] = {  Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ) };

    // update //
    valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] ] = ( average8( val, ix_C ) + average8( valn, ix_C ) ) * (real)0.5;
}


template <typename T>
__HOST__ __DEVICE__
void  LBM_L2F_kernel(
    const int   i,
    const int   j,
    const int   k,
    const int   ii,
    const int   jj,
    const int   kk,
          T*    valF,
    const T*    val,
    const int   offsetsF[],
    const int   offsets[],
    const int   lvL,
    const real  vis,
    const real  c_ref,
    const real  dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf = DefAMR::NN_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);
    constexpr int  nQ = LBM_velocity_model::nQ;

#if 1
    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix_C[8] = {  Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) };

    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) { fs[idv] = average8( &val[nn_leaf*idv], ix_C ); }
#else
    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        const int _id = nn_leaf*idv;
        fs[idv] = (   val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                  ) * (real)0.125;
    }
#endif


    // cal //
//    const real dtL = dt0 * pow((real)0.5, lvL);
    const real dtL = dt0 * powf((real)0.5, (real)lvL);
    const real dtX = dtL*(real)0.5;
    const real tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const real tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );
//    const real coefAMR = dtX/dtL * (tauX - (real)1.0)/(tauL - (real)1.0);
    const real coefAMR = dtX/dtL * (tauX)/(tauL);

    real  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );

    // update //
    for (int idv=0; idv<nQ; idv++) {
        valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] + nn_leaf*idv ] = fsAMR[idv];
    }
}


template <typename T>
__HOST__ __DEVICE__
void  LBM_L2F_cumulant_kernel(
    const int   i,
    const int   j,
    const int   k,
    const int   ii,
    const int   jj,
    const int   kk,
          T*    valF,
    const T*    val,
    const int   offsetsF[],
    const int   offsets[],
    const int   lvL,
    const real  vis,
    const real  c_ref,
    const real  dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf = DefAMR::NN_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);
    constexpr int  nQ = LBM_velocity_model::nQ;

#if 1
    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix_C[8] = {  Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) };

    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) { fs[idv] = average8( &val[nn_leaf*idv], ix_C ); }
#else
    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        const int _id = nn_leaf*idv;
        fs[idv] = (   val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                    + val[ _id + Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2),  offsets  ) ] 
                  ) * (real)0.125;
    }
#endif


    // cal //
//    const real dtL = dt0 * pow((real)0.5, lvL);
    const real dtL = dt0 * powf((real)0.5, (real)lvL);
    const real dtX = dtL*(real)0.5;
    const real tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const real tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );
//    const real coefAMR = dtX/dtL * (tauX - (real)1.0)/(tauL - (real)1.0);
//    const real coefAMR = dtX/dtL * (tauX)/(tauL);

    const real omega = 1.0 - dtX/dtL * (tauX)/(tauL); 

    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    real  fIs[27];
    FuncCumulantLBM::fs_cumulant_lbm(fIs, fs, omega, rhos, us,vs,ws);

    // update //
    for (int idv=0; idv<nQ; idv++) {
        valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] + nn_leaf*idv ] = fIs[idv];
    }
}


template <typename T>
__HOST__ __DEVICE__
void  LBM_L2FA_kernel(
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
    const int   lvL,
    const real  vis,
    const real  c_ref,
    const real  dt0
    )
{
    constexpr int  nx_leaf = DefAMR::NX_LEAF;
    constexpr int  nn_leaf = DefAMR::NN_LEAF;
    constexpr int  nx_leafC  = (int)(nx_leaf/2);
    constexpr int  nQ = LBM_velocity_model::nQ;

    // ii,jj,kk : 0,1 //
    // i,j,k    : 0,nx_leaf-1 //
    const int  ix_C[8] = {  Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k  )/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j  )/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i  )/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ),
                            Index::id(  ii*nx_leafC + (int)((i+1)/2),  jj*nx_leafC + (int)((j+1)/2),  kk*nx_leafC + (int)((k+1)/2), offsets  ) };

    real  fs[nQ];
    for (int idv=0; idv<nQ; idv++) {
        fs[idv] = average8A( &val[nn_leaf*idv], &valn[nn_leaf*idv], ix_C );
    }


    // cal //
//    const real dtL = dt0 * pow(0.5, lvL);
    const real dtL = dt0 * powf((real)0.5, (real)lvL);
    const real dtX = dtL*(real)0.5;
    const real tauL = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtL) );
    const real tauX = FuncLBM::relaxation_time( FuncLBM::kvis_lbm(vis, c_ref, dtX) );

//    const real coefAMR = dtX/dtL * (tauX - (real)1.0)/(tauL - (real)1.0);
    const real coefAMR = dtX/dtL * (tauX)/(tauL);

    real  fsAMR[nQ];
    FuncLBMAMR::feqAMR_L2X(
        fsAMR,
        fs,
        coefAMR
        );

    // update //
    for (int idv=0; idv<nQ; idv++) {
        valF[ Index::id( i,j,k ) + offsetsF[Index::idv(ii,jj,kk)] + nn_leaf*idv ] = fsAMR[idv];
    }

}


};


#endif
