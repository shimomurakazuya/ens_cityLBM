#include "LBMCalculationGPU.h"
#include "FuncLBMKernel.h"
#include "FuncAMRInterpolationCell.h"
#include "FuncScalarAdvection.h"
#include "FuncThermalConvection.h"

#include "Index.h"
#include "IndexLBM.h"
#include "defineLBM.h"


#ifdef USE_NVCC
namespace  LBMCalculationGPU {


__global__
void  to_euler_variables_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets,
    const real*  f_lbm,
          real*  rho,
          real*  u,
          real*  v,
          real*  w
    )
{
    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];
    const int  offset = mesh_offsets[idl];

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        const int  id          = Index   ::id0      (i,j,k, offset);
        const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);


        FuncLBMKernel::to_euler_variables(
            id,
            id_lbm_base,
            f_lbm,
            rho,
            u, v, w
            );
    }
    }
    }
}


__global__
void  stream_collision_sgs_gpu_wo_wall(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
          real*  sgs_vis,
    const real*  T,
    const real*  Tn,
    const real*  u,
    const real*  v,
    const real*  w,
    const real   kvis,
    const real   c_ref,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        int  ids[27];
        int  ids_lbm_base[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids         [idv_leaf] = Index   ::id      (i+iv,j+jv,k+kv, offset3d);
                    ids_lbm_base[idv_leaf] = IndexLBM::id_base0(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }


        FuncLBMKernel::stream_collision_cumulant_sgs_wo_wall(
                ids, ids_lbm_base,
                f_lbm, fn_lbm,
                sgs_vis,
                T, Tn,
                u, v, w,
                kvis,
                c_ref, dt
                );
    }
    }
    }
}


__global__
void  stream_collision_sgs_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
          real*  sgs_vis,
    const real*  T,
    const real*  Tn,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   kvis,
    const real   c_ref,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;


        int  ids[27];
        int  ids_lbm_base[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids         [idv_leaf] = Index   ::id      (i+iv,j+jv,k+kv, offset3d);
                    ids_lbm_base[idv_leaf] = IndexLBM::id_base0(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }

        FuncLBMKernel::stream_collision_cumulant_sgs(
                ids,
                ids_lbm_base,
                f_lbm, fn_lbm,
                sgs_vis,
                T, Tn,
                u,v,w,
                lv_obs,
                u_obs, v_obs, w_obs,
                kvis,
                c_ref, dt
                );
    }
    }
    }
}


__global__
void  scalar_advection_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real   dx,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        int  ids[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids[idv_leaf] = Index::id(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }

        FuncScalarAdvection::scalar_advection(
            ids,
            scalar, scalar_n,
            u,v,w,
            sgs_vis,
            lv_obs,
            dx, dt
            );
    }
    }
    }
}


__global__
void  scalar_advection_w_source_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  scalar,
          real*  scalar_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real*  sc_scalar,
    const real   dx,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        int  ids[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids[idv_leaf] = Index::id(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }

        FuncScalarAdvection::scalar_advection(
            ids,
            scalar, scalar_n,
            u,v,w,
            sgs_vis,
            lv_obs,
            sc_scalar,
            dx, dt
            );
    }
    }
    }
}


__global__
void  temperature_advection_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  T,
          real*  T_n,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  sgs_vis,
    const real*  lv_obs,
    const real*  T_obs,
    const real*  Tn_obs,
    const int*   bcTypes_T,
    const real   time_weight0,
    const real   coef_heatf,
    const real   dx,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        int  ids[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids[idv_leaf] = Index::id(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }

        FuncThermalConvection::cal_thermal_convection(
            ids,
            T, T_n,
            u,v,w,
            sgs_vis,
            lv_obs,
            T_obs, Tn_obs, bcTypes_T, time_weight0,
            coef_heatf,
            dx, dt
            );
    }
    }
    }
}


__global__
void  LBM_L2F_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  f_lbmF,
    const real*  f_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d_lbm [27];
    int  offset3dF_lbm[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d_lbm [idv] = mesh_offsets3x3x3[ idl          *nQ + idv] * nQ;
        offset3dF_lbm[idv] = mesh_offsets3x3x3[ id_child[idl]*nQ + idv] * nQ;
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int offset3d_lbm [27];
    __shared__ int offset3dF_lbm[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl          *nQ + _idx] * nQ;
            offset3dF_lbm[_idx] = mesh_offsets3x3x3[ id_child[idl]*nQ + _idx] * nQ;
        }
    }
    __syncthreads();
#endif


    const int  ii = (int(blockIdx.y%2));
    const int  jj = (int(blockIdx.y/2))%2;
    const int  kk = (int(blockIdx.y/4));

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        FuncAMRInterpolationCell::LBM_L2F_kernel(
                i,j,k, ii,jj,kk,
                f_lbmF,
                f_lbm,
                offset3dF_lbm, offset3d_lbm,
                lv,
                kvis, c_ref, dt0
                );
    }
    }
    }
}


__global__
void  LBM_L2FA_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  f_lbmF,
    const real*  f_lbm,
    const real*  fn_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d_lbm[27];
    int  offset3dF_lbm[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d_lbm [idv] = mesh_offsets3x3x3[ idl          *nQ + idv] * nQ;
        offset3dF_lbm[idv] = mesh_offsets3x3x3[ id_child[idl]*nQ + idv] * nQ;
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d_lbm[27];
    __shared__ int  offset3dF_lbm[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl          *nQ + _idx] * nQ;
            offset3dF_lbm[_idx] = mesh_offsets3x3x3[ id_child[idl]*nQ + _idx] * nQ;
        }
    }
    __syncthreads();
#endif


    const int  ii = (int(blockIdx.y%2));
    const int  jj = (int(blockIdx.y/2))%2;
    const int  kk = (int(blockIdx.y/4));

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        FuncAMRInterpolationCell::LBM_L2FA_kernel(
                i,j,k, ii,jj,kk,
                f_lbmF,
                f_lbm, fn_lbm,
                offset3dF_lbm, offset3d_lbm,
                lv,
                kvis, c_ref, dt0
                );
    }
    }
    }
}


__global__
void  Val_L2F_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  scalarF,
    const real*  scalar
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    int  offset3dF[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d [idv] = mesh_offsets3x3x3[ idl          *nQ + idv];
        offset3dF[idv] = mesh_offsets3x3x3[ id_child[idl]*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d[27];
    __shared__ int  offset3dF[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d [_idx] = mesh_offsets3x3x3[ idl          *nQ + _idx];
            offset3dF[_idx] = mesh_offsets3x3x3[ id_child[idl]*nQ + _idx];
        }
    }
    __syncthreads();
#endif

    const int  ii = (int(blockIdx.y%2));
    const int  jj = (int(blockIdx.y/2))%2;
    const int  kk = (int(blockIdx.y/4));

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;


        FuncAMRInterpolationCell::L2F_kernel(
                i,j,k, ii,jj,kk,
                scalarF,
                scalar,
                offset3dF, offset3d,
                lv
                );
    }
    }
    }
}


__global__
void  Val_L2FA_gpu(
    const int*   id_tasks,
    const int*   id_child,
    const int*   mesh_offsets3x3x3,
    const int    lv,
          real*  scalarF,
    const real*  scalar,
    const real*  scalarn
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    int  offset3dF[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d [idv] = mesh_offsets3x3x3[ idl          *nQ + idv];
        offset3dF[idv] = mesh_offsets3x3x3[ id_child[idl]*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d[27];
    __shared__ int  offset3dF[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d [_idx] = mesh_offsets3x3x3[ idl          *nQ + _idx];
            offset3dF[_idx] = mesh_offsets3x3x3[ id_child[idl]*nQ + _idx];
        }
    }
    __syncthreads();
#endif

    const int  ii = (int(blockIdx.y%2));
    const int  jj = (int(blockIdx.y/2))%2;
    const int  kk = (int(blockIdx.y/4));

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;


        FuncAMRInterpolationCell::L2FA_kernel(
                i,j,k, ii,jj,kk,
                scalarF,
                scalar,
                scalarn,
                offset3dF, offset3d,
                lv
                );
    }
    }
    }
}


__global__
void  LBM_F2L_gpu(
    const int*   id_tasks,
    const int*   id_parent,
    const int*   mesh_offsets3x3x3,
    const int    lvF,
    const int*   node_pos_x,
    const int*   node_pos_y,
    const int*   node_pos_z,
          real*  f_lbmC,
    const real*  f_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

    const int  idv0 = Index::idv0();
//    const int  offset_lbm  = mesh_offsets3x3x3[ idl           *nQ + idv0 ] * nQ;
    const int  offsetC_lbm = mesh_offsets3x3x3[ id_parent[idl]*nQ + idv0 ] * nQ;
    const int  lv_offset[3] = { node_pos_x[idl]%2, node_pos_y[idl]%2, node_pos_z[idl]%2 }; // [3] : 0 or 1 (left or right side in a coarse mesh) //

//    int  offset3d_lbm[27];
//    for (int idv=0; idv<nQ; idv++) {
//        offset3d_lbm [idv] = mesh_offsets3x3x3[ idl*nQ + idv] * nQ;
//    }
#ifndef USE_SHARED_MEMORY__
    int  offset3d_lbm[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d_lbm [idv] = mesh_offsets3x3x3[ idl*nQ + idv] * nQ;
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d_lbm[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl*nQ + _idx] * nQ;
        }
    }
    __syncthreads();
#endif

    const int  mx = int(DefAMR::NX_LEAF/2) / blockDim.x;
    const int  my = int(DefAMR::NX_LEAF/2) / blockDim.y;
    const int  mz = int(DefAMR::NX_LEAF/2) / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        FuncAMRInterpolationCell::LBM_L2C_kernel(
                i,j,k,
                f_lbmC, f_lbm,
                offsetC_lbm,
                offset3d_lbm,
                lv_offset,
                lvF,
                kvis, c_ref, dt0
                );
    }
    }
    }
}


__global__
void  LBM_F2LA_gpu(
    const int*   id_tasks,
    const int*   id_parent,
    const int*   mesh_offsets3x3x3,
    const int    lvF,
    const int*   node_pos_x,
    const int*   node_pos_y,
    const int*   node_pos_z,
          real*  f_lbmC,
          real*  fn_lbmC,
    const real*  f_lbm,
    const real   kvis,
    const real   c_ref,
    const real   dt0
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

    const int  idv0 = Index::idv0();
//    const int  offset_lbm  = mesh_offsets3x3x3[ idl           *nQ + idv0 ] * nQ;
    const int  offsetC_lbm = mesh_offsets3x3x3[ id_parent[idl]*nQ + idv0 ] * nQ;
    const int  lv_offset[3] = { node_pos_x[idl]%2, node_pos_y[idl]%2, node_pos_z[idl]%2 }; // [3] : 0 or 1 (left or right side in a coarse mesh) //

//    int  offset3d_lbm[27];
//    for (int idv=0; idv<nQ; idv++) {
//        offset3d_lbm [idv] = mesh_offsets3x3x3[ idl*nQ + idv] * nQ;
//    }
#ifndef USE_SHARED_MEMORY__
    int  offset3d_lbm[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d_lbm [idv] = mesh_offsets3x3x3[ idl*nQ + idv] * nQ;
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d_lbm[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl*nQ + _idx] * nQ;
        }
    }
    __syncthreads();
#endif

    const int  mx = int(DefAMR::NX_LEAF/2) / blockDim.x;
    const int  my = int(DefAMR::NX_LEAF/2) / blockDim.y;
    const int  mz = int(DefAMR::NX_LEAF/2) / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        FuncAMRInterpolationCell::LBM_L2CA_kernel(
                i,j,k,
                f_lbmC, fn_lbmC, f_lbm,
                offsetC_lbm,
                offset3d_lbm,
                lv_offset,
                lvF,
                kvis, c_ref, dt0
                );
    }
    }
    }
}


__global__
void  Val_F2L_gpu(
    const int*   id_tasks,
    const int*   id_parent,
    const int*   mesh_offsets3x3x3,
    const int    lvF,
    const int*   node_pos_x,
    const int*   node_pos_y,
    const int*   node_pos_z,
          real*  scalarC,
    const real*  scalar
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

    const int  idv0 = Index::idv0();
//    const int  offset  = mesh_offsets3x3x3[ idl           *nQ + idv0 ];
    const int  offsetC = mesh_offsets3x3x3[ id_parent[idl]*nQ + idv0 ];
    const int  lv_offset[3] = { node_pos_x[idl]%2, node_pos_y[idl]%2, node_pos_z[idl]%2 }; // [3] : 0 or 1 (left or right side in a coarse mesh) //

//    int  offset3d[27];
//    for (int idv=0; idv<nQ; idv++) {
//        offset3d [idv] = mesh_offsets3x3x3[ idl*nQ + idv];
//    }
#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d [idv] = mesh_offsets3x3x3[ idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) {
            offset3d [_idx] = mesh_offsets3x3x3[ idl*nQ + _idx];
        }
    }
    __syncthreads();
#endif

    const int  mx = int(DefAMR::NX_LEAF/2) / blockDim.x;
    const int  my = int(DefAMR::NX_LEAF/2) / blockDim.y;
    const int  mz = int(DefAMR::NX_LEAF/2) / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        FuncAMRInterpolationCell::L2C_kernel(
                i,j,k,
                scalarC, scalar,
                offsetC,
                offset3d,
                lv_offset,
                lvF
                );
    }
    }
    }
}


__global__
void  boundary_conditions_in_wall_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets,
    const real*  f_lbm,
          real*  fn_lbm,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
    const real*  lv_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real   c_ref,
    const real   dx,
    const real   dt
    )
{
    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];
	const int  offset = mesh_offsets[idl];

    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

	    const int  id          = Index   ::id0      (i,j,k, offset);
	    const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);

        FuncLBMKernel::flbm_in_wall(
            id,
            id_lbm_base,
            f_lbm, fn_lbm,
            lv_obs,
            u_obs, v_obs, w_obs
            );
    }
    }
    }
}


__global__
void  boundary_conditions_inflow_outflow_gpu(
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
    // solid wall //
    const real*  lv_obs,
    const real*  rho_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight,
    const real*  viscosity_weight,
    const real   c_ref,
    const real   dx,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d    [_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif


    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        int  ids[27], ids_lbm_base[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids         [idv_leaf] = Index   ::id      (i+iv,j+jv,k+kv, offset3d);
                    ids_lbm_base[idv_leaf] = IndexLBM::id_base0(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }

        const int  id          = ids         [ Index::idv0() ];
        const int  id_lbm_base = ids_lbm_base[ Index::idv0() ];

        FuncLBMKernel::flbm_inflow_outflow(
            ids,
            ids_lbm_base,
            f_lbm, fn_lbm,
            rho,u,v,w,
            lv_obs,
            rho_obs,
            u_obs, v_obs, w_obs,
            bcTypes_f,
            dirichlet_weight, viscosity_weight
            );
    }
    }
    }
}


__global__
void  boundary_condition_data_assimilation_gpu(
    const real   time_weight0,
    const int*   id_tasks,
    const int*   mesh_offsets3x3x3,
    const real*  f_lbm,
          real*  fn_lbm,
    const real*  rho,
    const real*  u,
    const real*  v,
    const real*  w,
          real*  Tn,
          real*  scalarn,
    // solid wall //
    const real*  lv_obs,
    const real*  rho_obs,
    const real*  u_obs,
    const real*  v_obs,
    const real*  w_obs,
    const real*  T_obs,
    const real*  scalar_obs,
    const real*  rhon_obs,
    const real*  un_obs,
    const real*  vn_obs,
    const real*  wn_obs,
    const real*  Tn_obs,
    const real*  scalarn_obs,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const real*  dirichlet_weight,
    const real*  viscosity_weight,
    const real   c_ref,
    const real   dx,
    const real   dt
    )
{
    constexpr int  nQ = LBM_velocity_model::nQ;

    const int  l   = blockIdx.x;
    const int  idl = id_tasks[l];

#ifndef USE_SHARED_MEMORY__
    int  offset3d[27];
    for (int idv=0; idv<nQ; idv++) {
        offset3d    [idv] = mesh_offsets3x3x3[idl*nQ + idv];
    }
#else // shared_memory //
    const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
    const int _tloop = int(27/_bDim) + 1;
    __shared__ int  offset3d[27];
    for (int _t=0; _t<_tloop; _t++) {
        const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
        if (_idx < 27) { offset3d    [_idx] = mesh_offsets3x3x3[idl*nQ + _idx]; }
    }
    __syncthreads();
#endif


    const int  mx = DefAMR::NX_LEAF / blockDim.x;
    const int  my = DefAMR::NX_LEAF / blockDim.y;
    const int  mz = DefAMR::NX_LEAF / blockDim.z;
    for (int _k=0; _k<mz; _k++) {
    for (int _j=0; _j<my; _j++) {
    for (int _i=0; _i<mx; _i++) {
        const int i = threadIdx.x + blockDim.x*_i;
        const int j = threadIdx.y + blockDim.y*_j;
        const int k = threadIdx.z + blockDim.z*_k;

        int  ids[27], ids_lbm_base[27];
        for (int kv=-1; kv<=1; kv++) {
            for (int jv=-1; jv<=1; jv++) {
                for (int iv=-1; iv<=1; iv++) {
                    const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                    ids         [idv_leaf] = Index   ::id      (i+iv,j+jv,k+kv, offset3d);
                    ids_lbm_base[idv_leaf] = IndexLBM::id_base0(i+iv,j+jv,k+kv, offset3d);
                }
            }
        }

        const int  id          = ids         [ Index::idv0() ];
        const int  id_lbm_base = ids_lbm_base[ Index::idv0() ];

        FuncLBMKernel::flbm_data_assimilation2(
            time_weight0,
            ids,
            ids_lbm_base,
            f_lbm, fn_lbm,
            Tn, scalarn,
            lv_obs,
            rho_obs,  u_obs,  v_obs,  w_obs,  T_obs,  scalar_obs,
            rhon_obs, un_obs, vn_obs, wn_obs, Tn_obs, scalarn_obs,
            bcTypes_f,
            dirichlet_weight, viscosity_weight
            );
    }
    }
    }
}


};
#endif
