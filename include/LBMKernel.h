#pragma once
#ifndef LBMKERNEL_H_
#define LBMKERNEL_H_

#include "definePrecision.h"
#include "defineAMR.h"

#include "foreach.h"
#include "FuncLBMKernel.h"

namespace LBMKernel {

template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  to_euler_variables(
    const int* id_tasks,
    const int  num_tasks,
    const int* mesh_offsets,
    const T0*  f_lbm,
          T1*  rho,
          T1*  u,
          T1*  v,
          T1*  w
    ) noexcept
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];
        const int  offset = mesh_offsets[idl];
    
        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
            const int  id          = Index   ::id0      (i,j,k, offset);
            const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);

            FuncLBMKernel::to_euler_variables(
                id,
                id_lbm_base,
                f_lbm,
                rho,
                u, v, w
                );
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1, typename T2>
inline __HD__
void  to_euler_variables_with_time_average(
    const int* id_tasks,
    const int  num_tasks,
    const int* mesh_offsets,
    const T0*  f_lbm,
          T1*  rho,
          T1*  u,
          T1*  v,
          T1*  w,
    const T1*  T,
          T1*  u_mean_1min,
          T1*  v_mean_1min,
          T1*  w_mean_1min,
          T1*  T_mean_1min,
          T1*  vel2_fluc_1min,
//          T1*  uu_fluc_1min,
//          T1*  vv_fluc_1min,
//          T1*  ww_fluc_1min,
    const T2   time_1min,
    const T2   dt
    ) noexcept
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];
        const int  offset = mesh_offsets[idl];
    
        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
            const int  id          = Index   ::id0      (i,j,k, offset);
            const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);

            FuncLBMKernel::to_euler_variables_with_time_average(
                id,
                id_lbm_base,
                f_lbm,
                rho,
                u, v, w, T,
                u_mean_1min,  v_mean_1min,  w_mean_1min,  T_mean_1min,
                vel2_fluc_1min,
//                uu_fluc_1min, vv_fluc_1min, ww_fluc_1min,
                time_1min,
                dt
                );
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  stream_collision_sgs(
    const int* id_tasks,
    const int  num_tasks,
    const int* mesh_offsets3x3x3,
    const T0*  f_lbm,
          T0*  fn_lbm,
          T0*  sgs_vis,
    const T0*  T,
    const T0*  Tn,
    const T0*  u,
    const T0*  v,
    const T0*  w,
    const T0*  lv_obj,
    const T0*  u_obj,
    const T0*  v_obj,
    const T0*  w_obj,
    const T0*  pad_obj,
    const T1   kvis,
    const T1   c_ref,
    const T1   dt
    ) noexcept
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

#ifndef __CUDA_ARCH__
        int  offset3d[27];
        for (int idv=0; idv<27; idv++) { offset3d[idv] = mesh_offsets3x3x3[idl*27 + idv]; }
#else // shared_memory //
        __shared__ int offset3d[27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*27 + _idx]; }
        }  __syncthreads();
#endif
        
    
        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
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
                    ids, ids_lbm_base,
                    f_lbm, fn_lbm,
                    sgs_vis,
                    T, Tn,
                    u,v,w,
                    lv_obj,
                    u_obj, v_obj, w_obj,
                    pad_obj,
                    kvis,
                    c_ref, dt
                    );
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  stream_collision_sgs_wo_wall(
    const int* id_tasks,
    const int  num_tasks,
    const int* mesh_offsets3x3x3,
    const T0*  f_lbm,
          T0*  fn_lbm,
          T0*  sgs_vis,
    const T0*  T,
    const T0*  Tn,
    const T0*  u,
    const T0*  v,
    const T0*  w,
    const T0*  pad_obj,
    const T1   kvis,
    const T1   c_ref,
    const T1   dt
    ) noexcept
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

#ifndef __CUDA_ARCH__
        int  offset3d[27];
        for (int idv=0; idv<27; idv++) { offset3d[idv] = mesh_offsets3x3x3[idl*27 + idv]; }
#else // shared_memory //
        __shared__ int offset3d[27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) { offset3d[_idx] = mesh_offsets3x3x3[idl*27 + _idx]; }
        }  __syncthreads();
#endif
        
    
        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
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
                    ids,
                    ids_lbm_base,
                    f_lbm, fn_lbm,
                    sgs_vis,
                    T, Tn,
                    u,v,w,
                    pad_obj,
                    kvis,
                    c_ref, dt
                    );
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1, typename T2>
inline __HD__
void  boundary_conditions_in_wall(
    const int* id_tasks,
    const int  num_tasks,
    const int* mesh_offsets,
          T0*  fn_lbm,
    const T1*  lv_obj,
    const T1*  u_obj,
    const T1*  v_obj,
    const T1*  w_obj,
    const T2   c_ref,
    const T2   dx,
    const T2   dt
    ) noexcept
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];
        const int  offset = mesh_offsets[idl];
    
        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
            const int  id          = Index   ::id0      (i,j,k, offset);
            const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);

            FuncLBMKernel::flbm_in_wall(
                id,
                id_lbm_base,
                fn_lbm,
                lv_obj,
                u_obj, v_obj, w_obj
                );
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1, typename T2>
inline __HD__
void  boundary_conditions_data_assimilation(
    const int    nn_max,
    const int    n_scalars,
    const real   time_weight0,
    const int*   id_tasks,
    const int    num_tasks,
    const int*   mesh_offsets,
    const T0*  f_lbm,
          T0*  fn_lbm,
    const T1*  u,
    const T1*  v,
    const T1*  w,
          T1*  Tn,
          T1*  scalarn,
    // solid wall //
    const T1*  lv_obj,
    const T1*  rho_obj,
    const T1*  u_obj,
    const T1*  v_obj,
    const T1*  w_obj,
    const T1*  T_obj,
    const T1*  scalar_obj,
    const T1*  rhon_obj,
    const T1*  un_obj,
    const T1*  vn_obj,
    const T1*  wn_obj,
    const T1*  Tn_obj,
    const T1*  scalarn_obj,
    // inflow & outflow bc //
    const int*   bcTypes_f,
    const T1*  dirichlet_weight,
    const T2   c_ref,
    const T2   dx,
    const T2   dt
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];
        const int  offset = mesh_offsets[idl];
    
        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
            const int  id          = Index   ::id0      (i,j,k, offset);
            const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);

#if 1
            FuncLBMKernel::flbm_data_assimilation2(
                time_weight0,
                n_scalars,
                nn_max,
                id,
                id_lbm_base,
                f_lbm, fn_lbm,
                Tn, scalarn,
                lv_obj,
                rho_obj,  u_obj,  v_obj,  w_obj,  T_obj,  scalar_obj,
                rhon_obj, un_obj, vn_obj, wn_obj, Tn_obj, scalarn_obj,
                bcTypes_f,
                dirichlet_weight
                );
#else
            FuncLBMKernel::flbm_data_assimilation2_force(
                time_weight0,
                n_scalars,
                nn_max,
                id,
                id_lbm_base,
                f_lbm, fn_lbm,
                u, v, w,
                Tn, scalarn,
                lv_obj,
                rho_obj,  u_obj,  v_obj,  w_obj,  T_obj,  scalar_obj,
                rhon_obj, un_obj, vn_obj, wn_obj, Tn_obj, scalarn_obj,
                bcTypes_f,
                dirichlet_weight
                );
#endif
        }
    }
}


}

#endif
