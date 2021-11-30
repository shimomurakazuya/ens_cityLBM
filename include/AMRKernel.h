#pragma once
#ifndef AMRKERNEL_H_
#define AMRKERNEL_H_


#include "definePrecision.h"
#include "defineAMR.h"

#include "foreach.h"
#include "FuncAMRInterpolationCell.h"


namespace  AMRKernel {


template<int NX_LEAF, typename T0>
inline __HD__
void  Val_L2F(
    const int* id_tasks,
    const int  num_tasks,
    const int* id_child,
    const int* mesh_offsets3x3x3,
    const int  lv,
          T0*  scalarF,
    const T0*  scalar
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

#ifndef __CUDA_ARCH__
        int  offset3d[27];
        int  offset3dF[27];
        for (int idv=0; idv<27; idv++) {
            offset3d [idv] = mesh_offsets3x3x3[ idl          *27 + idv];
            offset3dF[idv] = mesh_offsets3x3x3[ id_child[idl]*27 + idv];
        }
#else // shared_memory //
        __shared__ int  offset3d [27];
        __shared__ int  offset3dF[27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) {
                offset3d [_idx] = mesh_offsets3x3x3[ idl          *27 + _idx];
                offset3dF[_idx] = mesh_offsets3x3x3[ id_child[idl]*27 + _idx];
            }
        }  __syncthreads();
#endif

        FOR_EACH1D_BLOCKIDY(_ijk8, 8) {
            const int  ii = (int(_ijk8%2));
            const int  jj = (int(_ijk8/2))%2;
            const int  kk = (int(_ijk8/4));

            FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
                FuncAMRInterpolationCell::L2F_kernel(
                    i,j,k, ii,jj,kk,
                    scalarF,
                    scalar,
                    offset3dF, offset3d,
                    lv
                    );

            } // FOR_EACH3D
        } // FOR_EACH1D_BLOCKIDY
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0>
inline __HD__
void  Val_L2FA(
    const int* id_tasks,
    const int  num_tasks,
    const int* id_child,
    const int* mesh_offsets3x3x3,
    const int  lv,
          T0*  scalarF,
    const T0*  scalar,
    const T0*  scalarn
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

#ifndef __CUDA_ARCH__
        int  offset3d[27];
        int  offset3dF[27];
        for (int idv=0; idv<27; idv++) {
            offset3d [idv] = mesh_offsets3x3x3[ idl          *27 + idv];
            offset3dF[idv] = mesh_offsets3x3x3[ id_child[idl]*27 + idv];
        }
#else // shared_memory //
        __shared__ int  offset3d [27];
        __shared__ int  offset3dF[27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) {
                offset3d [_idx] = mesh_offsets3x3x3[ idl          *27 + _idx];
                offset3dF[_idx] = mesh_offsets3x3x3[ id_child[idl]*27 + _idx];
            }
        }  __syncthreads();
#endif

        FOR_EACH1D_BLOCKIDY(_ijk8, 8) {
            const int  ii = (int(_ijk8%2));
            const int  jj = (int(_ijk8/2))%2;
            const int  kk = (int(_ijk8/4));

            FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
                FuncAMRInterpolationCell::L2FA_kernel(
                    i,j,k, ii,jj,kk,
                    scalarF,
                    scalar, scalarn,
                    offset3dF, offset3d,
                    lv
                    );

            } // FOR_EACH3D
        } // FOR_EACH1D_BLOCKIDY
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0>
inline __HD__
void  Val_F2L(
    const int* id_tasks,
    const int  num_tasks,
    const int* id_parent,
    const int* mesh_offsets3x3x3,
    const int  lvF,
    const int* node_pos_x,
    const int* node_pos_y,
    const int* node_pos_z,
          T0*  scalarC,
    const T0*  scalar
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

        const int  idv0 = Index::idv0();
        const int  offsetC = mesh_offsets3x3x3[ id_parent[idl]*27 + idv0 ];
        const int  lv_offset[3] = { node_pos_x[idl]%2, node_pos_y[idl]%2, node_pos_z[idl]%2 }; // [3] : 0 or 1 (left or right side in a coarse mesh) //

#ifndef __CUDA_ARCH__
        int  offset3d[27];
        for (int idv=0; idv<27; idv++) {
            offset3d [idv] = mesh_offsets3x3x3[ idl*27 + idv];
        }
#else // shared_memory //
        __shared__ int  offset3d [27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) {
                offset3d [_idx] = mesh_offsets3x3x3[ idl*27 + _idx];
            }
        }  __syncthreads();
#endif

        FOR_EACH3D(i, j, k, (NX_LEAF/2),(NX_LEAF/2),(NX_LEAF/2)) {
            FuncAMRInterpolationCell::L2C_kernel(
                i,j,k,
                scalarC, scalar,
                offsetC,
                offset3d,
                lv_offset,
                lvF
                );
        
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  LBM_L2F(
    const int* id_tasks,
    const int  num_tasks,
    const int* id_child,
    const int* mesh_offsets3x3x3,
    const int    lv,
          T0*  f_lbmF,
    const T0*  f_lbm,
    const T1   kvis,
    const T1   c_ref,
    const T1   dt0
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

#ifndef __CUDA_ARCH__
        int  offset3d_lbm [27];
        int  offset3dF_lbm[27];
        for (int idv=0; idv<27; idv++) {
            offset3d_lbm [idv] = mesh_offsets3x3x3[ idl          *27 + idv] * 27;
            offset3dF_lbm[idv] = mesh_offsets3x3x3[ id_child[idl]*27 + idv] * 27;
        }
#else // shared_memory //
        __shared__ int offset3d_lbm [27];
        __shared__ int offset3dF_lbm[27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) {
                offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl          *27 + _idx] * 27;
                offset3dF_lbm[_idx] = mesh_offsets3x3x3[ id_child[idl]*27 + _idx] * 27;
            }
        }  __syncthreads();
#endif

        FOR_EACH1D_BLOCKIDY(_ijk8, 8) {
            const int  ii = (int(_ijk8%2));
            const int  jj = (int(_ijk8/2))%2;
            const int  kk = (int(_ijk8/4));

            FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
                FuncAMRInterpolationCell::LBM_L2F_kernel(
                    i,j,k, ii,jj,kk,
                    f_lbmF,
                    f_lbm,
                    offset3dF_lbm, offset3d_lbm,
                    lv,
                    kvis, c_ref, dt0
                    );

            } // FOR_EACH3D
        } // FOR_EACH1D_BLOCKIDY
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  LBM_L2FA(
    const int* id_tasks,
    const int  num_tasks,
    const int* id_child,
    const int* mesh_offsets3x3x3,
    const int    lv,
          T0*  f_lbmF,
    const T0*  f_lbm,
    const T0*  fn_lbm,
    const T1   kvis,
    const T1   c_ref,
    const T1   dt0
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

#ifndef __CUDA_ARCH__
        int  offset3d_lbm [27];
        int  offset3dF_lbm[27];
        for (int idv=0; idv<27; idv++) {
            offset3d_lbm [idv] = mesh_offsets3x3x3[ idl          *27 + idv] * 27;
            offset3dF_lbm[idv] = mesh_offsets3x3x3[ id_child[idl]*27 + idv] * 27;
        }
#else // shared_memory //
        __shared__ int offset3d_lbm [27];
        __shared__ int offset3dF_lbm[27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) {
                offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl          *27 + _idx] * 27;
                offset3dF_lbm[_idx] = mesh_offsets3x3x3[ id_child[idl]*27 + _idx] * 27;
            }
        }  __syncthreads();
#endif

        FOR_EACH1D_BLOCKIDY(_ijk8, 8) {
            const int  ii = (int(_ijk8%2));
            const int  jj = (int(_ijk8/2))%2;
            const int  kk = (int(_ijk8/4));

            FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
                FuncAMRInterpolationCell::LBM_L2FA_kernel(
                    i,j,k, ii,jj,kk,
                    f_lbmF,
                    f_lbm, fn_lbm,
                    offset3dF_lbm, offset3d_lbm,
                    lv,
                    kvis, c_ref, dt0
                    );

            } // FOR_EACH3D
        } // FOR_EACH1D_BLOCKIDY
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  LBM_F2L(
    const int* id_tasks,
    const int  num_tasks,
    const int* id_parent,
    const int* mesh_offsets3x3x3,
    const int  lvF,
    const int* node_pos_x,
    const int* node_pos_y,
    const int* node_pos_z,
          T0*  f_lbmC,
    const T0*  f_lbm,
    const T1   kvis,
    const T1   c_ref,
    const T1   dt0
    )
{
    FOR_EACH1D_BLOCKIDX(l, num_tasks) {
        const int  idl = id_tasks[l];

        const int  idv0 = Index::idv0();
        const int  offsetC_lbm = mesh_offsets3x3x3[ id_parent[idl]*27 + idv0 ] * 27;
        const int  lv_offset[3] = { node_pos_x[idl]%2, node_pos_y[idl]%2, node_pos_z[idl]%2 }; // [3] : 0 or 1 (left or right side in a coarse mesh) //
#ifndef __CUDA_ARCH__
        int  offset3d_lbm [27];
        for (int idv=0; idv<27; idv++) {
            offset3d_lbm [idv] = mesh_offsets3x3x3[ idl*27 + idv] * 27;
        }
#else // shared_memory //
        __shared__ int offset3d_lbm [27];
        const int _bDim  = blockDim.x*blockDim.y*blockDim.z;
        const int _tloop = int(27/_bDim) + 1;
        for (int _t=0; _t<_tloop; _t++) {
            const int _idx = threadIdx.x + blockDim.x*threadIdx.y + blockDim.x*blockDim.y*threadIdx.z + _bDim*_t;
            if (_idx < 27) {
                offset3d_lbm [_idx] = mesh_offsets3x3x3[ idl*27 + _idx] * 27;
            }
        }  __syncthreads();
#endif

        FOR_EACH3D(i, j, k, NX_LEAF,NX_LEAF,NX_LEAF) {
            FuncAMRInterpolationCell::LBM_L2C_kernel(
                i,j,k,
                f_lbmC, f_lbm,
                offsetC_lbm,
                offset3d_lbm,
                lv_offset,
                lvF,
                kvis, c_ref, dt0
                );
        
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


} // namespace  NSKernel


#endif
