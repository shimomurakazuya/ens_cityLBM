#pragma once
#ifndef NSKERNEL_H_
#define NSKERNEL_H_


#include "definePrecision.h"
#include "defineAMR.h"

#include "foreach.h"
#include "IndexMG.h"
#include "NSOperator.h"
#include "FuncDifferentialOperator.h"
#include "FuncHeatFluxModel.h"


namespace  NSKernel {


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  TemperatureAdvDifforg(
    const int*   id_tasks,
    const int    num_tasks,
    const int*   mesh_offsets3x3x3,
          T0*    Tn,
    const T0*    T, // operator //
    const T0*    u,
    const T0*    v,
    const T0*    w,
    const T0*    sgs_vis,
    const T0*    lv_obj,
    const T0*    T_obj,
    const T0*    Tn_obj,
    const int*   bcTypes_T,
    const T1     time_weight0,
    const T1     coef_heatf,
    const T1     dx,
    const T1     dt
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            auto val_weight0  = [&](auto a0, auto a1){ return a0*time_weight0 + a1*((T0)1.0 - time_weight0); };
        
            if ( FuncObj::is_obj( lv_obj[ID(0,0,0)] ) ) {
                if      ( bcTypes_T[ID(0,0,0)] == BCTypes::BCDirichlet ) { Tn[ID(0,0,0)] = T_obj[ID(0,0,0)]; }
                else if ( bcTypes_T[ID(0,0,0)] == BCTypes::BCRegion    ) { Tn[ID(0,0,0)] = NSOperator::average_neighbor<NX_LEAF>( offset3d, i,j,k, T ); }
                else                                                     { Tn[ID(0,0,0)] = NSOperator::average_w_obj   <NX_LEAF>( offset3d, i,j,k, T, lv_obj ); }
                SKIP_FOR();
            }


            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID( 0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID( 0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID( 0,  0,  0)], w[ID( 0,  0,  1)] };
        
            const T0  T3x[] = { T[ID(-1,  0,  0)], T[ID( 0,  0,  0)], T[ID( 1,  0,  0)] };
            const T0  T3y[] = { T[ID( 0, -1,  0)], T[ID( 0,  0,  0)], T[ID( 0,  1,  0)] };
            const T0  T3z[] = { T[ID( 0,  0, -1)], T[ID( 0,  0,  0)], T[ID( 0,  0,  1)] };

            const T0  lv_obj3x[] = { lv_obj[ID(-1,  0,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 1,  0,  0)] };
            const T0  lv_obj3y[] = { lv_obj[ID( 0, -1,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  1,  0)] };
            const T0  lv_obj3z[] = { lv_obj[ID( 0,  0, -1)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  0,  1)] };

            auto Tobj_weight0 = [&](int id)          { return val_weight0(T_obj[id], Tn_obj[id]); };
            const T0  Tobj3x[] = { Tobj_weight0(ID(-1,  0,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 1,  0,  0)) };
            const T0  Tobj3y[] = { Tobj_weight0(ID( 0, -1,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  1,  0)) };
            const T0  Tobj3z[] = { Tobj_weight0(ID( 0,  0, -1)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  0,  1)) };

            const int   bcTypesT3x[] = { bcTypes_T[ID(-1,  0,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 1,  0,  0)] };
            const int   bcTypesT3y[] = { bcTypes_T[ID( 0, -1,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  1,  0)] };
            const int   bcTypesT3z[] = { bcTypes_T[ID( 0,  0, -1)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  0,  1)] };
        
            // adv
            auto Tadv1st = [&](auto f3[], auto u3[], auto lv_obj3[]) {
                const T0  fxs[] = { (f3[1] - f3[0]),
                                    (f3[2] - f3[1]) };
            
                const T0  ufx_upwind = (u3[1] > (T0)0.0) ? u3[1]*fxs[0] : u3[1]*fxs[1];
            
                const T0  dft  = - ufx_upwind;
            
                return  dft;
            };

            auto Tadv2nd = [&](auto f3[], auto u3[], auto lv_obj3[]) {
                const T0  us[] = { ( FuncObj::is_obj( lv_obj3[0] ) || FuncObj::is_obj( lv_obj3[1] ) ) ? (T0)0.0 : (u3[0] + u3[1])*(T0)0.5,
                                   ( FuncObj::is_obj( lv_obj3[2] ) || FuncObj::is_obj( lv_obj3[1] ) ) ? (T0)0.0 : (u3[1] + u3[2])*(T0)0.5 };
            
                const T0  fxs[] = { f3[1] - f3[0],
                                    f3[2] - f3[1] };
            
                const T0  ufx_central = ( us[0]*fxs[0] + us[1]*fxs[1] ) * (T0)0.5;
            
                const T0  dft  = - ( ufx_central );
                const T0  dftt =  (u3[1]*u3[1])*(T0)0.5 * ( fxs[1] - fxs[0] );
            
            
                return  dft + (T0)0.5*dftt;
            };

            const T0  dT_adv1st = (   Tadv1st(T3x, u3, lv_obj3x)
                                    + Tadv1st(T3y, v3, lv_obj3y)
                                    + Tadv1st(T3z, w3, lv_obj3z) );

            const T0  dT_adv2nd = (   Tadv2nd(T3x, u3, lv_obj3x)
                                    + Tadv2nd(T3y, v3, lv_obj3y)
                                    + Tadv2nd(T3z, w3, lv_obj3z) );

            constexpr T0  ww   = (T0)0.02;
            const T0 dT_adv = ((T0)1.0-ww)*dT_adv1st + ww*dT_adv2nd;


            // diff
            auto Tdiff = [&](auto f3[], auto lv_obj3[], auto Tobj3[], auto bcTypesT3[], auto coef_heatf) {
                const T0  _T3[] = { (FuncObj::is_obj( lv_obj3[0] ) && bcTypesT3[0] == BCTypes::BCDirichlet) ? Tobj3[0] : f3[0],
                                    (FuncObj::is_obj( lv_obj3[1] ) && bcTypesT3[1] == BCTypes::BCDirichlet) ? Tobj3[1] : f3[1],
                                    (FuncObj::is_obj( lv_obj3[2] ) && bcTypesT3[2] == BCTypes::BCDirichlet) ? Tobj3[2] : f3[2]  };
            
                const T0  __T3[] = { (FuncObj::is_obj( lv_obj3[0] ) && bcTypesT3[0] == BCTypes::BCNeumann) ? _T3[1] : _T3[0],
                                     (FuncObj::is_obj( lv_obj3[1] ) && bcTypesT3[1] == BCTypes::BCNeumann) ? _T3[1] : _T3[1],
                                     (FuncObj::is_obj( lv_obj3[2] ) && bcTypesT3[2] == BCTypes::BCNeumann) ? _T3[1] : _T3[2]  };
            
                // dirichlet : wall is placed at the cell surface (dx -> 0.5*dx)
                const T0  _fx[] = { __T3[1] - __T3[0],
                                    __T3[2] - __T3[1] };
                const T0  fx[] = { FuncObj::is_obj( lv_obj3[0] ) ? _fx[0]*(T0)2.0 : _fx[0],
                                   FuncObj::is_obj( lv_obj3[2] ) ? _fx[1]*(T0)2.0 : _fx[1]  };
            
                const T0  dft = coef_heatf*(fx[1] - fx[0]);
            
                return  dft;
            };

            const T0 coef_heatf_sgs   = (sgs_vis[ID(0,0,0)]/air_property::Pr);
            const T0 _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_total = ( _coef_heatf_total < (T0)1.0/(T0)6.0/dt ) ? _coef_heatf_total : (T0)1.0/(T0)6.0/dt;

            const T0  dT_diff = dt * ( Tdiff(T3x, lv_obj3x, Tobj3x, bcTypesT3x, coef_heatf_total)
                                     + Tdiff(T3y, lv_obj3y, Tobj3y, bcTypesT3y, coef_heatf_total)
                                     + Tdiff(T3z, lv_obj3z, Tobj3z, bcTypesT3z, coef_heatf_total) );
        
            const T0 _Tn =  T[ID(0,0,0)] + dT_adv + dT_diff;

            Tn[ID(0,0,0)] =  _Tn;
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  TemperatureAdvDiff1st(
    const int*   id_tasks,
    const int    num_tasks,
    const int*   mesh_offsets3x3x3,
          T0*    Tn,
    const T0*    T, // operator //
    const T0*    u,
    const T0*    v,
    const T0*    w,
    const T0*    sgs_vis,
    const T0*    lv_obj,
    const T0*    T_obj,
    const T0*    Tn_obj,
    const int*   bcTypes_T,
    const T1     time_weight0,
    const T1     coef_heatf,
    const T1     dx,
    const T1     dt
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            auto val_weight0  = [&](auto a0, auto a1){ return a0*time_weight0 + a1*((T0)1.0 - time_weight0); };
        
            if ( FuncObj::is_obj( lv_obj[ID(0,0,0)] ) ) {
                if      ( bcTypes_T[ID(0,0,0)] == BCTypes::BCDirichlet ) { Tn[ID(0,0,0)] = T_obj[ID(0,0,0)]; }
                else if ( bcTypes_T[ID(0,0,0)] == BCTypes::BCRegion    ) { Tn[ID(0,0,0)] = NSOperator::average_neighbor<NX_LEAF>( offset3d, i,j,k, T ); }
                else                                                     { Tn[ID(0,0,0)] = NSOperator::average_w_obj   <NX_LEAF>( offset3d, i,j,k, T, lv_obj ); }
                SKIP_FOR();
            }


            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID( 0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID( 0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID( 0,  0,  0)], w[ID( 0,  0,  1)] };
        
            const T0  T3x[] = { T[ID(-1,  0,  0)], T[ID( 0,  0,  0)], T[ID( 1,  0,  0)] };
            const T0  T3y[] = { T[ID( 0, -1,  0)], T[ID( 0,  0,  0)], T[ID( 0,  1,  0)] };
            const T0  T3z[] = { T[ID( 0,  0, -1)], T[ID( 0,  0,  0)], T[ID( 0,  0,  1)] };

            const T0  lv_obj3x[] = { lv_obj[ID(-1,  0,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 1,  0,  0)] };
            const T0  lv_obj3y[] = { lv_obj[ID( 0, -1,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  1,  0)] };
            const T0  lv_obj3z[] = { lv_obj[ID( 0,  0, -1)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  0,  1)] };

            auto is_obj = [&](int id) { return FuncObj::is_obj(lv_obj[id]); };
            const bool  obj3x[] = { is_obj(ID(-1,  0,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 1,  0,  0)) };
            const bool  obj3y[] = { is_obj(ID( 0, -1,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  1,  0)) };
            const bool  obj3z[] = { is_obj(ID( 0,  0, -1)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  0,  1)) };

            auto Tobj_weight0 = [&](int id)          { return val_weight0(T_obj[id], Tn_obj[id]); };
            const T0  Tobj3x[] = { Tobj_weight0(ID(-1,  0,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 1,  0,  0)) };
            const T0  Tobj3y[] = { Tobj_weight0(ID( 0, -1,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  1,  0)) };
            const T0  Tobj3z[] = { Tobj_weight0(ID( 0,  0, -1)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  0,  1)) };

            const int   bcTypesT3x[] = { bcTypes_T[ID(-1,  0,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 1,  0,  0)] };
            const int   bcTypesT3y[] = { bcTypes_T[ID( 0, -1,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  1,  0)] };
            const int   bcTypesT3z[] = { bcTypes_T[ID( 0,  0, -1)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  0,  1)] };
        
            // adv
            auto adv_op = [&](auto... args) { return FuncDifferentialOperator::advection1d_fdm_euler_neumann_st3(args...); };
            const T0  dT_adv =  adv_op(T3x[0], T3x[1], T3x[2], u3[0], u3[1], u3[2], obj3x[0], obj3x[1], obj3x[2], 1.0, 1.0) // c*uc*dt/dx = 1*uc*1/1
                              + adv_op(T3y[0], T3y[1], T3y[2], v3[0], v3[1], v3[2], obj3y[0], obj3y[1], obj3y[2], 1.0, 1.0)
                              + adv_op(T3z[0], T3z[1], T3z[2], w3[0], w3[1], w3[2], obj3z[0], obj3z[1], obj3z[2], 1.0, 1.0);


            auto Txx_op = [&](T0& Txx_fl, T0& Txx_ob, const T0 T3[], const bool obj3[], const T0 Tobj3[]) {
                T0 Tx_fl[2], Tx_ob[2];
                if (obj3[0]) { Tx_fl[0] = 0.0;
                               Tx_ob[0] = T3[1] - Tobj3[0]; }
                else         { Tx_fl[0] = T3[1] - T3[0];
                               Tx_ob[0] = 0.0; }

                if (obj3[2]) { Tx_fl[1] = 0.0;
                               Tx_ob[1] = Tobj3[2] - T3[1]; }
                else         { Tx_fl[1] = T3[2] - T3[1];
                               Tx_ob[1] = 0.0; }

                Txx_fl = Tx_fl[1] - Tx_fl[0];
                Txx_ob = Tx_ob[1] - Tx_ob[0];
            };


            T0  Txx_fl, Tyy_fl, Tzz_fl; 
            T0  Txx_ob, Tyy_ob, Tzz_ob;
            Txx_op(Txx_fl, Txx_ob, T3x, obj3x, Tobj3x);
            Txx_op(Tyy_fl, Tyy_ob, T3y, obj3y, Tobj3y);
            Txx_op(Tzz_fl, Tzz_ob, T3z, obj3z, Tobj3z);

            const T0 coef_heatf_obj    = (coef_heatf) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_sgs    = (sgs_vis[ID(0,0,0)]/air_property::Pr);
            const T0 _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_total  = ( _coef_heatf_total < (T0)1.0/(T0)6.0/dt ) ? _coef_heatf_total : (T0)1.0/(T0)6.0/dt;

            const T0 dT_diff = (  (Txx_fl + Tyy_fl + Tzz_fl)*coef_heatf_total 
                                + (Txx_ob + Tyy_ob + Tzz_ob)*coef_heatf_obj   ) * dt;

            const T0 _Tn =  T[ID(0,0,0)] + dT_adv + dT_diff;

            Tn[ID(0,0,0)] =  _Tn;
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  TemperatureAdvDiffweno(
    const int*   id_tasks,
    const int    num_tasks,
    const int*   mesh_offsets3x3x3,
          T0*    Tn,
    const T0*    T, // operator //
    const T0*    u,
    const T0*    v,
    const T0*    w,
    const T0*    sgs_vis,
    const T0*    lv_obj,
    const T0*    T_obj,
    const T0*    Tn_obj,
    const int*   bcTypes_T,
    const T1     time_weight0,
    const T1     coef_heatf,
    const T1     dx,
    const T1     dt
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            auto val_weight0  = [&](auto a0, auto a1){ return a0*time_weight0 + a1*((T0)1.0 - time_weight0); };
        
        
            if ( FuncObj::is_obj( lv_obj[ID(0,0,0)] ) ) {
                if      ( bcTypes_T[ID(0,0,0)] == BCTypes::BCDirichlet ) { Tn[ID(0,0,0)] = T_obj[ID(0,0,0)]; }
                else if ( bcTypes_T[ID(0,0,0)] == BCTypes::BCRegion    ) { Tn[ID(0,0,0)] = NSOperator::average_neighbor<NX_LEAF>( offset3d, i,j,k, T ); }
                else                                                     { Tn[ID(0,0,0)] = NSOperator::average_w_obj   <NX_LEAF>( offset3d, i,j,k, T, lv_obj ); }
                SKIP_FOR();
            }


            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID( 0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID( 0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID( 0,  0,  0)], w[ID( 0,  0,  1)] };
        
            const T0  T3x[] = { T[ID(-1,  0,  0)], T[ID( 0,  0,  0)], T[ID( 1,  0,  0)] };
            const T0  T3y[] = { T[ID( 0, -1,  0)], T[ID( 0,  0,  0)], T[ID( 0,  1,  0)] };
            const T0  T3z[] = { T[ID( 0,  0, -1)], T[ID( 0,  0,  0)], T[ID( 0,  0,  1)] };

            const T0  T7x[] = { T[ID(-3,  0,  0)], T[ID(-2,  0,  0)], T[ID(-1,  0,  0)], T[ID( 0,  0,  0)], T[ID( 1,  0,  0)], T[ID( 2,  0,  0)], T[ID( 3,  0,  0)] };
            const T0  T7y[] = { T[ID( 0, -3,  0)], T[ID( 0, -2,  0)], T[ID( 0, -1,  0)], T[ID( 0,  0,  0)], T[ID( 0,  1,  0)], T[ID( 0,  2,  0)], T[ID( 0,  3,  0)] };
            const T0  T7z[] = { T[ID( 0,  0, -3)], T[ID( 0,  0, -2)], T[ID( 0,  0, -1)], T[ID( 0,  0,  0)], T[ID( 0,  0,  1)], T[ID( 0,  0,  2)], T[ID( 0,  0,  3)] };

            auto is_obj = [&](int id) { return FuncObj::is_obj(lv_obj[id]); };
            const bool  obj3x[] = { is_obj(ID(-1,  0,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 1,  0,  0)) };
            const bool  obj3y[] = { is_obj(ID( 0, -1,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  1,  0)) };
            const bool  obj3z[] = { is_obj(ID( 0,  0, -1)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  0,  1)) };

            const bool  obj7x[] = { is_obj(ID(-3,  0,  0)), is_obj(ID(-2,  0,  0)), is_obj(ID(-1,  0,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 1,  0,  0)), is_obj(ID( 2,  0,  0)), is_obj(ID( 3,  0,  0)) };
            const bool  obj7y[] = { is_obj(ID( 0, -3,  0)), is_obj(ID( 0, -2,  0)), is_obj(ID( 0, -1,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  1,  0)), is_obj(ID( 0,  2,  0)), is_obj(ID( 0,  3,  0)) };
            const bool  obj7z[] = { is_obj(ID( 0,  0, -3)), is_obj(ID( 0,  0, -2)), is_obj(ID( 0,  0, -1)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  0,  1)), is_obj(ID( 0,  0,  2)), is_obj(ID( 0,  0,  3)) };

            auto Tobj_weight0 = [&](int id)          { return val_weight0(T_obj[id], Tn_obj[id]); };
            const T0  Tobj3x[] = { Tobj_weight0(ID(-1,  0,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 1,  0,  0)) };
            const T0  Tobj3y[] = { Tobj_weight0(ID( 0, -1,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  1,  0)) };
            const T0  Tobj3z[] = { Tobj_weight0(ID( 0,  0, -1)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  0,  1)) };

            const int   bcTypesT3x[] = { bcTypes_T[ID(-1,  0,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 1,  0,  0)] };
            const int   bcTypesT3y[] = { bcTypes_T[ID( 0, -1,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  1,  0)] };
            const int   bcTypesT3z[] = { bcTypes_T[ID( 0,  0, -1)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  0,  1)] };
        
            // adv
            auto adv_op = [&](auto... args) { return FuncDifferentialOperator::advection1d_fdm_euler_neumann_st7(args...); };
            const T0  dT_adv =  adv_op(T7x[0], T7x[1], T7x[2], T7x[3], T7x[4], T7x[5], T7x[6], u3[1], obj7x[0], obj7x[1], obj7x[2], obj7x[3], obj7x[4], obj7x[5], obj7x[6], 1.0, 1.0) // c*uc*dt/dx = 1*uc*1/1
                              + adv_op(T7y[0], T7y[1], T7y[2], T7y[3], T7y[4], T7y[5], T7y[6], v3[1], obj7y[0], obj7y[1], obj7y[2], obj7y[3], obj7y[4], obj7y[5], obj7y[6], 1.0, 1.0)
                              + adv_op(T7z[0], T7z[1], T7z[2], T7z[3], T7z[4], T7z[5], T7z[6], w3[1], obj7z[0], obj7z[1], obj7z[2], obj7z[3], obj7z[4], obj7z[5], obj7z[6], 1.0, 1.0);


            // diff
//            auto Tdiff = [&](auto f3[], auto lv_obj3[], auto Tobj3[], auto bcTypesT3[], auto coef_heatf) {
//                const T0  _T3[] = { (FuncObj::is_obj( lv_obj3[0] ) && bcTypesT3[0] == BCTypes::BCDirichlet) ? Tobj3[0] : f3[0],
//                                    (FuncObj::is_obj( lv_obj3[1] ) && bcTypesT3[1] == BCTypes::BCDirichlet) ? Tobj3[1] : f3[1],
//                                    (FuncObj::is_obj( lv_obj3[2] ) && bcTypesT3[2] == BCTypes::BCDirichlet) ? Tobj3[2] : f3[2]  };
//            
//                const T0  __T3[] = { (FuncObj::is_obj( lv_obj3[0] ) && bcTypesT3[0] == BCTypes::BCNeumann) ? _T3[1] : _T3[0],
//                                     (FuncObj::is_obj( lv_obj3[1] ) && bcTypesT3[1] == BCTypes::BCNeumann) ? _T3[1] : _T3[1],
//                                     (FuncObj::is_obj( lv_obj3[2] ) && bcTypesT3[2] == BCTypes::BCNeumann) ? _T3[1] : _T3[2]  };
//            
//                // dirichlet : wall is placed at the cell surface (dx -> 0.5*dx)
//                const T0  _fx[] = { __T3[1] - __T3[0],
//                                    __T3[2] - __T3[1] };
//                const T0  fx[] = { FuncObj::is_obj( lv_obj3[0] ) ? _fx[0]*(T0)2.0 : _fx[0],
//                                   FuncObj::is_obj( lv_obj3[2] ) ? _fx[1]*(T0)2.0 : _fx[1]  };
//            
//                const T0  dft = coef_heatf*(fx[1] - fx[0]);
//            
//                return  dft;
//            };
//
//            const T0 coef_heatf_sgs   = (sgs_vis[ID(0,0,0)]/air_property::Pr);
//            const T0 _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
//            const T0 coef_heatf_total = ( _coef_heatf_total < (T0)1.0/(T0)6.0/dt ) ? _coef_heatf_total : (T0)1.0/(T0)6.0/dt;
//
//            const T0  dT_diff = dt * ( Tdiff(T3x, lv_obj3x, Tobj3x, bcTypesT3x, coef_heatf_total)
//                                     + Tdiff(T3y, lv_obj3y, Tobj3y, bcTypesT3y, coef_heatf_total)
//                                     + Tdiff(T3z, lv_obj3z, Tobj3z, bcTypesT3z, coef_heatf_total) );
        
            auto Txx_op = [&](auto& Txx_fl, auto& Txx_ob, auto T3[], auto obj3[], auto Tobj3[]) {
                T0 Tx_fl[2], Tx_ob[2];
                if (obj3[0]) { Tx_fl[0] = 0.0;
                               Tx_ob[0] = T3[1] - Tobj3[0]; }
                else         { Tx_fl[0] = T3[1] - T3[0];
                               Tx_ob[0] = 0.0; }

                if (obj3[2]) { Tx_fl[1] = 0.0;
                               Tx_ob[1] = Tobj3[2] - T3[1]; }
                else         { Tx_fl[1] = T3[2] - T3[1];
                               Tx_ob[1] = 0.0; }

                Txx_fl = Tx_fl[1] - Tx_fl[0];
                Txx_ob = Tx_ob[1] - Tx_ob[0];
            };


            T0  Txx_fl, Tyy_fl, Tzz_fl; 
            T0  Txx_ob, Tyy_ob, Tzz_ob;
            Txx_op(Txx_fl, Txx_ob, T3x, obj3x, Tobj3x);
            Txx_op(Tyy_fl, Tyy_ob, T3y, obj3y, Tobj3y);
            Txx_op(Tzz_fl, Tzz_ob, T3z, obj3z, Tobj3z);

            const T0 coef_heatf_obj    = (coef_heatf) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_sgs    = (sgs_vis[ID(0,0,0)]/air_property::Pr);
            const T0 _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_total  = ( _coef_heatf_total < (T0)1.0/(T0)6.0/dt ) ? _coef_heatf_total : (T0)1.0/(T0)6.0/dt;

            const T0 dT_diff = (  (Txx_fl + Tyy_fl + Tzz_fl)*coef_heatf_total 
                                + (Txx_ob + Tyy_ob + Tzz_ob)*coef_heatf_obj   ) * dt;

            const T0 _Tn =  T[ID(0,0,0)] + dT_adv + dT_diff;

            Tn[ID(0,0,0)] =  _Tn;
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  TemperatureAdvDiff1st_hflux(
    const int*   id_tasks,
    const int    num_tasks,
    const int*   mesh_offsets3x3x3,
          T0*    Tn,
    const T0*    T, // operator //
    const T0*    u,
    const T0*    v,
    const T0*    w,
    const T0*    sgs_vis,
    const T0*    lv_obj,
    const T0*    T_obj,
    const T0*    Tn_obj,
    const T0*    hflux_z_obj,
    const T0*    hfluxn_z_obj,
    const T1     time_weight0,
    const T1     coef_heatf,
    const T1     dx,
    const T1     dt
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            auto val_weight0  = [&](auto a0, auto a1){ return a0*time_weight0 + a1*((T0)1.0 - time_weight0); };
        
            if ( FuncObj::is_obj( lv_obj[ID(0,0,0)] ) ) {
                Tn[ID(0,0,0)] = NSOperator::average_w_obj   <NX_LEAF>( offset3d, i,j,k, T, lv_obj );
                SKIP_FOR();
            }


            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID( 0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID( 0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID( 0,  0,  0)], w[ID( 0,  0,  1)] };
        
            const T0  T3x[] = { T[ID(-1,  0,  0)], T[ID( 0,  0,  0)], T[ID( 1,  0,  0)] };
            const T0  T3y[] = { T[ID( 0, -1,  0)], T[ID( 0,  0,  0)], T[ID( 0,  1,  0)] };
            const T0  T3z[] = { T[ID( 0,  0, -1)], T[ID( 0,  0,  0)], T[ID( 0,  0,  1)] };

            const T0  lv_obj3x[] = { lv_obj[ID(-1,  0,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 1,  0,  0)] };
            const T0  lv_obj3y[] = { lv_obj[ID( 0, -1,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  1,  0)] };
            const T0  lv_obj3z[] = { lv_obj[ID( 0,  0, -1)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  0,  1)] };

            auto is_obj = [&](int id) { return FuncObj::is_obj(lv_obj[id]); };
            const bool  obj3x[] = { is_obj(ID(-1,  0,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 1,  0,  0)) };
            const bool  obj3y[] = { is_obj(ID( 0, -1,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  1,  0)) };
            const bool  obj3z[] = { is_obj(ID( 0,  0, -1)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  0,  1)) };

//            auto Tobj_weight0 = [&](int id)          { return val_weight0(T_obj[id], Tn_obj[id]); };
//            const T0  Tobj3x[] = { Tobj_weight0(ID(-1,  0,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 1,  0,  0)) };
//            const T0  Tobj3y[] = { Tobj_weight0(ID( 0, -1,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  1,  0)) };
//            const T0  Tobj3z[] = { Tobj_weight0(ID( 0,  0, -1)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  0,  1)) };
//
//            const int   bcTypesT3x[] = { bcTypes_T[ID(-1,  0,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 1,  0,  0)] };
//            const int   bcTypesT3y[] = { bcTypes_T[ID( 0, -1,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  1,  0)] };
//            const int   bcTypesT3z[] = { bcTypes_T[ID( 0,  0, -1)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  0,  1)] };
        
            auto hflux_weight0 = [&](int id)  { return val_weight0(hflux_z_obj[id], hfluxn_z_obj[id]); };
            const T0  hflux2z[] = { hflux_weight0(ID( 0,  0,  0)), hflux_weight0(ID( 0,  0,  1)) };

            // adv
            auto adv_op = [&](auto... args) { return FuncDifferentialOperator::advection1d_fdm_euler_neumann_st3(args...); };
            const T0  dT_adv =  adv_op(T3x[0], T3x[1], T3x[2], u3[0], u3[1], u3[2], obj3x[0], obj3x[1], obj3x[2], 1.0, 1.0) // c*uc*dt/dx = 1*uc*1/1
                              + adv_op(T3y[0], T3y[1], T3y[2], v3[0], v3[1], v3[2], obj3y[0], obj3y[1], obj3y[2], 1.0, 1.0)
                              + adv_op(T3z[0], T3z[1], T3z[2], w3[0], w3[1], w3[2], obj3z[0], obj3z[1], obj3z[2], 1.0, 1.0);


            auto Txx_op = [&](T0& Txx_fl, const T0 T3[], const bool obj3[]) {
                T0 Tx_fl[2];
                Tx_fl[0] = obj3[0] ? 0.0 : T3[1] - T3[0];
                Tx_fl[1] = obj3[2] ? 0.0 : T3[2] - T3[1];

                Txx_fl = Tx_fl[1] - Tx_fl[0];
            };


            T0  Txx_fl, Tyy_fl, Tzz_fl; 
            Txx_op(Txx_fl, T3x, obj3x);
            Txx_op(Tyy_fl, T3y, obj3y);
            Txx_op(Tzz_fl, T3z, obj3z);

            // hflux : direction 0 to 1
            T0 hflux_zm = obj3z[0] ? hflux2z[0] : 0.0; // K m s^-1
            T0 hflux_zp = obj3z[2] ? hflux2z[1] : 0.0;
            T0 hflux_z = ( hflux_zm - hflux_zp ) / dx; // K s^-1

            // avoid heat flux in narrow cell
            auto true2int = [](bool a) -> int { return (a == true) ? 1 : 0; };
            int true_obj = true2int(obj3x[0]) + true2int(obj3x[2]) + true2int(obj3y[0]) + true2int(obj3y[2]) + true2int(obj3z[0]) + true2int(obj3z[2]);
            if (true_obj >= 3) { hflux_z = 0.0; }

//            if ( (obj3x[0] && obj3x[2]) || (obj3y[0] && obj3y[2]) || (obj3z[0] && obj3z[2]) ) { hflux_z = 0.0; }

            const T0 coef_heatf_obj    = (coef_heatf) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_sgs    = (sgs_vis[ID(0,0,0)]/air_property::Pr);
            const T0 _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_total  = ( _coef_heatf_total < (T0)1.0/(T0)6.0/dt ) ? _coef_heatf_total : (T0)1.0/(T0)6.0/dt;

            const T0 dT_diff = (  (Txx_fl + Tyy_fl + Tzz_fl)*coef_heatf_total 
                                +  hflux_z  ) * dt;

            const T0 _Tn =  T[ID(0,0,0)] + dT_adv + dT_diff;

            Tn[ID(0,0,0)] =  _Tn;
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
inline __HD__
void  TemperatureAdvDiff1st_hflux_les(
    const int*   id_tasks,
    const int    num_tasks,
    const int*   mesh_offsets3x3x3,
          T0*    Tn,
    const T0*    T, // operator //
    const T0*    u,
    const T0*    v,
    const T0*    w,
    const T0*    sgs_vis,
    const T0*    lv_obj,
    const T0*    T_obj,
    const T0*    Tn_obj,
    const T1     time_weight0,
    const T1     coef_heatf,
    const T1     dx,
    const T1     dt,
    const T1     c_ref,
    const T0*    u_av,
    const T0*    v_av,
    const T0*    w_av,
    const T0*    T_av
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            auto val_weight0  = [&](auto a0, auto a1){ return a0*time_weight0 + a1*((T0)1.0 - time_weight0); };
        
            if ( FuncObj::is_obj( lv_obj[ID(0,0,0)] ) ) {
                Tn[ID(0,0,0)] = NSOperator::average_w_obj   <NX_LEAF>( offset3d, i,j,k, T, lv_obj );
                SKIP_FOR();
            }


            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID( 0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID( 0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID( 0,  0,  0)], w[ID( 0,  0,  1)] };
        
            const T0  T3x[] = { T[ID(-1,  0,  0)], T[ID( 0,  0,  0)], T[ID( 1,  0,  0)] };
            const T0  T3y[] = { T[ID( 0, -1,  0)], T[ID( 0,  0,  0)], T[ID( 0,  1,  0)] };
            const T0  T3z[] = { T[ID( 0,  0, -1)], T[ID( 0,  0,  0)], T[ID( 0,  0,  1)] };

            const T0  lv_obj3x[] = { lv_obj[ID(-1,  0,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 1,  0,  0)] };
            const T0  lv_obj3y[] = { lv_obj[ID( 0, -1,  0)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  1,  0)] };
            const T0  lv_obj3z[] = { lv_obj[ID( 0,  0, -1)], lv_obj[ID( 0,  0,  0)], lv_obj[ID( 0,  0,  1)] };

            auto is_obj = [&](int id) { return FuncObj::is_obj(lv_obj[id]); };
            const bool  obj3x[] = { is_obj(ID(-1,  0,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 1,  0,  0)) };
            const bool  obj3y[] = { is_obj(ID( 0, -1,  0)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  1,  0)) };
            const bool  obj3z[] = { is_obj(ID( 0,  0, -1)), is_obj(ID( 0,  0,  0)), is_obj(ID( 0,  0,  1)) };

            auto Tobj_weight0 = [&](int id)          { return val_weight0(T_obj[id], Tn_obj[id]); };
//            const T0  Tobj3x[] = { Tobj_weight0(ID(-1,  0,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 1,  0,  0)) };
//            const T0  Tobj3y[] = { Tobj_weight0(ID( 0, -1,  0)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  1,  0)) };
//            const T0  Tobj3z[] = { Tobj_weight0(ID( 0,  0, -1)), Tobj_weight0(ID( 0,  0,  0)), Tobj_weight0(ID( 0,  0,  1)) };
//
//            const int   bcTypesT3x[] = { bcTypes_T[ID(-1,  0,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 1,  0,  0)] };
//            const int   bcTypesT3y[] = { bcTypes_T[ID( 0, -1,  0)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  1,  0)] };
//            const int   bcTypesT3z[] = { bcTypes_T[ID( 0,  0, -1)], bcTypes_T[ID( 0,  0,  0)], bcTypes_T[ID( 0,  0,  1)] };
        
            // adv
            auto adv_op = [&](auto... args) { return FuncDifferentialOperator::advection1d_fdm_euler_neumann_st3(args...); };
            const T0  dT_adv =  adv_op(T3x[0], T3x[1], T3x[2], u3[0], u3[1], u3[2], obj3x[0], obj3x[1], obj3x[2], 1.0, 1.0) // c*uc*dt/dx = 1*uc*1/1
                              + adv_op(T3y[0], T3y[1], T3y[2], v3[0], v3[1], v3[2], obj3y[0], obj3y[1], obj3y[2], 1.0, 1.0)
                              + adv_op(T3z[0], T3z[1], T3z[2], w3[0], w3[1], w3[2], obj3z[0], obj3z[1], obj3z[2], 1.0, 1.0);


            auto Txx_op = [&](T0& Txx_fl, const T0 T3[], const bool obj3[]) {
                T0 Tx_fl[2];
                Tx_fl[0] = obj3[0] ? 0.0 : T3[1] - T3[0];
                Tx_fl[1] = obj3[2] ? 0.0 : T3[2] - T3[1];

                Txx_fl = Tx_fl[1] - Tx_fl[0];
            };


            T0  Txx_fl, Tyy_fl, Tzz_fl; 
            Txx_op(Txx_fl, T3x, obj3x);
            Txx_op(Tyy_fl, T3y, obj3y);
            Txx_op(Tzz_fl, T3z, obj3z);

            // hflux : direction 0 to 1

            auto hflux_zdir = [&](int _k){
                const T0 us_av = u_av[ID(0,0,_k)];
                const T0 vs_av = v_av[ID(0,0,_k)];
                const T0 Ts_av = T_av[ID(0,0,_k)];
//                const T0 us_av = NSOperator::xy_average_w_obj   <NX_LEAF>( offset3d, i,j,k+_k, u_av, lv_obj );
//                const T0 vs_av = NSOperator::xy_average_w_obj   <NX_LEAF>( offset3d, i,j,k+_k, v_av, lv_obj );
//                const T0 Ts_av = NSOperator::xy_average_w_obj   <NX_LEAF>( offset3d, i,j,k+_k, T_av, lv_obj );

                const T0 u_chi = sqrt( us_av*us_av + vs_av*vs_av ) * c_ref;
                const T0 T_chi = Ts_av;
                const T0 Tb    = Tobj_weight0(ID(0,0, _k-1));

                return FuncHeatFluxModel::q3w_HeatFlux(u_chi, T_chi, Tb, (T0)(0.5*dx), (T0)fluid_property::Temperature_bouyancy);
            };


//            T0 hflux_zm = obj3z[0] ? hflux_zdir(0, hflux2z[0]) : 0.0; // K m s^-1
//            T0 hflux_zp = obj3z[2] ? hflux_zdir(1, hflux2z[1]) : 0.0;
//            T0 hflux_z = ( hflux_zm - hflux_zp ) / dx; // K s^-1

            T0 hflux_zm = obj3z[0] ? hflux_zdir(0) : 0.0; // K m s^-1
            T0 hflux_z = ( hflux_zm ) / dx; // K s^-1

            // avoid heat flux in narrow cell
            auto true2int = [](bool a) -> int { return (a == true) ? 1 : 0; };
            int true_obj = true2int(obj3x[0]) + true2int(obj3x[2]) + true2int(obj3y[0]) + true2int(obj3y[2]) + true2int(obj3z[0]) + true2int(obj3z[2]);
            if (true_obj >= 3) { hflux_z = 0.0; }

//            if ( (obj3x[0] && obj3x[2]) || (obj3y[0] && obj3y[2]) || (obj3z[0] && obj3z[2]) ) { hflux_z = 0.0; }

            const T0 coef_heatf_obj    = (coef_heatf) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_sgs    = (sgs_vis[ID(0,0,0)]/air_property::Pr);
            const T0 _coef_heatf_total = (coef_heatf + coef_heatf_sgs) / (dx*dx); // because dx of the fxx set to 1 //
            const T0 coef_heatf_total  = ( _coef_heatf_total < (T0)1.0/(T0)6.0/dt ) ? _coef_heatf_total : (T0)1.0/(T0)6.0/dt;

            const T0 dT_diff = (  (Txx_fl + Tyy_fl + Tzz_fl)*coef_heatf_total 
                                +  hflux_z  ) * dt;

            const T0 _Tn =  T[ID(0,0,0)] + dT_adv + dT_diff;

            Tn[ID(0,0,0)] =  _Tn;
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
__HD__
void  ScalarConvectionorg(
    const int*  id_tasks,
    const int   num_tasks,
    const int*  mesh_offsets3x3x3,
    const int   nn_max,
    const int   n_scalars,
          T0*   scalar_n,
    const T0*   scalar,
    const T0*   u,
    const T0*   v,
    const T0*   w,
    const T0*   lv_obj,
    const T0*   sc_scalar,
    const bool  flag_source,
    const T1    dx,
    const T1    dt
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            if ( FuncObj::is_obj(lv_obj[ID(0,0,0)]) ) {
                for(int n=0; n<n_scalars; n++) {
                    int idx = n * nn_max;

                    scalar_n[idx + ID(0, 0, 0)] = NSOperator::average_neighbor_w_obj<NX_LEAF>(offset3d, i,j,k, &scalar[idx], lv_obj);
                }
                SKIP_FOR();
            }
        
            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID(0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID(0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID(0,  0,  0)], w[ID( 0,  0,  1)] };

            const T0  lv_obj3x[] = { lv_obj[ID(-1,  0,  0)], lv_obj[ID(0,  0,  0)], lv_obj[ID( 1,  0,  0)] };
            const T0  lv_obj3y[] = { lv_obj[ID( 0, -1,  0)], lv_obj[ID(0,  0,  0)], lv_obj[ID( 0,  1,  0)] };
            const T0  lv_obj3z[] = { lv_obj[ID( 0,  0, -1)], lv_obj[ID(0,  0,  0)], lv_obj[ID( 0,  0,  1)] };

            auto adv1st = [&](const T0 f3[], const T0 vel3[], const T0 lv_obj3[]){
                const T0  us[] = { ( FuncObj::is_obj( lv_obj3[0] ) ) ? (T0)0.0 : (vel3[0] + vel3[1])*(T0)0.5,
                                   ( FuncObj::is_obj( lv_obj3[2] ) ) ? (T0)0.0 : (vel3[1] + vel3[2])*(T0)0.5 };
            
                const T0  fs[] = { ( us[0] >= (T0)0.0 ) ? f3[0] : f3[1],
                                   ( us[1] >= (T0)0.0 ) ? f3[1] : f3[2] };
            
                const T0  dft  = - ( us[1]*fs[1] - us[0]*fs[0] );
            
                return  dft;
            };

            auto advection_f = [&](T0* fn, const T0* f, const T0* sc_f) {
                const T0  f3x[] = { f[ID(-1,  0,  0)], f[ID(0,  0,  0)], f[ID( 1,  0,  0)] };
                const T0  f3y[] = { f[ID( 0, -1,  0)], f[ID(0,  0,  0)], f[ID( 0,  1,  0)] };
                const T0  f3z[] = { f[ID( 0,  0, -1)], f[ID(0,  0,  0)], f[ID( 0,  0,  1)] };


                // convection //
                const T0  df_conv =  adv1st(f3x, u3, lv_obj3x) // c*uc*dt/dx = 1*uc*1/1
                                   + adv1st(f3y, v3, lv_obj3y)
                                   + adv1st(f3z, w3, lv_obj3z);

                // source term
                const T0  source  = flag_source ? (sc_f[ID(0,0,0)] * dt) : (T0)0.0;

                const T0 _fn =  f[ID(0,0,0)] + df_conv + source;

                fn[ID(0,0,0)] =  _fn;
            };

            // update
            for(int n=0; n<n_scalars; n++) {
                int idx = n * nn_max;

                advection_f(&scalar_n[idx], &scalar[idx], &sc_scalar[idx]);
            }
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


template<int NX_LEAF, typename T0, typename T1>
__HD__
void  ScalarConvectionweno(
    const int*  id_tasks,
    const int   num_tasks,
    const int*  mesh_offsets3x3x3,
    const int   nn_max,
    const int   n_scalars,
          T0*   scalar_n,
    const T0*   scalar,
    const T0*   u,
    const T0*   v,
    const T0*   w,
    const T0*   lv_obj,
    const T0*   sc_scalar,
    const bool  flag_source,
    const T1    dx,
    const T1    dt
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
            // index
            auto ID = [&](int _i, int _j, int _k) { return IndexMG::id<NX_LEAF>(i+_i, j+_j, k+_k, offset3d); };

            if ( FuncObj::is_obj(lv_obj[ID(0,0,0)]) ) {
                for(int n=0; n<n_scalars; n++) {
                    int idx = n * nn_max;

                    scalar_n[idx + ID(0, 0, 0)] = NSOperator::average_neighbor_w_obj<NX_LEAF>(offset3d, i,j,k, &scalar[idx], lv_obj);
                }
                SKIP_FOR();
            }
        
            const T0  u3[] = { u[ID(-1,  0,  0)], u[ID(0,  0,  0)], u[ID( 1,  0,  0)] };
            const T0  v3[] = { v[ID( 0, -1,  0)], v[ID(0,  0,  0)], v[ID( 0,  1,  0)] };
            const T0  w3[] = { w[ID( 0,  0, -1)], w[ID(0,  0,  0)], w[ID( 0,  0,  1)] };

            const T0  lv_obj3x[] = { lv_obj[ID(-1,  0,  0)], lv_obj[ID(0,  0,  0)], lv_obj[ID( 1,  0,  0)] };
            const T0  lv_obj3y[] = { lv_obj[ID( 0, -1,  0)], lv_obj[ID(0,  0,  0)], lv_obj[ID( 0,  1,  0)] };
            const T0  lv_obj3z[] = { lv_obj[ID( 0,  0, -1)], lv_obj[ID(0,  0,  0)], lv_obj[ID( 0,  0,  1)] };

            auto fluxweno = [&](auto... args){ return FuncDifferentialOperator::flux_euler_weno(args..., 1.0, 1.0); };
//            auto fluxweno = [&](T0 f0, T0 f1, T0 f2, T0 f3, T0 f4, T0 f5, T0 vel){ return FuncDifferentialOperator::flux_euler_weno(f0, f1, f2, f3, f4, f5, vel, 1.0, 1.0); };

            auto advweno = [&](const T0 f7[], const T0 vel3[], const T0 lv_obj3[]){
                const T0  us[] = { ( FuncObj::is_obj( lv_obj3[0] ) ) ? (T0)0.0 : (vel3[0] + vel3[1])*(T0)0.5,
                                   ( FuncObj::is_obj( lv_obj3[2] ) ) ? (T0)0.0 : (vel3[1] + vel3[2])*(T0)0.5 };
            
                const T0 ufs[] = { fluxweno(f7[0], f7[1], f7[2], f7[3], f7[4], f7[5], us[0]),
                                   fluxweno(f7[1], f7[2], f7[3], f7[4], f7[5], f7[6], us[1]) };
            
                const T0  dft  = - ( ufs[1] - ufs[0] );
            
                return  dft;
            };

            auto advection_f = [&](T0* fn, const T0* f, const T0* sc_f) {
                const T0  f7x[] = { f[ID(-3,  0,  0)], f[ID(-2,  0,  0)], f[ID(-1,  0,  0)], f[ID(0,  0,  0)], f[ID( 1,  0,  0)], f[ID( 2,  0,  0)], f[ID( 3,  0,  0)] };
                const T0  f7y[] = { f[ID( 0, -3,  0)], f[ID( 0, -2,  0)], f[ID( 0, -1,  0)], f[ID(0,  0,  0)], f[ID( 0,  1,  0)], f[ID( 0,  2,  0)], f[ID( 0,  3,  0)] };
                const T0  f7z[] = { f[ID( 0,  0, -3)], f[ID( 0,  0, -2)], f[ID( 0,  0, -1)], f[ID(0,  0,  0)], f[ID( 0,  0,  1)], f[ID( 0,  0,  2)], f[ID( 0,  0,  3)] };


                // convection //
                const T0  df_conv =  advweno(f7x, u3, lv_obj3x) // c*uc*dt/dx = 1*uc*1/1
                                   + advweno(f7y, v3, lv_obj3y)
                                   + advweno(f7z, w3, lv_obj3z);

                // source term
                const T0  source  = flag_source ? (sc_f[ID(0,0,0)] * dt) : (T0)0.0;

                const T0 _fn =  f[ID(0,0,0)] + df_conv + source;

                fn[ID(0,0,0)] =  _fn;
            };

            // update
            for(int n=0; n<n_scalars; n++) {
                int idx = n * nn_max;

                advection_f(&scalar_n[idx], &scalar[idx], &sc_scalar[idx]);
            }
        } // FOR_EACH3D
    } // FOR_EACH1D_BLOCKIDX
}


} // namespace  NSKernel


#endif
