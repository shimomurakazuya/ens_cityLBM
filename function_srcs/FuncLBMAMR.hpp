#include "FuncLBMAMR.h"
#include "FuncLBM.h"
#include "IndexLBM.h"


namespace  FuncLBMAMR {


inline
__HOST__ __DEVICE__
real  feqAMR_L2X(
    const real  fs,
    const int   idv,
    const real  rhos,
    const real  us,
    const real  vs,
    const real  ws,
    const real  coefAMR
    )
{
    int  iv, jv, kv;
    IndexLBM::lbm_index(iv, jv, kv, idv);

    const real  fs_eq  = FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);
    const real  fs_deq = fs - fs_eq;

    return  fs_eq + coefAMR*fs_deq;
}


inline
__HOST__ __DEVICE__
void  feqAMR_L2X(
          real  fsAMR[],
    const real  fs[],
    const real  coefAMR
    )
{
    real  rhos, us, vs, ws;
    FuncLBM::velocity_lbm(fs, rhos, us, vs, ws);

    for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
            for (int iv=-1; iv<=1; iv++) {
                const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                const real  fs_eq  = FuncLBM::feq_D3Q27(rhos, us, vs, ws, iv, jv, kv);
                const real  fs_deq = fs[idv_leaf] - fs_eq;

                fsAMR[idv_leaf] = fs_eq + coefAMR*fs_deq;
            }
        }
    }
}


};
