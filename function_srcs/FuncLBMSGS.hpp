#include "FuncLBMSGS.h"
#include "defineSGS.h"
#include "FuncLBM.h"


namespace  FuncLBMSGS {


inline
__HOST__ __DEVICE__
real
sgs_viscosity_D3Q27(
    const real  Fcs,
    const real  Csgs,
    const real  rho,
    const real  vis,
    const real  fs_deq[]
    )
{
    const real  SS   = SS_D3Q27(rho, vis, fs_deq);

    // Coherent-structure Smagorinsky model
    return  Csgs * pow(Fcs, real(1.5)) * sqrt((real)2.0*SS);
}


inline
__HOST__ __DEVICE__
real
SS_D3Q27(
    const real  rho,
    const real  vis,
    const real  fs_deq[]
    )
{
    const real  tau = FuncLBM::relaxation_time(vis);
    const real  A   = -(real)1.5/(rho*tau);

    const real  Q11 =  fs_deq[0] + fs_deq[2] + fs_deq[3] + fs_deq[5] + fs_deq[6] + fs_deq[8] + fs_deq[9] + fs_deq[11] + fs_deq[12] + fs_deq[14] + fs_deq[15] + fs_deq[17] + fs_deq[18] + fs_deq[20] + fs_deq[21] + fs_deq[23] + fs_deq[24] + fs_deq[26] ;
    const real  Q22 =  fs_deq[0] + fs_deq[1] + fs_deq[2] + fs_deq[6] + fs_deq[7] + fs_deq[8] + fs_deq[9] + fs_deq[10] + fs_deq[11] + fs_deq[15] + fs_deq[16] + fs_deq[17] + fs_deq[18] + fs_deq[19] + fs_deq[20] + fs_deq[24] + fs_deq[25] + fs_deq[26] ;
    const real  Q33 =  fs_deq[0] + fs_deq[1] + fs_deq[2] + fs_deq[3] + fs_deq[4] + fs_deq[5] + fs_deq[6] + fs_deq[7] + fs_deq[8] + fs_deq[18] + fs_deq[19] + fs_deq[20] + fs_deq[21] + fs_deq[22] + fs_deq[23] + fs_deq[24] + fs_deq[25] + fs_deq[26] ;
    const real  Q12 =  fs_deq[0] - fs_deq[2] - fs_deq[6] + fs_deq[8] + fs_deq[9] - fs_deq[11] - fs_deq[15] + fs_deq[17] + fs_deq[18] - fs_deq[20] - fs_deq[24] + fs_deq[26] ;
    const real  Q13 =  fs_deq[0] - fs_deq[2] + fs_deq[3] - fs_deq[5] + fs_deq[6] - fs_deq[8] - fs_deq[18] + fs_deq[20] - fs_deq[21] + fs_deq[23] - fs_deq[24] + fs_deq[26] ;
    const real  Q23 =  fs_deq[0] + fs_deq[1] + fs_deq[2] - fs_deq[6] - fs_deq[7] - fs_deq[8] - fs_deq[18] - fs_deq[19] - fs_deq[20] + fs_deq[24] + fs_deq[25] + fs_deq[26] ;

    const real  SS = ( ( Q11*Q11 + Q22*Q22 + Q33*Q33 ) + (real)2.0*( Q12*Q12 + Q23*Q23 + Q13*Q13 ) ) * (A*A);

    return  SS;
}


};
