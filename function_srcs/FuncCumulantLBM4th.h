#pragma once
#ifndef FUNCCUMULANTLBM4TH_H_
#define FUNCCUMULANTLBM4TH_H_


#include "defineCUDA.h"
#include <cmath>
#include "definePrecision.h"
#include "FuncMath.h"


namespace  FuncCumulantLBM4th {


inline
__HOST__ __DEVICE__
void  C2CI4th(real CIs[], const real Cs[], const real kAs[], const real omega, const real rho, const real uc, const real vc, const real wc)
{
	const real u = uc;
    const real v = vc;
    const real w = wc;

    const real  Lim = 0.01;
//    const real  Lim = 1.0e6;

    const real absC120p102 =  FuncMath::abs(kAs[7] + kAs[19]) ;
    const real absC210p012 =  FuncMath::abs(kAs[5] + kAs[21]) ;
    const real absC201p021 =  FuncMath::abs(kAs[11] + kAs[15]) ;
    const real absC120m102 =  FuncMath::abs(kAs[7] - kAs[19]) ;
    const real absC210m012 =  FuncMath::abs(kAs[5] - kAs[21]) ;
    const real absC201m021 =  FuncMath::abs(kAs[11] - kAs[15]) ;
    const real absC111p111 =  FuncMath::abs(kAs[13]) ;

    const real  ep  = 1.0e-6;
    const real  omega1_lim = (omega <= 7.0/4.0) ? 7.0/4.0+ep :
                             (omega >= 2.0)     ? 2.0-ep :
                                                  omega;

    const real  omega1 = omega;
//    const real  omega2 = 1.0;

    const real  omega3 =  8*(2*pow(omega1_lim,2)-3*omega1_lim-2)/(7*pow(omega1_lim,2)-14*omega1_lim-8);
    const real  omega4 =  8*(omega1_lim-2)*(4*omega1_lim-7)/(9*pow(omega1_lim,2)-50*omega1_lim+56);
    const real  omega5 =  24*(3*pow(omega1_lim,3)-13*pow(omega1_lim,2)+12*omega1_lim+4)/(29*pow(omega1_lim,3)-130*pow(omega1_lim,2)+152*omega1_lim+48);

    const real  A =  ((-3)*pow(omega1_lim,2)+2*omega1_lim+4)/(5*pow(omega1_lim,2)-7*omega1_lim+2);
    const real  B =  ((-14)*pow(omega1_lim,2)+28*omega1_lim+4)/(15*pow(omega1_lim,2)-21*omega1_lim+6);

    // org //
//    const real  omega3 =  1.0;
//    const real  omega4 =  1.0;
//    const real  omega5 =  1.0;
//    const real  A = 0.0;
//    const real  B = 0.0;

    const real  omega120p102 =  absC120p102*((-omega3)+1)/(absC120p102+Lim*rho)+omega3;
    const real  omega210p012 =  absC210p012*((-omega3)+1)/(absC210p012+Lim*rho)+omega3;
    const real  omega201p021 =  absC201p021*((-omega3)+1)/(absC201p021+Lim*rho)+omega3;
    const real  omega120m102 =  absC120m102*((-omega4)+1)/(absC120m102+Lim*rho)+omega4;
    const real  omega210m012 =  absC210m012*((-omega4)+1)/(absC210m012+Lim*rho)+omega4;
    const real  omega201m021 =  absC201m021*((-omega4)+1)/(absC201m021+Lim*rho)+omega4;
    const real  omega111p111 =  absC111p111*((-omega5)+1)/(absC111p111+Lim*rho)+omega5;


    const real Dxus =  (omega1*(-2*Cs[2] + Cs[6] + Cs[18]) - Cs[2] - Cs[6] - Cs[18] + kAs[0])/(2*rho) ;
    const real Dyvs =  (omega1*Cs[2] - 2*omega1*Cs[6] + omega1*Cs[18] - Cs[2] - Cs[6] - Cs[18] + kAs[0])/(2*rho) ;
    const real Dzws =  -(-omega1*Cs[2] - omega1*Cs[6] + 2*omega1*Cs[18] + Cs[2] + Cs[6] + Cs[18] - kAs[0])/(2*rho) ;
    const real DxvDyus =  -3*omega1*Cs[4]/rho ;
    const real DxwDzus =  -3*omega1*Cs[10]/rho ;
    const real DywDzvs =  -3*omega1*Cs[12]/rho ;


    CIs[ 0 ] =  ((-omega1)+1)*Cs[(0)];
    CIs[ 1 ] =  ((-omega1)+1)*Cs[(1)];
    CIs[ 2 ] =  rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))/2+rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))/2-0.5*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-(omega1-1)*(Cs[(2)]-Cs[(6)])/3-(omega1-1)*(Cs[(2)]-Cs[(18)])/3+kAs[(0)]/3;
    CIs[ 3 ] =  ((-omega1)+1)*Cs[(3)];
    CIs[ 4 ] =  ((-omega1)+1)*Cs[(4)];
    CIs[ 5 ] =  (-(omega210m012-1))*(Cs[(5)]-Cs[(21)])/2-(omega210p012-1)*(Cs[(5)]+Cs[(21)])/2;
    CIs[ 6 ] =  (-rho)*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))/2-0.5*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(Cs[(2)]-Cs[(6)])/3-(omega1-1)*(Cs[(2)]-Cs[(18)])/3+kAs[(0)]/3;
    CIs[ 7 ] =  (-(omega120m102-1))*(Cs[(7)]-Cs[(19)])/2-(omega120p102-1)*(Cs[(7)]+Cs[(19)])/2;
    CIs[ 8 ] =  A*rho*(Dyvs+Dzws)*(omega1-2)/(3*omega1);
    CIs[ 9 ] =  ((-omega1)+1)*Cs[(9)];
    CIs[ 10 ] =  ((-omega1)+1)*Cs[(10)];
    CIs[ 11 ] =  (-(omega201m021-1))*(Cs[(11)]-Cs[(15)])/2-(omega201p021-1)*(Cs[(11)]+Cs[(15)])/2;
    CIs[ 12 ] =  ((-omega1)+1)*Cs[(12)];
    CIs[ 13 ] =  ((-omega111p111)+1)*Cs[(13)];
    CIs[ 14 ] =  B*DywDzvs*rho*(omega1-2)/(6*omega1);
    CIs[ 15 ] =  (omega201m021-1)*(Cs[(11)]-Cs[(15)])/2-(omega201p021-1)*(Cs[(11)]+Cs[(15)])/2;
    CIs[ 16 ] =  B*DxwDzus*rho*(omega1-2)/(6*omega1);
    CIs[ 17 ] =  0;
    CIs[ 18 ] =  rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))/2-rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-0.5*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-(omega1-1)*(Cs[(2)]-Cs[(6)])/3+2*(omega1-1)*(Cs[(2)]-Cs[(18)])/3+kAs[(0)]/3;
    CIs[ 19 ] =  (omega120m102-1)*(Cs[(7)]-Cs[(19)])/2-(omega120p102-1)*(Cs[(7)]+Cs[(19)])/2;
    CIs[ 20 ] =  A*rho*(Dxus+Dzws)*(omega1-2)/(3*omega1);
    CIs[ 21 ] =  (omega210m012-1)*(Cs[(5)]-Cs[(21)])/2-(omega210p012-1)*(Cs[(5)]+Cs[(21)])/2;
    CIs[ 22 ] =  B*DxvDyus*rho*(omega1-2)/(6*omega1);
    CIs[ 23 ] =  0;
    CIs[ 24 ] =  A*rho*(Dxus+Dyvs)*(omega1-2)/(3*omega1);
    CIs[ 25 ] =  0;
    CIs[ 26 ] =  0;
}


inline
__HOST__ __DEVICE__
void  kappa2kappaI4th_inline(real kAIs[], const real kAs[], const real omega, const real rho, const real uc, const real vc, const real wc)
{
	const real u = uc;
    const real v = vc;
    const real w = wc;

    const real  Lim = 0.01;
//    const real  Lim = 1.0e6;

    const real absC120p102 =  FuncMath::abs(kAs[7] + kAs[19]) ;
    const real absC210p012 =  FuncMath::abs(kAs[5] + kAs[21]) ;
    const real absC201p021 =  FuncMath::abs(kAs[11] + kAs[15]) ;
    const real absC120m102 =  FuncMath::abs(kAs[7] - kAs[19]) ;
    const real absC210m012 =  FuncMath::abs(kAs[5] - kAs[21]) ;
    const real absC201m021 =  FuncMath::abs(kAs[11] - kAs[15]) ;
    const real absC111p111 =  FuncMath::abs(kAs[13]) ;

    const real  ep  = 1.0e-6;
    const real  omega1_lim = (omega <= 7.0/4.0) ? 7.0/4.0+ep :
                             (omega >= 2.0)     ? 2.0-ep :
                                                  omega;

    const real  omega1 = omega;
//    const real  omega2 = 1.0;

    const real  omega3 =  8*(2*pow(omega1_lim,2)-3*omega1_lim-2)/(7*pow(omega1_lim,2)-14*omega1_lim-8);
    const real  omega4 =  8*(omega1_lim-2)*(4*omega1_lim-7)/(9*pow(omega1_lim,2)-50*omega1_lim+56);
    const real  omega5 =  24*(3*pow(omega1_lim,3)-13*pow(omega1_lim,2)+12*omega1_lim+4)/(29*pow(omega1_lim,3)-130*pow(omega1_lim,2)+152*omega1_lim+48);

    const real  A =  ((-3)*pow(omega1_lim,2)+2*omega1_lim+4)/(5*pow(omega1_lim,2)-7*omega1_lim+2);
    const real  B =  ((-14)*pow(omega1_lim,2)+28*omega1_lim+4)/(15*pow(omega1_lim,2)-21*omega1_lim+6);

    // org //
//    const real  omega3 =  1.0;
//    const real  omega4 =  1.0;
//    const real  omega5 =  1.0;
//    const real  A = 0.0;
//    const real  B = 0.0;

    const real  omega120p102 =  absC120p102*((-omega3)+1)/(absC120p102+Lim*rho)+omega3;
    const real  omega210p012 =  absC210p012*((-omega3)+1)/(absC210p012+Lim*rho)+omega3;
    const real  omega201p021 =  absC201p021*((-omega3)+1)/(absC201p021+Lim*rho)+omega3;
    const real  omega120m102 =  absC120m102*((-omega4)+1)/(absC120m102+Lim*rho)+omega4;
    const real  omega210m012 =  absC210m012*((-omega4)+1)/(absC210m012+Lim*rho)+omega4;
    const real  omega201m021 =  absC201m021*((-omega4)+1)/(absC201m021+Lim*rho)+omega4;
    const real  omega111p111 =  absC111p111*((-omega5)+1)/(absC111p111+Lim*rho)+omega5;


    const real Dxus =  (omega1*(-2*kAs[2] + kAs[6] + kAs[18]) + kAs[0] - kAs[2] - kAs[6] - kAs[18])/(2*rho) ;
    const real Dyvs =  (omega1*kAs[2] - 2*omega1*kAs[6] + omega1*kAs[18] + kAs[0] - kAs[2] - kAs[6] - kAs[18])/(2*rho) ;
    const real Dzws =  (omega1*kAs[2] + omega1*kAs[6] - 2*omega1*kAs[18] + kAs[0] - kAs[2] - kAs[6] - kAs[18])/(2*rho) ;
    const real DxvDyus =  -3*omega1*kAs[4]/rho ;
    const real DxwDzus =  -3*omega1*kAs[10]/rho ;
    const real DywDzvs =  -3*omega1*kAs[12]/rho ;


    kAIs[ 0 ] =  kAs[(0)];
    kAIs[ 1 ] =  (-kAs[(1)]);
    kAIs[ 2 ] =  rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))/2+rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))/2-0.5*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-(omega1-1)*(kAs[(2)]-kAs[(6)])/3-(omega1-1)*(kAs[(2)]-kAs[(18)])/3+kAs[(0)]/3;
    kAIs[ 3 ] =  (-kAs[(3)]);
    kAIs[ 4 ] =  ((-omega1)+1)*kAs[(4)];
    kAIs[ 5 ] =  (-(omega210m012-1))*(kAs[(5)]-kAs[(21)])/2-(omega210p012-1)*(kAs[(5)]+kAs[(21)])/2;
    kAIs[ 6 ] =  (-rho)*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))/2-0.5*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(kAs[(2)]-kAs[(6)])/3-(omega1-1)*(kAs[(2)]-kAs[(18)])/3+kAs[(0)]/3;
    kAIs[ 7 ] =  (-(omega120m102-1))*(kAs[(7)]-kAs[(19)])/2-(omega120p102-1)*(kAs[(7)]+kAs[(19)])/2;
    kAIs[ 8 ] =  (12*A*pow(rho,2)*(Dyvs+Dzws)*(omega1-2)+omega1*(72*pow((omega1-1),2)*pow(kAs[(4)],2)-(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])))/(36*omega1*rho);
    kAIs[ 9 ] =  (-kAs[(9)]);
    kAIs[ 10 ] =  ((-omega1)+1)*kAs[(10)];
    kAIs[ 11 ] =  (-(omega201m021-1))*(kAs[(11)]-kAs[(15)])/2-(omega201p021-1)*(kAs[(11)]+kAs[(15)])/2;
    kAIs[ 12 ] =  ((-omega1)+1)*kAs[(12)];
    kAIs[ 13 ] =  ((-omega111p111)+1)*kAs[(13)];
    kAIs[ 14 ] =  (B*DywDzvs*pow(rho,2)*(omega1-2)+omega1*(omega1-1)*(12*(omega1-1)*kAs[(4)]*kAs[(10)]-(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*kAs[(12)]))/(6*omega1*rho);
    kAIs[ 15 ] =  (omega201m021-1)*(kAs[(11)]-kAs[(15)])/2-(omega201p021-1)*(kAs[(11)]+kAs[(15)])/2;
    kAIs[ 16 ] =  (B*DxwDzus*pow(rho,2)*(omega1-2)+omega1*(omega1-1)*(12*(omega1-1)*kAs[(4)]*kAs[(12)]+(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])*kAs[(10)]))/(6*omega1*rho);
    kAIs[ 17 ] =  (4*(omega1-1)*(omega111p111-1)*kAs[(4)]*kAs[(13)]+(omega1-1)*((omega120m102-1)*(kAs[(7)]-kAs[(19)])+(omega120p102-1)*(kAs[(7)]+kAs[(19)]))*kAs[(10)]+(omega1-1)*((omega210m012-1)*(kAs[(5)]-kAs[(21)])+(omega210p012-1)*(kAs[(5)]+kAs[(21)]))*kAs[(12)]+((-(omega201m021-1))*(kAs[(11)]-kAs[(15)])+(omega201p021-1)*(kAs[(11)]+kAs[(15)]))*((-3)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])/12+((omega201m021-1)*(kAs[(11)]-kAs[(15)])+(omega201p021-1)*(kAs[(11)]+kAs[(15)]))*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])/12)/rho;
    kAIs[ 18 ] =  rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))/2-rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-0.5*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-(omega1-1)*(kAs[(2)]-kAs[(6)])/3+2*(omega1-1)*(kAs[(2)]-kAs[(18)])/3+kAs[(0)]/3;
    kAIs[ 19 ] =  (omega120m102-1)*(kAs[(7)]-kAs[(19)])/2-(omega120p102-1)*(kAs[(7)]+kAs[(19)])/2;
    kAIs[ 20 ] =  (12*A*pow(rho,2)*(Dxus+Dzws)*(omega1-2)+omega1*(72*pow((omega1-1),2)*pow(kAs[(10)],2)+(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])+4*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])))/(36*omega1*rho);
    kAIs[ 21 ] =  (omega210m012-1)*(kAs[(5)]-kAs[(21)])/2-(omega210p012-1)*(kAs[(5)]+kAs[(21)])/2;
    kAIs[ 22 ] =  (B*DxvDyus*pow(rho,2)*(omega1-2)+omega1*(omega1-1)*(12*(omega1-1)*kAs[(10)]*kAs[(12)]-(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])+4*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*kAs[(4)]))/(6*omega1*rho);
    kAIs[ 23 ] =  (4*(omega1-1)*(omega111p111-1)*kAs[(10)]*kAs[(13)]+(omega1-1)*((-(omega120m102-1))*(kAs[(7)]-kAs[(19)])+(omega120p102-1)*(kAs[(7)]+kAs[(19)]))*kAs[(4)]+(omega1-1)*((omega201m021-1)*(kAs[(11)]-kAs[(15)])+(omega201p021-1)*(kAs[(11)]+kAs[(15)]))*kAs[(12)]+((-(omega210m012-1))*(kAs[(5)]-kAs[(21)])+(omega210p012-1)*(kAs[(5)]+kAs[(21)]))*((-3)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])/12+((omega210m012-1)*(kAs[(5)]-kAs[(21)])+(omega210p012-1)*(kAs[(5)]+kAs[(21)]))*((-3)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(kAs[(2)]-kAs[(6)])-4*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])/12)/rho;
    kAIs[ 24 ] =  (12*A*pow(rho,2)*(Dxus+Dyvs)*(omega1-2)+omega1*(72*pow((omega1-1),2)*pow(kAs[(12)],2)-(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])+4*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])))/(36*omega1*rho);
    kAIs[ 25 ] =  (4*(omega1-1)*(omega111p111-1)*kAs[(12)]*kAs[(13)]+(omega1-1)*((-(omega201m021-1))*(kAs[(11)]-kAs[(15)])+(omega201p021-1)*(kAs[(11)]+kAs[(15)]))*kAs[(10)]+(omega1-1)*((-(omega210m012-1))*(kAs[(5)]-kAs[(21)])+(omega210p012-1)*(kAs[(5)]+kAs[(21)]))*kAs[(4)]+((-(omega120m102-1))*(kAs[(7)]-kAs[(19)])+(omega120p102-1)*(kAs[(7)]+kAs[(19)]))*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])/12+((omega120m102-1)*(kAs[(7)]-kAs[(19)])+(omega120p102-1)*(kAs[(7)]+kAs[(19)]))*((-3)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(kAs[(2)]-kAs[(6)])-4*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])/12)/rho;
    kAIs[ 26 ] =  ((-108)*omega1*rho*((-8)*pow((omega111p111-1),2)*pow(kAs[(13)],2)+((omega120m102-1)*(kAs[(7)]-kAs[(19)])-(omega120p102-1)*(kAs[(7)]+kAs[(19)]))*((omega120m102-1)*(kAs[(7)]-kAs[(19)])+(omega120p102-1)*(kAs[(7)]+kAs[(19)]))+((omega201m021-1)*(kAs[(11)]-kAs[(15)])-(omega201p021-1)*(kAs[(11)]+kAs[(15)]))*((omega201m021-1)*(kAs[(11)]-kAs[(15)])+(omega201p021-1)*(kAs[(11)]+kAs[(15)]))+((omega210m012-1)*(kAs[(5)]-kAs[(21)])-(omega210p012-1)*(kAs[(5)]+kAs[(21)]))*((omega210m012-1)*(kAs[(5)]-kAs[(21)])+(omega210p012-1)*(kAs[(5)]+kAs[(21)])))+omega1*(3456*pow((omega1-1),3)*kAs[(4)]*kAs[(10)]*kAs[(12)]+pow((omega1-1),2)*((-432)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-432*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+432.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+288*(omega1-1)*(kAs[(2)]-kAs[(6)])+288*(omega1-1)*(kAs[(2)]-kAs[(18)])-288*kAs[(0)])*pow(kAs[(12)],2)+pow((omega1-1),2)*((-432)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+864*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+432.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+288*(omega1-1)*(kAs[(2)]-kAs[(6)])-576*(omega1-1)*(kAs[(2)]-kAs[(18)])-288*kAs[(0)])*pow(kAs[(4)],2)+pow((omega1-1),2)*(864*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-432*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+432.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-576*(omega1-1)*(kAs[(2)]-kAs[(6)])+288*(omega1-1)*(kAs[(2)]-kAs[(18)])-288*kAs[(0)])*pow(kAs[(10)],2)+2*(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])+4*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)]))-144*(omega1-1)*(B*DxvDyus*pow(rho,2)*(omega1-2)+omega1*(omega1-1)*(12*(omega1-1)*kAs[(10)]*kAs[(12)]+(3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*((-omega1)+1)*(kAs[(2)]-kAs[(6)])+4*((-omega1)+1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])*kAs[(4)]))*kAs[(4)]-144*(omega1-1)*(B*DxwDzus*pow(rho,2)*(omega1-2)+omega1*(omega1-1)*(12*(omega1-1)*kAs[(4)]*kAs[(12)]-((-6)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+4*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])*kAs[(10)]))*kAs[(10)]-144*(omega1-1)*(B*DywDzvs*pow(rho,2)*(omega1-2)+omega1*(omega1-1)*(12*(omega1-1)*kAs[(4)]*kAs[(10)]+(3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*((-omega1)+1)*(kAs[(2)]-kAs[(6)])-2*((-omega1)+1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])*kAs[(12)]))*kAs[(12)]+(12*A*pow(rho,2)*(Dxus+Dyvs)*(omega1-2)+omega1*(72*pow((omega1-1),2)*pow(kAs[(12)],2)-(3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*((-omega1)+1)*(kAs[(2)]-kAs[(6)])+4*((-omega1)+1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])*((-6)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+4*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])))*(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])-2*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)])-(12*A*pow(rho,2)*(Dxus+Dzws)*(omega1-2)+omega1*(72*pow((omega1-1),2)*pow(kAs[(10)],2)+(3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*((-omega1)+1)*(kAs[(2)]-kAs[(6)])+4*((-omega1)+1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])*(3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))+3*rho*((-omega1)+2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*((-omega1)+1)*(kAs[(2)]-kAs[(6)])-2*((-omega1)+1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])))*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])+(12*A*pow(rho,2)*(Dyvs+Dzws)*(omega1-2)+omega1*(72*pow((omega1-1),2)*pow(kAs[(4)],2)+((-3)*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))+2*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])*(6*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-3*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))+3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-4*(omega1-1)*(kAs[(2)]-kAs[(6)])+2*(omega1-1)*(kAs[(2)]-kAs[(18)])-2*kAs[(0)])))*(3*rho*(omega1-2)*(Dxus*pow(u,2)-Dyvs*pow(v,2))-6*rho*(omega1-2)*(Dxus*pow(u,2)-Dzws*pow(w,2))-3.0*rho*(Dxus*pow(u,2)+Dyvs*pow(v,2)+Dzws*pow(w,2))-2*(omega1-1)*(kAs[(2)]-kAs[(6)])+4*(omega1-1)*(kAs[(2)]-kAs[(18)])+2*kAs[(0)]))/(216*omega1*pow(rho,2));

}


};


#endif
