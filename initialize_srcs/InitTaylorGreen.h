#pragma once
#ifndef INITTAYLORGREEN_H_
#define INITTAYLORGREEN_H_

#include <iostream>
#include <cmath>


class  InitTaylorGreen {
private:
    double  U0_;
    double  L0_;
    double  rho0_;
    double  ReNum_;

public:
    InitTaylorGreen (const double ReNum) :
        U0_(1.0),
        L0_(1.0),
        rho0_(1.0),
        ReNum_(ReNum)
    {}

    ~InitTaylorGreen () {}

public:
    double U0()     const { return U0_; }
    double L0()     const { return L0_; }
    double rho0()   const { return rho0_; }
    double ReNum()  const { return ReNum_; }

    double k0()     const { return 2.0*M_PI/L0(); }
    double vis()    const { return U0()*L0()/ReNum(); }

    double tc()     const { return 1.0/(2.0*vis()*k0()*k0()); }

    double ua(const double x, const double y, const double t)
    {
        return - U0() * cos(k0()*x) * sin(k0()*y) * exp(-1.0*t/tc());
    }

    double va(const double x, const double y, const double t)
    {
        return   U0() * sin(k0()*x) * cos(k0()*y) * exp(-1.0*t/tc());
    }

    double dpa(const double x, const double y, const double t)
    {
        return - (rho0()*U0()*U0()/4.0) * ( cos(2.0*k0()*x) + cos(2.0*k0()*y) ) * exp(-2.0*t/tc());
    }

private:

};


#endif
