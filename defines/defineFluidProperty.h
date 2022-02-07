#pragma once
#ifndef DEFINEFLUIDPROPERTYS_H_
#define DEFINEFLUIDPROPERTYS_H_


#include "definePrecision.h"


namespace  world_property {
    constexpr float Earths_rotation = 7.292*1.0e-5;

    #if defined(CORIOLIS_FORCE)
    constexpr float Geodetic_latitude = (33.0 + 37.0/60.0)/180.0*M_PI; // oklahoama state
    #else
    constexpr float Geodetic_latitude = 0.0;
    #endif
}


namespace  fluid_property {
    constexpr real  TemperatureSA = 298.15; // K //
    constexpr real  PressureSA    = 101325.0; // Pa //

    constexpr real  PressureOffset = PressureSA;

    constexpr real  GravitationalAcceleration = -9.80665;
//    constexpr real  GravitationalAcceleration = -1.0;
//    constexpr real  GravitationalAcceleration = 0.0;

//    constexpr real  Temperature0 = 273.15;
    constexpr real  Temperature0 = 273.15 + 30.0;
//    constexpr real  Temperature0 = 300.0;
//    constexpr real  Temperature0 = 312.0;
//    constexpr real  Temperature0 = 288;
    constexpr real  Temperature_bouyancy = Temperature0;
//    constexpr real  Temperature_bouyancy = 0.5;
};


#if 1
namespace  air_property {
    constexpr real  DensitySA = 1.205; // kg / m3 //

    // KViscosity: by DEFINEFLUIDPROP_KVIS
    // or, as defualt, set physical property of the air
    #ifdef DEFINEFLUIDPROP_KVIS
    constexpr real  KViscosity = DEFINEFLUIDPROP_KVIS;
    #else
    constexpr real  KViscosity = 1.512 * 1.0e-5;
    #endif
    constexpr real  Viscosity  = KViscosity*DensitySA;

    constexpr real  Re_tau = 1.0/KViscosity;
    constexpr real  ReNum = 1.0/KViscosity;

    constexpr real  HeatConduction = 0.0257;

    constexpr real  Pr  = 0.71;
    constexpr real  xi = KViscosity / Pr;

//    constexpr real  beta = 1.0/273.15;
    constexpr real  beta = 1.0/fluid_property::Temperature0;

    constexpr real  Cd_pad = 0.125;
//    constexpr real  Cd_pad = 0.10;
};
#endif


#if 0 // natural convection2d //
namespace  air_property {
    constexpr real  DensitySA = 1.0; // kg / m3 //
//    constexpr real  DensitySA = 1.205; // kg / m3 //

//    constexpr real  Viscosity  = 0.00266458251889485; // Ra = 1.0e3
    constexpr real  Viscosity  = 0.000842614977317636; // Ra = 1.0e4
//    constexpr real  Viscosity  = 0.000266458251889485; // Ra = 1.0e5
//    constexpr real  Viscosity  = 0.0000842614977317636; // Ra = 1.0e6
//    constexpr real  Viscosity  = 0.0000266458251889485; // Ra = 1.0e7
//    constexpr real  Viscosity  = 0.00000842614977317636; // Ra = 1.0e8

    constexpr real  KViscosity = Viscosity/DensitySA;

    constexpr real  Pr = 0.71;
    constexpr real  xi = Viscosity / Pr;

    constexpr real  Re_tau = 1.0/KViscosity;
    constexpr real  ReNum = 1.0/KViscosity;

    constexpr real  beta = 0.00102040816326531;

};
#endif


#if 0 // natural convection3d //
namespace  air_property {
    constexpr real  DensitySA = 1.205; // kg / m3 //

    constexpr real  KViscosity = 1.512 * 1.0e-5;
    constexpr real  Viscosity  = KViscosity*DensitySA;

    constexpr real  Pr = 0.71;
    constexpr real  xi = Viscosity / Pr;

    constexpr real  Re_tau = 1.0/KViscosity;
    constexpr real  ReNum = 1.0/KViscosity;

//    constexpr real  beta = 0.00102040816326531;
    constexpr real  beta = 1.0/273.15;

};
#endif


// channel flow //
#if 0
#define CAL_CHANNEL_FLOW__
namespace  air_property {
    constexpr real  DensitySA = 1.0; // kg / m3 //

    // cavity flow //
//    constexpr real ReNum = 1000.0;

    // channel flow //
    constexpr real Re_tau  = 180.0;
//    constexpr real Re_tau  = 395.0;
//    constexpr real Re_tau  = 590.0;
    constexpr real ReNum = Re_tau;

    constexpr real  KViscosity = 1.0/Re_tau;
    constexpr real  Viscosity  = KViscosity*DensitySA;

    constexpr real  HeatConduction = 0.0;
};
#endif


#if 0 // channel flow with heat transfer 3d //
namespace  air_property {
    constexpr real  DensitySA = 1.205; // kg / m3 //

    constexpr real  Re_tau = 150.0;
    constexpr real  ReNum  = Re_tau;

    constexpr real  KViscosity = 1.0/Re_tau;
    constexpr real  Viscosity  = KViscosity*DensitySA;

    constexpr real  Pr = 0.71;
    constexpr real  xi = Viscosity / Pr;

//    constexpr real  beta = 1.0/273.15;
    constexpr real  beta = 3.56 * 1.0e-3;

};
#endif


#endif
