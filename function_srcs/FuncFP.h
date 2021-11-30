#pragma once
#ifndef FUNCFP_H_
#define FUNCFP_H_


#include "definePrecision.h"


namespace FP {
  
template<typename T0, typename T1> inline __HD__ T0 r(T1 a) { return static_cast<T0>(a); }

#ifdef USE_NVCC
  template<> inline __HD__ float  r<float,  fp16  >(fp16   a) { return FP16toFP32(a); }
  template<> inline __HD__ double r<double, fp16  >(fp16   a) { return static_cast<double>(FP16toFP32(a)); }
  template<> inline __HD__ fp16   r<fp16,   float >(float  a) { return FP32toFP16(a); }
  template<> inline __HD__ fp16   r<fp16,   double>(double a) { return FP64toFP16(a); }
  
  template<> inline __HD__ double r<double, double>(double  a) { return a; }
  template<> inline __HD__ float  r<float , float >(float   a) { return a; }
  template<> inline __HD__ fp16   r<fp16  , fp16  >(fp16    a) { return a; }
#endif


#ifdef USE_NVCC
  inline __HD__ double2 r2(double* a) { return *reinterpret_cast<double2*>(a); }
  inline __HD__ float2  r2(float*  a) { return *reinterpret_cast<float2* >(a); }
  inline __HD__ double2 r2(const double* a) { return r2( const_cast<double*>(a) ); }
  inline __HD__ float2  r2(const float*  a) { return r2( const_cast<float* >(a) ); }

  inline __HD__ float2  r2(fp16*  a)  {
    #if   defined(FLOAT16_CAL)
      return __half22float2( *reinterpret_cast<__nv_half2*>(a) );
    #elif defined(BFLOAT16_CAL)
      return make_float2(a[0],a[1]);
  //    return __bfloat1622float2( *reinterpret_cast<__nv_bfloat162*>(a) );
    #endif
  }

  inline __HD__ float2  r2(const fp16*  a)  { return r2( const_cast<fp16*  >(a) ); }


  inline __HD__ void w2(double* a, double2 b) {
      double2* a2 = reinterpret_cast<double2*>(a);
      a2[0] = b;
  }
  
  inline __HD__ void w2(float* a, float2 b) {
      float2* a2 = reinterpret_cast<float2*>(a);
      a2[0] = b;
  }


  inline __HD__ void w2(fp16* a, float2 b) {
    #if   defined(FLOAT16_CAL)
      auto* a2 = reinterpret_cast<__nv_half2*>(a);
      a2[0] = __float22half2_rn(b);
    #elif defined(BFLOAT16_CAL)
      a[0] = b.x;
      a[1] = b.y;
  
  //    auto* a2 = reinterpret_cast<__nv_bfloat162*>(a);
  //    a2[0] = __float22bfloat162_rn(b);
    #endif
  }


  template<typename T>
  inline __HD__ void w2(double* a, T b1, T b2) {
      double2* a2 = reinterpret_cast<double2*>(a);
      a2[0] = make_double2(b1, b2);
  }
  
  template<typename T>
  inline __HD__ void w2(float* a, T b1, T b2) {
      float2* a2 = reinterpret_cast<float2*>(a);
      a2[0] = make_float2(b1, b2);
  }
  
  
  template<typename T>
  inline __HD__ void w2(fp16* a, T b1, T b2) {
      auto b = make_float2(b1, b2);
      w2(a, b);
  }
#endif


//template<> inline __HD__ double2 r2<double, double>(double* a) {
//    double2* a2 = reinterpret_cast<double2*>a;
//    return a2[0];
//}


//template<> inline __HD__ double2 r2<double2, double>(double* a) { return *reinterpret_cast<double2*>(a); }
//template<> inline __HD__ float2  r2<float2 , float >(float*  a) { return *reinterpret_cast<float2* >(a); }
//
//template<> inline __HD__ float2  r2<float2 , fp16  >(fp16*  a)  {
//    auto b162 = *reinterpret_cast<__nv_bfloat162*>(a);
//    return __bfloat1622float2(b162);
//}
//
//
//inline __HD__ void w2<double*, double2>(double* a, double2 b) {
//    double2* a2 = reinterpret_cast<double2*>(a);
//    a2[0] = b;
//}
//
//inline __HD__ void w2<float*, float2>(float* a, float2 b) {
//    float2* a2 = reinterpret_cast<float2*>(a);
//    a2[0] = b;
//}
//
//inline __HD__ void w2<fp16*, float2>(fp16* a, float2 b) {
//    auto* a2 = reinterpret_cast<__nv_bfloat162*>(a);
//    a2[0] = __float22bfloat162_rn(b);
//}
//

} // FuncFP


#endif
