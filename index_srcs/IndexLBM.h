#pragma once
#ifndef INDEXLBM_H_
#define INDEXLBM_H_


#include "defineCUDA.h"
#include <iostream>


namespace  IndexLBM {


inline
__HOST__ __DEVICE__
void lbm_index(int& iv, int& jv, int& kv, const int idv);


inline
__HOST__ __DEVICE__
int  idv_lbm(const int iv, const int jv, const int kv);


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int  idv_lbm();


inline
__HOST__ __DEVICE__
int  idv_lbm0();


inline
__HOST__ __DEVICE__
int  idv_lbm_leaf(const int iv, const int jv, const int kv);


inline
__HOST__ __DEVICE__
int id0(
    const int  id_lbm_base,
    const int  iv, // lbm velocity index //
    const int  jv,
    const int  kv
    );


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int id0(const int  id_lbm_base);


inline
__HOST__ __DEVICE__
int id0(
    const int  ids_lbm_base[],
    const int  iv, // lbm velocity index //
    const int  jv,
    const int  kv
    );


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int id0(const int  ids_lbm_base[]);


inline
__HOST__ __DEVICE__
int id(
    const int  ids_lbm_base[],
    const int  in, // neighbor access //
    const int  jn,
    const int  kn,
    const int  iv, // lbm velocity index //
    const int  jv,
    const int  kv
    );


template<int in, int jn, int kn, int iv, int jv, int kv>
__HOST__ __DEVICE__
int id(const int  ids_lbm_base[]);


inline
__HOST__ __DEVICE__
int id(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offsets[],
    const int   iv, // lbm velocity index //
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
int id_base0(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offsets[]
    );


inline
__HOST__ __DEVICE__
int id0(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offset0,
    const int   iv, // lbm velocity index //
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
int id0(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offsets[],
    const int   iv, // lbm velocity index //
    const int   jv,
    const int   kv
    );


inline
__HOST__ __DEVICE__
int id0_base0(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offset0
    );


inline
__HOST__ __DEVICE__
int  local_id(
    const int   ix,
    const int   iy,
    const int   iz
    );


inline
__HOST__ __DEVICE__
int  id_mesh(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offsets[]
    );


inline
__HOST__ __DEVICE__
int  id0_mesh(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offset
    );


};


#include "IndexLBM.hpp"


#endif
