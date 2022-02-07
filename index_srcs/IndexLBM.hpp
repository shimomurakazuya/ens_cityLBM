#include "IndexLBM.h"
#include "defineAMR.h"
#include "defineLBM.h"


namespace  IndexLBM {


inline
__HOST__ __DEVICE__
void lbm_index(int& iv, int& jv, int& kv, const int idv)
{
    iv = int(idv%3)   - 1;
    jv = int(idv/3)%3 - 1;
    kv = int(idv/9)   - 1;
}


inline
__HOST__ __DEVICE__
int  idv_lbm(const int iv, const int jv, const int kv)
// iv,jv,kv : indexed velocity = -1,0,1 //
{
    return  ( (iv+1) + 3*(jv+1) + 9*(kv+1) );
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int  idv_lbm()
{
    return  ( (iv+1) + 3*(jv+1) + 9*(kv+1) );
}

inline
__HOST__ __DEVICE__
int  idv_lbm0()
{
    return 13;
//    return idv_lbm(0,0,0);
}


inline
__HOST__ __DEVICE__
int  idv_lbm_leaf(const int iv, const int jv, const int kv)
{
    return  idv_lbm(iv,jv,kv) * DefAMR::NN_LEAF;
}


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
    )
{
    return  ids_lbm_base[ idv_lbm(in,jn,kn) ] + idv_lbm(iv,jv,kv) * DefAMR::NN_LEAF;
}


template<int in, int jn, int kn, int iv, int jv, int kv>
__HOST__ __DEVICE__
int id(const int  ids_lbm_base[])
{
    return  ids_lbm_base[ idv_lbm(in,jn,kn) ] + idv_lbm(iv,jv,kv) * DefAMR::NN_LEAF;
}


inline
__HOST__ __DEVICE__
int id0(
    const int  id_lbm_base,
    const int  iv, // lbm velocity index //
    const int  jv,
    const int  kv
    )
{
    return  id_lbm_base + ( (iv+1) + 3*(jv+1) + 9*(kv+1) ) * DefAMR::NN_LEAF;
//    return  id_lbm_base + idv_lbm(iv,jv,kv) * DefAMR::NN_LEAF;
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int id0(const int  id_lbm_base)
{
    return  id_lbm_base + ( (iv+1) + 3*(jv+1) + 9*(kv+1) ) * DefAMR::NN_LEAF;
//    return  id_lbm_base + idv_lbm(iv,jv,kv) * DefAMR::NN_LEAF;
}


inline
__HOST__ __DEVICE__
int id0(
    const int  ids_lbm_base[],
    const int  iv, // lbm velocity index //
    const int  jv,
    const int  kv
    )
{
    return  ids_lbm_base[ idv_lbm0() ] + idv_lbm(iv,jv,kv) * DefAMR::NN_LEAF;
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int id0(const int  ids_lbm_base[])
{
    return  ids_lbm_base[ idv_lbm0() ] + idv_lbm<iv,jv,kv>() * DefAMR::NN_LEAF;
}


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
    )
{
    return  id_mesh(i,j,k, offsets)  +  idv_lbm_leaf(iv, jv, kv);
}


inline
__HOST__ __DEVICE__
int id_base0(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offsets[]
    )
{
//    constexpr int iv = -1;
//    constexpr int jv = -1;
//    constexpr int kv = -1;
//    return  id_mesh(i,j,k, offsets)  +  idv_lbm_leaf(iv, jv, kv);
    return  id_mesh(i,j,k, offsets);
}


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
    )
{
    return  id0_mesh(i,j,k, offset0)  +  idv_lbm_leaf(iv, jv, kv);
}


inline
__HOST__ __DEVICE__
int id0_base0(
    const int   i, // geometory index //
    const int   j,
    const int   k,
    const int   offset0
    )
{
//    constexpr int iv = -1;
//    constexpr int jv = -1;
//    constexpr int kv = -1;
//    return  id0_mesh(i,j,k, offset0)  +  idv_lbm_leaf(iv, jv, kv);
    return  id0_mesh(i,j,k, offset0);
}


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
    )
{
    const int offset0 = offsets[ idv_lbm0() ];
    return id0(i,j,k, offset0, iv,jv,kv);
}


inline
__HOST__ __DEVICE__
int  local_id(
    const int   ix,
    const int   iy,
    const int   iz
    )
{
    constexpr int   nx_leaf = DefAMR::NX_LEAF;

    return  ix + nx_leaf*iy + nx_leaf*nx_leaf*iz;
}


inline
__HOST__ __DEVICE__
int  id_mesh(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offsets[]
    )
{
    constexpr int   nx_leaf = DefAMR::NX_LEAF;
    constexpr int   nQ = LBM_velocity_model::nQ;

    // inside //
    const int  ix_tmp = (ix+nx_leaf)%nx_leaf;
    const int  iy_tmp = (iy+nx_leaf)%nx_leaf;
    const int  iz_tmp = (iz+nx_leaf)%nx_leaf;
    const int  id_tmp = local_id(ix_tmp, iy_tmp, iz_tmp);

    // outside //
    const int  nix_tmp = (ix < 0) ?      -1 :
                         (ix < nx_leaf) ? 0 :
                                          1 ;

    const int  niy_tmp = (iy < 0) ?      -1 :
                         (iy < nx_leaf) ? 0 :
                                          1 ;

    const int  niz_tmp = (iz < 0) ?      -1 :
                         (iz < nx_leaf) ? 0 :
                                          1 ;

    const int  offset_tmp  = offsets[ idv_lbm(nix_tmp, niy_tmp, niz_tmp) ];

    // global offset //
    return  id_tmp + offset_tmp * nQ;
}


inline
__HOST__ __DEVICE__
int  id0_mesh(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offset
    )
{
    constexpr int   nQ = LBM_velocity_model::nQ;
    return  local_id(ix,iy,iz) + offset * nQ;
}


};
