#include "Index.h"
#include "defineAMR.h"


namespace  Index {


inline
__HOST__ __DEVICE__
int  idv(const int iv, const int jv, const int kv)
{
    return  ( (iv+1) + 3*(jv+1) + 9*(kv+1) );
}


template<int iv, int jv, int kv>
__HOST__ __DEVICE__
int  idv()
{
    return  ( (iv+1) + 3*(jv+1) + 9*(kv+1) );
}


inline
__HOST__ __DEVICE__
int  idv0()
{
    return  13;
//    return  ( (1) + 3*(1) + 9*(1) );
//    return  idv(0,0,0);
}


inline
__HOST__ __DEVICE__
int  id(
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
int id(
    const int  ids[],
    const int  in, // neighbor access //
    const int  jn,
    const int  kn
    )
{
    return  ids[ idv(in,jn,kn) ];
}


template<int in, int jn, int kn>
__HOST__ __DEVICE__
int id(const int  ids[])
{
    return  ids[ idv<in,jn,kn>() ];
}


inline
__HOST__ __DEVICE__
int id0(const int  ids[])
{
    return  ids[ idv0() ];
}


inline
__HOST__ __DEVICE__
int  id0(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offset0
    )
{
    return  id(ix,iy,iz) + offset0;
}


inline
__HOST__ __DEVICE__
int  id0(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offsets3d[]
    )
{
    const int  offset0 = offsets3d[ idv0() ];
    return  id0(ix,iy,iz, offset0);
}


inline
__HOST__ __DEVICE__
int  id(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   offsets3d[]
    )
{
    constexpr int   nx_leaf = DefAMR::NX_LEAF;

    // inside //
    const int  ix_tmp = (ix+nx_leaf)%nx_leaf;
    const int  iy_tmp = (iy+nx_leaf)%nx_leaf;
    const int  iz_tmp = (iz+nx_leaf)%nx_leaf;
    const int  id_tmp = id(ix_tmp, iy_tmp, iz_tmp);

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

    const int  offset_tmp  = offsets3d[ idv(nix_tmp, niy_tmp, niz_tmp) ];

    // global offset //
    return  id_tmp + offset_tmp;
}


};
