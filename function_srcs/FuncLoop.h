#pragma once
#ifndef FUNCLOOP_H_
#define FUNCLOOP_H_


#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <functional>
#include <array>


namespace  FuncLoop {


inline
int  id(
    const int   ix,
    const int   iy,
    const int   iz,
    const int   nx,
    const int   ny,
    const int   nz
    )
{
    return  ix + nx*iy + nx*ny*iz;
}


inline
std::array<int, 3>  strides(
    const int   nx,
    const int   ny,
    const int   nz
    )
{
    const std::array<int, 3>  ar = { 1, nx, nx*ny };
    return  ar;
}


struct Loop1d {
    const int  ibegin_;
    const int  iend_;

    Loop1d(int ibegin, int iend) : ibegin_(ibegin), iend_(iend) {}

    template<typename Function>
	void  for_each (Function f)
    {
		for (int i=ibegin_; i<iend_; i++) {
			f(i);
		}
    }

    template<typename Function>
	void  for_each_omp (Function f)
    {
#pragma omp parallel for
		for (int i=ibegin_; i<iend_; i++) {
			f(i);
		}
    }

    template<typename Function>
    void  for_each_acc (Function f)
    {
//        for_each(f);
        for_each_omp(f);
    }

    template<typename T, typename Function>
	T  reduce (T init, Function f)
    {
        T total = (T)0;
		for (int i=ibegin_; i<iend_; i++) {
			total += f(i);
		}
        return  total;
    }
};


struct Loop1dT {
    Loop1dT(){}

    template<int ibegin, int iend, typename Function>
	void  for_each (Function f)  const
    {
		for (int i=ibegin; i<iend; i++) {
			f(i);
		}
    }
};


struct Loop3d {
    const int  ibegin_;
    const int  iend_;

    const int  jbegin_;
    const int  jend_;

    const int  kbegin_;
    const int  kend_;

    Loop3d(int ibegin, int iend,
           int jbegin, int jend,
           int kbegin, int kend) :
        ibegin_(ibegin), iend_(iend),
        jbegin_(jbegin), jend_(jend),
        kbegin_(kbegin), kend_(kend)
	{}


    template<typename Function>
    void  for_each (Function f)
    {
        for (int k=kbegin_; k<kend_; k++) {
        for (int j=jbegin_; j<jend_; j++) {
        for (int i=ibegin_; i<iend_; i++) {
            f(i, j, k);
        }
        }
        }
    }

    template<typename Function>
    void  for_each_omp (Function f)
    {
#pragma omp parallel for collapse(3)
        for (int k=kbegin_; k<kend_; k++) {
        for (int j=jbegin_; j<jend_; j++) {
        for (int i=ibegin_; i<iend_; i++) {
            f(i, j, k);
        }
        }
        }
    }

};


struct  Loop3dT {

    template <int ibegin, int iend, int jbegin, int jend, int kbegin, int kend, typename Function>
    void  for_each (Function f)  const
    {
        for (int k=kbegin; k<kend; k++) {
        for (int j=jbegin; j<jend; j++) {
        for (int i=ibegin; i<iend; i++) {
            f(i, j, k);
        }
        }
        }
    }

    template <int ibegin, int iend, int jbegin, int jend, int kbegin, int kend, typename Function>
    void  for_each_omp (Function f)  const
    {
#pragma omp parallel for collapse(3)
        for (int k=kbegin; k<kend; k++) {
        for (int j=jbegin; j<jend; j++) {
        for (int i=ibegin; i<iend; i++) {
            f(i, j, k);
        }
        }
        }
    }

};


struct Loop4d {
    const int  ibegin_;
    const int  iend_;

    const int  jbegin_;
    const int  jend_;

    const int  kbegin_;
    const int  kend_;

    const int  lbegin_;
    const int  lend_;

    Loop4d(int ibegin, int iend,
           int jbegin, int jend,
           int kbegin, int kend,
           int lbegin, int lend) :
        ibegin_(ibegin), iend_(iend),
        jbegin_(jbegin), jend_(jend),
        kbegin_(kbegin), kend_(kend),
        lbegin_(lbegin), lend_(lend)
	{}


    template<typename Function>
    void  for_each (Function f)
    {
        for (int l=lbegin_; l<lend_; l++) {
        for (int k=kbegin_; k<kend_; k++) {
        for (int j=jbegin_; j<jend_; j++) {
        for (int i=ibegin_; i<iend_; i++) {
            f(i, j, k, l);
        }
        }
        }
        }
    }

    template<typename Function>
    void  for_each_omp (Function f)
    {
#pragma omp parallel for collapse(4)
        for (int l=lbegin_; l<lend_; l++) {
        for (int k=kbegin_; k<kend_; k++) {
        for (int j=jbegin_; j<jend_; j++) {
        for (int i=ibegin_; i<iend_; i++) {
            f(i, j, k, l);
        }
        }
        }
        }
    }

};


struct Loop4dLT {
    const int  lbegin_;
    const int  lend_;

    Loop4dLT(int lbegin, int lend) :
        lbegin_(lbegin), lend_(lend)
	{}


    template <int ibegin, int iend, int jbegin, int jend, int kbegin, int kend, typename Function>
    void  for_each (Function f)
    {
        for (int l=lbegin_; l<lend_; l++) {
        for (int k=kbegin; k<kend; k++) {
        for (int j=jbegin; j<jend; j++) {
        for (int i=ibegin; i<iend; i++) {
                f(i, j, k, l);
        }
        }
        }
        }
    }

    template <int ibegin, int iend, int jbegin, int jend, int kbegin, int kend, typename Function>
    void  for_each_omp (Function f)
    {
#pragma omp parallel for collapse(4)
        for (int l=lbegin_; l<lend_; l++) {
        for (int k=kbegin; k<kend; k++) {
        for (int j=jbegin; j<jend; j++) {
        for (int i=ibegin; i<iend; i++) {
                f(i, j, k, l);
        }
        }
        }
        }
    }

    template <int ibegin, int iend, int jbegin, int jend, int kbegin, int kend, typename Function>
    void  for_each_omp_simd (Function f)
    {
#pragma omp parallel for collapse(3)
        for (int l=lbegin_; l<lend_; l++) {
        for (int k=kbegin; k<kend; k++) {
        for (int j=jbegin; j<jend; j++) {
#pragma omp simd
        for (int i=ibegin; i<iend; i++) {
                f(i, j, k, l);
        }
        }
        }
        }
    }

};


struct  LoopD3Q27 {

    template <typename Function>
    void  for_each (Function f)  const
    {
        for (int kv=-1; kv<=1; kv++) {
        for (int jv=-1; jv<=1; jv++) {
        for (int iv=-1; iv<=1; iv++) {
            f(iv, jv, kv);
        }
        }
        }
    }

};


struct  LoopAMRL2F {
    const int  lbegin_;
    const int  lend_;

    LoopAMRL2F(int lbegin, int lend) :
        lbegin_(lbegin), lend_(lend)
	{}

//    template <int nx_leaf, typename Function>
//    void  for_each_omp (Function f)
//    {
//        constexpr int  nx_leafC  = nx_leaf/2;
//
//#pragma omp parallel for collapse(7)
//        for (int l=lbegin_; l<lend_; l++) {
//        for (int kk=0; kk<=1; kk++) {
//        for (int jj=0; jj<=1; jj++) {
//        for (int ii=0; ii<=1; ii++) {
//        for (int k=0; k<nx_leafC; k++) {
//        for (int j=0; j<nx_leafC; j++) {
//        for (int i=0; i<nx_leafC; i++) {
//
//            for (int kkk=0; kkk<=1; kkk++) {
//            for (int jjj=0; jjj<=1; jjj++) {
//            for (int iii=0; iii<=1; iii++) {
//                const int ix = 2*i + iii;
//                const int jx = 2*j + jjj;
//                const int kx = 2*k + kkk;
//
//                f(ix,jx,kx, ii,jj,kk, l);
//            }
//            }
//            }
//
//        }
//        }
//        }
//        }
//        }
//        }
//        }
//    }

    template <int nx_leaf, typename Function>
    void  for_each_omp (Function f)
    {
#pragma omp parallel for collapse(7)
        for (int l=lbegin_; l<lend_; l++) {
        for (int kk=0; kk<=1; kk++) {
        for (int jj=0; jj<=1; jj++) {
        for (int ii=0; ii<=1; ii++) {

        for (int k=0; k<nx_leaf; k++) {
        for (int j=0; j<nx_leaf; j++) {
        for (int i=0; i<nx_leaf; i++) {
                f(i,j,k, ii,jj,kk, l);
            }
            }
            }

        }
        }
        }
        }
    }

};


struct  LoopAMRF2L {
    const int  lbegin_;
    const int  lend_;

    LoopAMRF2L(int lbegin, int lend) :
        lbegin_(lbegin), lend_(lend)
	{}

    template <int nx_leaf, typename Function>
    void  for_each_omp (Function f)
    {
        constexpr int  nx_leafC  = nx_leaf/2;

#pragma omp parallel for collapse(4)
        for (int l=lbegin_; l<lend_; l++) {
        for (int k=0; k<nx_leafC; k++) {
        for (int j=0; j<nx_leafC; j++) {
        for (int i=0; i<nx_leafC; i++) {
            f(i,j,k, l);
        }
        }
        }
        }
    }

};


};


#endif
