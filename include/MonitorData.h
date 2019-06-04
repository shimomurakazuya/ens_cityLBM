#pragma once
#ifndef MONITORDATA_H_
#define MONITORDATA_H_

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <fstream>
#include <mpi.h>

#include <limits>

#include "defineAMR.h"
#include "Vector3d.h"
#include "Parser.h"
#include "Grid.h"
#include "Tree.h"
#include "Parameters.h"
#include "MeshValue.h"

#include "FuncAMRMesh.h"

class MonitorData{
    private:
        float x_;
        float y_;
        float z_;
        int init_rank_;
        int fileno_;

        float u_;
        float v_;
        float w_;
        float levelset_obj_;
        float scalar_;
        float T_;

        int iscalc_;
        int l_;
        int i_;
        int j_;
        int k_;

        float dist_coef_[3];

    public:
        void CoordinatesToCell(const Tree& tree,const MeshValue* meshValues,const int rank);
        void Approximatecalculation(const Tree& tree,const Parameters& parameters,const MeshValue* meshValues);


        MonitorData(){
            x_=y_=z_=std::numeric_limits<float>::quiet_NaN();
            l_=i_=j_=k_=-1;
            init_rank_ =-1;
            fileno_    =-1;
		    iscalc_ = 0;
            u_=v_=w_=std::numeric_limits<float>::quiet_NaN();
            levelset_obj_=std::numeric_limits<float>::quiet_NaN();
            scalar_=std::numeric_limits<float>::quiet_NaN();
            T_=std::numeric_limits<float>::quiet_NaN();
        }

	void set(float x, float y, float z)             { x_=x; y_=y; z_=z; }
	void set(int fileno, float x, float y, float z) { set(x,y,z); fileno_ = fileno; }

	float x() { return x_; }
	float y() { return y_; }
	float z() { return z_; }

	void
    setvalue(
        float u, float v, float w,
        float levelset_obj,
        float scalar,
        float T,
        int cflg
        )
    {
		u_ = u;
        v_ = v;
        w_ = w;
        levelset_obj_ = levelset_obj;
        scalar_       = scalar;
        T_            = T;
        iscalc_=cflg;
	}

	void
    setvalue(
        float u, float v, float w,
        float levelset_obj,
        float scalar,
        float T,
        int cflg,
        int rank
        )
    {
        setvalue(u,v,w, levelset_obj, scalar, T, cflg);
        init_rank_ = rank;
	}

	float u() const { return u_; }
	float v() const { return v_; }
	float w() const { return w_; }
    float levelset_obj() const { return levelset_obj_; }
    float scalar()       const { return scalar_; }
    float T()            const { return T_; }

	int iscalc() const { return iscalc_; }
	int rank()   const { return init_rank_; }
	int fileno() const { return fileno_; }

};


#endif
