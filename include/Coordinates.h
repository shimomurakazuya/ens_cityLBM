#pragma once
#ifndef COORDINATES_H_
#define COORDINATES_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"
#include "Tree.h"


class  Coordinates {
private:
    bool  is_initialize_;
    MemType memType_{MemType::NullPtr};

    // values are located at nodes //
    real*   x_ = nullptr;
    real*   y_ = nullptr;
    real*   z_ = nullptr;

    real    dx_;
    real    dy_;
    real    dz_;

public:
    Coordinates () {}

    ~Coordinates () {
        release();
    }

public:
    MemType  memType() const { return memType_; }

          real*  x()        { return  x_; }
          real*  y()        { return  y_; }
          real*  z()        { return  z_; }

    const real*  x()  const { return  x_; }
    const real*  y()  const { return  y_; }
    const real*  z()  const { return  z_; }

    real   x(const int i)  const { return  x()[i]; }
    real   y(const int i)  const { return  y()[i]; }
    real   z(const int i)  const { return  z()[i]; }

    real   dx()  const { return  dx_; }
    real   dy()  const { return  dy_; }
    real   dz()  const { return  dz_; }

    void   reset_dx(const int strides[])
    {
        dx_ = x_[strides[0]] - x_[0];
        dy_ = y_[strides[1]] - y_[0];
        dz_ = z_[strides[2]] - z_[0];
    }

    const real*  xyz(const int d) const
    {
        return  (d == 0) ? x_ :
                (d == 1) ? y_ :
                (d == 2) ? z_ :
                nullptr;
    }

    real  xyz(const int d, const int i) const { return  xyz(d)[i]; }

public:
    void init(const int  nn_max, const enum MemType  memType);

    void set_uniform(const int  lv, const Grid& grid, const Tree& tree);

    void set_x_from_xp(real* x, const real* xp, const int nx, const int ny, const int nz) const;


    void copy (const int nn_max, const Coordinates&  other);

    void reallocate(const MemType memType, const int nn_max)
    {
        release();

        memType_ = memType;
        allocate(nn_max);
    }

private:
    void  set_xyz(const int  offset, const Vector3d<int>&  nx, const Vector3d<real>&  x0, const real  dx);

    void  allocate(const int nn_max);
    void  release ();
};


#endif
