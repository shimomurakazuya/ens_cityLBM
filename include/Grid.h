#pragma once
#ifndef GRID_H_
#define GRID_H_


#include <iostream>
#include <fstream>
#include "definePrecision.h"
#include "Vector3d.h"
#include "stIOGrid.h"


class  Grid {
private:
    Vector3d<int>   nx_;
    Vector3d<real>  offset_;
    real            dx_;

public:
    Grid () {}
    Grid (const Vector3d<int>& nx, const Vector3d<real>& offset, const real dx) : nx_(nx), offset_(offset), dx_(dx) {}
    Grid (const Grid& grid) { init(grid.nx(), grid.offset(), grid.dx()); }

    ~Grid () {}

public:
    const Vector3d<int >&  nx    ()  const { return  nx_; }
    const Vector3d<real>&  offset()  const { return  offset_; }
          real             dx    ()  const { return  dx_; }

public:
    Vector3d<real>  length() const
    {
        Vector3d<real>  nx_tmp(nx_.x(), nx_.y(), nx_.z());

        return  nx_tmp*dx();
    }

    void  init(const Vector3d<int> nx, const Vector3d<real> offset, real dx)
    {
        nx_     = nx;
        offset_ = offset;
        dx_     = dx;
    }

    void  copy(const Grid& grid);

public:
    void  writeGrid(const std::string  filename) const;
    void  readGrid (const std::string  filename);

    void  cout(const std::string& str) const;
    void  cout() const { cout(""); }

    std::string  filename(const int lv, const int rank, const int step) const ;
private:

};


#endif
