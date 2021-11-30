#pragma once
#ifndef GRID_H_
#define GRID_H_


#include <iostream>
#include <fstream>
#include "definePrecision.h"
#include "defineFilenames.h"
#include "Vector3d.h"
#include "stIOGrid.h"


class  Grid {
private:
    Vector3d<int>   nx_;
    Vector3d<double>  offset_;
    double            dx_;

public:
    Grid () {}
    Grid (const Vector3d<int>& nx, const Vector3d<double>& offset, const double dx) : nx_(nx), offset_(offset), dx_(dx) {}
    Grid (const Grid& grid) { init(grid.nx(), grid.offset(), grid.dx()); }

    ~Grid () {}

public:
    const Vector3d<int >&  nx    ()  const { return  nx_; }
    const Vector3d<double>&  offset()  const { return  offset_; }
          double             dx    ()  const { return  dx_; }

public:
    Vector3d<double>  length() const
    {
        Vector3d<double>  nx_tmp(nx_.x(), nx_.y(), nx_.z());

        return  nx_tmp*dx();
    }

    void  init(const Vector3d<int> nx, const Vector3d<double> offset, double dx)
    {
        nx_     = nx;
        offset_ = offset;
        dx_     = dx;
    }

    void  copy(const Grid& grid)
    {
        nx_     = grid.nx();
        offset_ = grid.offset();
        dx_     = grid.dx();
    }

public:
    void  writeGrid(const std::string  filename) const
    {
        const stIOGrid  stIOgrid{
                            nx_.x(), nx_.y(), nx_.z(),
                            offset_.x(), offset_.y(), offset_.z(),
                            dx_
                            };
    
        // write //
        std::ofstream  fout;
        fout.open(filename, std::ios::binary);
        fout.write( ( char * ) &stIOgrid, sizeof( stIOGrid ) );
        fout.close();
    }

    void  readGrid (const std::string  filename)
    {
        stIOGrid  stIOgrid;
    
        // read //
        std::ifstream  fin;
        fin.open(filename, std::ios::binary);
        fin.read( ( char * ) &stIOgrid, sizeof( stIOGrid ) );
        fin.close();
    
    
        // update //
        nx_.    set_xyz( stIOgrid.nx,       stIOgrid.ny,       stIOgrid.nz       );
        offset_.set_xyz( stIOgrid.offset_x, stIOgrid.offset_y, stIOgrid.offset_z );
    
        dx_ = stIOgrid.dx;
    }

    void  cout(const std::string& str) const
    {
        const Vector3d<double>  dx(dx_,dx_,dx_);
    
        std::cout << str << std::endl;
        offset_. cout("offset = ");
        nx_.     cout("nx     = ");
        dx.      cout("dx     = ");
        length().cout("length = ");
    }

    void  cout() const { cout(""); }

    std::string  filename(const int lv, const int rank, const int step) const
    {
        return    Foldernames::io_folder + "/"
                + Filenames::grid_name0
                + "-lv"   + std::to_string(lv)
                + "-rank" + std::to_string(rank)
                + "-step" + std::to_string(step)
                + ".dat";
    }

private:

};


#endif
