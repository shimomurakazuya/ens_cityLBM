#pragma once
#ifndef FUNCMAPDATA_H_
#define FUNCMAPDATA_H_


#include <fstream>
#include <string>
#include <cmath>
#include "definePrecision.h"
#include <vector>
#include <cstdint>
#include <utility> // std::pair
#include <functional> // std::function
#include "runtime_error.hpp"


struct MapData {
    public:
    using value_type = int; // onodera .ssv format
    using id_type = std::intptr_t;

    private:
    std::vector<value_type> val_ {0};
    id_type mx_ = 0, my_ = 0;
    real dx_ = NAN;
    real x0_ = NAN, y0_ = NAN;
    real west_ = NAN, east_ = NAN;
    real south_ = NAN, north_ = NAN;

    public: 
    // <-- init
    void load(const std::string& fname) {
        std::ifstream fin(fname);
        if(!fin) { throw STD_RUNTIME_ERROR(std::string("failed to open file ") + fname); }
        fin >> mx_ >> my_;
        val_.resize(mx_*my_);
        //for(id_type j=my_; --j>=0;) { for(id_type i=0; i<mx_; i++) { fin >> val_.at(i + mx_*j); } }
        for(id_type j=0; j<my_; j++) { for(id_type i=0; i<mx_; i++) { fin >> val_.at(i + mx_*j); } }
    }

    void set_geo(real dx, real x0, real y0, real west, real east, real south, real north) {
        dx_ = dx; x0_ = x0; y0_ = y0; west_ = west; east_ = east; south_ = south; north_ = north;
    }
    // -->

    // <-- accessor
    value_type get_height(const id_type ix, const id_type iy) const { 
        return check_index_in_map(ix, iy) ? val_.at(ix + iy * mx_) : 0.0;
    }
    id_type mx() const { return mx_; }
    id_type my() const { return my_; }

    void levelset_map(
        real* lv_obj,
        std::function<double(int)> x,
        std::function<double(int)> y,
        std::function<double(int)> z,
        const int   nx_leaf, 
        const int   lv,
        const real  dx) const;

    double  get_levelset(
        const double xc,
        const double yc,
        const double zc,
        const int    lv,
        const double dx) const;
    // -->
    
    private:
    std::pair<id_type, id_type> convert_geometory_to_index(double x, double y) const {
        return std::make_pair( std::floor((x - x0_)/dx_), std::floor((y - y0_)/dx_ ) );
    }

    bool check_index_in_map(const id_type i, const id_type j) const { return (0 <= i && i < mx_) && (0 <= j && j < my_); }
};


// deleted by refactoring
namespace FuncMapData {
void  levelset_map(
          real* lv_obj,
    const real* x,
    const real* y,
    const real* z,
    const int   nx_leaf,
    const int   lv,
    const real  dx,
    const int*  height_map,
    const int   mx,
    const int   my,
    const real  dx_map,
    const real  map_offset_x,
    const real  map_offset_y,
    const real  domain_min_x,
    const real  domain_min_y,
    const real  domain_max_x,
    const real  domain_max_y
    ) = delete;


double  get_levelset(
    const double xc,
    const double yc,
    const double zc,
    const int    lv,
    const double dx,
    const int*   height_map,
    const int    mx,
    const int    my,
    const double dx_map,
    const double map_offset_x,
    const double map_offset_y,
    const double domain_min_x,
    const double domain_min_y,
    const double domain_max_x,
    const double domain_max_y
    ) = delete;


int  get_height(
    const int* height_map,
    const int mx,
    const int my,
    const int ix,
    const int iy
    ) = delete;


bool  check_index_in_map(
    const int  ix,
    const int  iy,
    const int  nx,
    const int  ny
    ) = delete;


void  convert_geometory_to_index(
          int& ix,
          int& iy,
    const double x,
    const double y,
    const double dx,
    const double map_offset_x,
    const double map_offset_y
    ) = delete;


};


#endif
