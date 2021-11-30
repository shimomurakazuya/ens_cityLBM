#ifndef COORDINATES_H_
#define COORDINATES_H_


#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <type_traits>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"
#include "Tree.h"
#include "FuncAllocate.h"
#include "FuncMath.h"
#include "range.hpp"
#include "defineCal.h"


class  Coordinates {
public:
    using real_type = double;
    using int_type = int; // ssize_t ??
private:
    MemType memType_{MemType::NullPtr};

    // i = (x - x0) / dx, etc; values are located at nodes //
    int*   i_ = nullptr;
    int*   j_ = nullptr;
    int*   k_ = nullptr;

    // x0 == global_domain_min_x, etc //
    real_type  x0_;
    real_type  y0_;
    real_type  z0_;


    real_type  dx_;
    real_type  dy_;
    real_type  dz_;

public:
    Coordinates () {}

    //~Coordinates () {
    //    release();
    //}

public:
    MemType  memType() const { return memType_; }

          auto*  i()        { return  i_; }
          auto*  j()        { return  j_; }
          auto*  k()        { return  k_; }

    const auto*  i()  const { return  i_; }
    const auto*  j()  const { return  j_; }
    const auto*  k()  const { return  k_; }

    __HD__ auto x0()  const { return  x0_; }
    __HD__ auto y0()  const { return  y0_; }
    __HD__ auto z0()  const { return  z0_; }

    __HD__ auto dx()  const { return  dx_; }
    __HD__ auto dy()  const { return  dy_; }
    __HD__ auto dz()  const { return  dz_; }

private: 
    __HD__ auto xnode_by_i(int_type i) const { return x0_ + dx_ * i; }
    __HD__ auto ynode_by_j(int_type j) const { return y0_ + dy_ * j; }
    __HD__ auto znode_by_k(int_type k) const { return z0_ + dz_ * k; }
    __HD__ auto xcell_by_i(int_type i) const { return x0_ + dx_ * (i + (real_type)0.5); }
    __HD__ auto ycell_by_j(int_type j) const { return y0_ + dy_ * (j + (real_type)0.5); }
    __HD__ auto zcell_by_k(int_type k) const { return z0_ + dz_ * (k + (real_type)0.5); }

public:
    __HD__ auto xnode(int id, int_type i_offset=0) const { return xnode_by_i(i_[id] + i_offset); }
    __HD__ auto ynode(int id, int_type j_offset=0) const { return ynode_by_j(j_[id] + j_offset); }
    __HD__ auto znode(int id, int_type k_offset=0) const { return znode_by_k(k_[id] + k_offset); }

    __HD__ auto xcell(int id, int_type i_offset=0) const { return xcell_by_i(i_[id] + i_offset); }
    __HD__ auto ycell(int id, int_type j_offset=0) const { return ycell_by_j(j_[id] + j_offset); }
    __HD__ auto zcell(int id, int_type k_offset=0) const { return zcell_by_k(k_[id] + k_offset); }

public:
    void init(const int  nn_max, const enum MemType  memType)
    {
        memType_ = memType;
    
        allocate(nn_max);
    }

    void set_uniform(const int  lv, const Grid& grid, const Tree& tree)
    {
        // initialize uniform mesh //
        const auto x0 = grid.offset();
        x0_ = x0.x();
        y0_ = x0.y();
        z0_ = x0.z();
        const auto dx = grid.dx() / (double)DefAMR::NX_LEAF;
        dx_ = dx;
        dy_ = dx;
        dz_ = dx;
        
        // memset i, j, k //
        int  id_offset = 0;
        for(const auto& node: tree.nodes()) {
            if(lv != node.level()) { continue; }
            const auto nx_leaf = DefAMR::NX_LEAF;
            static_assert(std::is_same<int_type, decltype(node.x())>::value, "type of Tree::Node::x() should be int_type");
            const auto i0 = node.x() * nx_leaf,
                       j0 = node.y() * nx_leaf,
                       k0 = node.z() * nx_leaf;

            // set_xyz
            for(auto k: util::irange(nx_leaf)) {
            for(auto j: util::irange(nx_leaf)) {
            for(auto i: util::irange(nx_leaf)) {
                size_t index = id_offset + i + j*nx_leaf + k*nx_leaf*nx_leaf;
                i_[index] = i0 + i;
                j_[index] = j0 + j;
                k_[index] = k0 + k;
            }}}
            id_offset += std::pow(nx_leaf, 3);
        }
    }

    void set_x0(real_type x0, real_type y0, real_type z0) {
        std::cout << __PRETTY_FUNCTION__ << ": force updating x0: "
            << "(" << x0_ << " " << y0_ << " " << z0_ << ")"
            << " --> (" << x0 << " " << y0 << " " << z0 << ")"
            << std::endl;
        x0_ = x0; y0_ = y0; z0_ = z0;
    }

    void set_dx(real_type dx, real_type dy, real_type dz) {
        std::cout << __PRETTY_FUNCTION__ << ": force updating dx: "
            << "(" << dx_ << " " << dy_ << " " << dz_ << ")"
            << " --> (" << dx << " " << dy << " " << dz << ")"
            << std::endl;
        dx_ = dx; dy_ = dy; dz_ = dz;
    }

    void copy (const int nn_max, const Coordinates&  other)
    {
        FuncAllocate::copy_values( i_, other.i_, nn_max, memType_ );
        FuncAllocate::copy_values( j_, other.j_, nn_max, memType_ );
        FuncAllocate::copy_values( k_, other.k_, nn_max, memType_ );
    
        x0_ = other.x0_;
        y0_ = other.y0_;
        z0_ = other.z0_;
        dx_ = other.dx_;
        dy_ = other.dy_;
        dz_ = other.dz_;
    }

    void reallocate(const MemType memType, const int nn_max)
    {
        release();

        memType_ = memType;
        allocate(nn_max);
    }

private:

    void  allocate(const int nn_max)
    {
        FuncAllocate::allocate_value(&i_, nn_max, memType_);
        FuncAllocate::allocate_value(&j_, nn_max, memType_);
        FuncAllocate::allocate_value(&k_, nn_max, memType_);
    }

    void  release ()
    {
        FuncAllocate::release_value(i_, memType_);
        FuncAllocate::release_value(j_, memType_);
        FuncAllocate::release_value(k_, memType_);
    }

};


#endif
