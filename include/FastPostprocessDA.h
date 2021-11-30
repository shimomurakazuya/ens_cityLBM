#ifndef FASTPOSTPROCESSDA_H_
#define FASTPOSTPROCESSDA_H_

#include "defineAMR.h"
#include "Tree.h"
#include "MeshValue.h"
#include "Parameters.h"
#include "defineFilenames.h"

#include "managed_vector.hpp"
#include "foreach.h"
#include "range.hpp"

#include <string>
#include <iostream>
#include <sys/stat.h>

class FastPostprocessDA {
private:
    const MPICommEnsemble comm_;
    static constexpr real epsilon = 0.01;

public:
    struct st_ids { ssize_t i, j, k, l, lv; float xx, yy, zz; };
    template<class T> using vector_t = util::managed_vector<T>;

private: // data to output
    vector_t<st_ids> ids_;
    vector_t<char>  calcflg_;
    vector_t<int>   d_rank_;
    vector_t<float> u_mean_, v_mean_, w_mean_;
    vector_t<float> vel2_fluc_;
    vector_t<float> x_, y_, z_;

public:
    const auto& calcflg()   const { return calcflg_; }

    const auto& u_mean() const { return u_mean_; }
    const auto& v_mean() const { return v_mean_; }
    const auto& w_mean() const { return w_mean_; }

    const auto& vel2_fluc() const { return vel2_fluc_; }

public:
    FastPostprocessDA() = delete;
    FastPostprocessDA(const MPICommEnsemble comm): comm_(comm) {}

    ~FastPostprocessDA() { }

    size_t monitorcount() const noexcept { return x_.size(); }


    void setupdata(const std::vector<float>& xx, const std::vector<float>& yy, const std::vector<float>& zz, const Tree& tree, const MeshValue* meshValues_ptr[])
    {
        if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    
        readdata(xx, yy, zz);
        mallocdata();
        setup_coordinate_to_cell(tree, meshValues_ptr);
    }


    void setupdata(const std::vector<float>& xx, const std::vector<float>& yy, const std::vector<float>& zz, const Tree& tree, const MeshValue *meshValues)
    { // compat
        const MeshValue* meshValues_ptr[DefAMR::LV_MAX];
        for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
            meshValues_ptr[lv] = &meshValues[lv]; // &(meshValues_new(lv))
        }
        setupdata(xx,yy,zz, tree, meshValues_ptr);
    }


    void InterpolateMonitorData(const Tree& tree, const Parameters& parameters, const MeshValue* meshValues_ptr[], const ValueTimeAverage* valueTimeAverage_ptr[])
    {
        if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
        if(comm_.is_rank0()) { std::cout << " monitor count: " << monitorcount() << std::endl; }
    
        // scan monitor data
        const auto n_monitor = this->x_.size();
        const st_ids* ids = this->ids_.data();
        const char* calcflg = this->calcflg_.data();
        const int* d_rank = this->d_rank_.data();
        float* u_mean = this->u_mean_.data();
        float* v_mean = this->v_mean_.data();
        float* w_mean = this->w_mean_.data();
        float* vel2_fluc = this->vel2_fluc_.data();
    
    
        const float* x = this->x_.data();
        const float* y = this->y_.data();
        const float* z = this->z_.data();
        const auto* mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
        for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
            const auto& meshValue = *(meshValues_ptr[lv]);
            const int   nn_max    = meshValue.nn_max();

            const auto& valueTimeAverage = *(valueTimeAverage_ptr[lv]);
            const auto* u_mean_mesh = valueTimeAverage.u_mean();
            const auto* v_mean_mesh = valueTimeAverage.v_mean();
            const auto* w_mean_mesh = valueTimeAverage.w_mean();
            const auto* vel2_fluc_mesh = valueTimeAverage.vel2_fluc();
    
            const auto c_ref = parameters.c_ref_lbm();
    
            foreach::exec_1d<foreach::opti>(
                n_monitor,
                [=] __HD__ () {
                    FOR_EACH1D_XX(im, n_monitor) {
                        if(calcflg[im] == 0) { return; } // skip noncal leaf
            
                        const auto id = ids[im];
                        if(id.lv != lv) { return; } // kernel called for each level
            
                        // mesh id
                        int offset3d[27];
                        for(int idv=0; idv<27; idv++) {
                            offset3d[idv] = mesh_offsets3x3x3[27*id.l + idv];
                        }
            
                        // interpolate val
                        float _u_mean=0, _v_mean=0, _w_mean=0, _vel2_fluc=0;
                        for(int kk=0; kk<2; kk++) {
                        for(int jj=0; jj<2; jj++) {
                        for(int ii=0; ii<2; ii++) {
                            const auto wx = id.xx*ii + ((ii+1)%2)*(1.f-id.xx);
                            const auto wy = id.yy*jj + ((jj+1)%2)*(1.f-id.yy);
                            const auto wz = id.zz*kk + ((kk+1)%2)*(1.f-id.zz);
                            const auto weight = wx*wy*wz;
            
                            const auto id_mesh = Index::id(id.i+ii, id.j+jj, id.k+kk, offset3d);
                            _u_mean += weight * u_mean_mesh[id_mesh];
                            _v_mean += weight * v_mean_mesh[id_mesh];
                            _w_mean += weight * w_mean_mesh[id_mesh];

                            _vel2_fluc += weight * vel2_fluc_mesh[id_mesh];
                        }}}
            
                        // write final
                        u_mean[im] = _u_mean * c_ref;
                        v_mean[im] = _v_mean * c_ref;
                        w_mean[im] = _w_mean * c_ref;
                        vel2_fluc[im] = _vel2_fluc * c_ref*c_ref;
                    } // FOR_EACH1D_XX
                } // [=] __HD__ ()
            ); // foreach::exec_1d
        } // for lv
    }


    void InterpolateMonitorData(const Tree& tree, const Parameters& parameters, const MeshValue *meshValues, const ValueTimeAverage *valueTimeAverage)
    { // compat
        const MeshValue*        meshValues_ptr      [DefAMR::LV_MAX];
        const ValueTimeAverage* valueTimeAverage_ptr[DefAMR::LV_MAX];
        for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
            meshValues_ptr      [lv] = &meshValues      [lv]; // &(meshValues_new(lv))
            valueTimeAverage_ptr[lv] = &valueTimeAverage[lv];
        }
        InterpolateMonitorData(tree, parameters, meshValues_ptr, valueTimeAverage_ptr);
    }


private:
    void readdata(const std::vector<float>& xx, const std::vector<float>& yy, const std::vector<float>& zz)
    {
        x_.clear(); for(auto&& _v : xx) { x_.push_back(_v); }
        y_.clear(); for(auto&& _v : yy) { y_.push_back(_v); }
        z_.clear(); for(auto&& _v : zz) { z_.push_back(_v); }
    }

    void mallocdata()
    {
        const auto nn = x_.size();
        ids_.resize(nn);
        calcflg_.resize(nn);
        d_rank_.resize(nn);
        u_mean_.resize(nn);
        v_mean_.resize(nn); 
        w_mean_.resize(nn);
        vel2_fluc_.resize(nn);
    
        for(auto&& ii: ids_) { 
            ii.i = ii.j = ii.k = ii.l = ii.lv = -1;
            ii.xx = ii.yy = ii.zz = NAN;
        }
        for(auto&& cc: calcflg_) { cc = 0; }
        for(auto&& dd: d_rank_) { dd = -1; }
    }


public: // to be private, but public is needed for cuda lambda...
    void setup_coordinate_to_cell(const Tree& tree, const MeshValue* meshValues_ptr[])
    {
        if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
        
        const auto rank = comm_.world().rank();
        const auto n_leaf = tree.number_of_nodes();
        const auto nx_leaf = DefAMR::NX_LEAF;
    
        for(auto&& l: util::irange(n_leaf)) {
            const auto& node = tree.nodes().at(l);
            if(node.nodeCalFlags().Cal() == false ) { continue; }
            const int lv = node.level();
    
            const auto& meshValue = *(meshValues_ptr[lv]);
            const auto& coordinates = meshValue.coordinates();
            const auto& offsets = node.neighbor_mesh_offsets();
            const auto offset0 = offsets.offset0();
    
            // coordinates
            const auto dx_cell = meshValue.coordinates().dx();
    
            const auto x0_leaf = coordinates.xnode(offset0);
            const auto y0_leaf = coordinates.ynode(offset0);
            const auto z0_leaf = coordinates.znode(offset0);
            const auto x1_leaf = coordinates.xnode(offset0, nx_leaf);
            const auto y1_leaf = coordinates.ynode(offset0, nx_leaf);
            const auto z1_leaf = coordinates.znode(offset0, nx_leaf);

    
            // search monitor in the leaf
            st_ids* ids = this->ids_.data();
            char* calcflg = this->calcflg_.data();
            int* d_rank = this->d_rank_.data();
            float* x = this->x_.data();
            float* y = this->y_.data();
            float* z = this->z_.data();
            const auto n_monitor = this->x_.size();
            foreach::exec_1d<foreach::opti>(
                n_monitor,
                [=] __HD__ () {
                    FOR_EACH1D_XX(im, n_monitor) {
                        // skip `out of range` or `already found`
                        if(im >= n_monitor) { SKIP_FOR(); }
                        if(calcflg[im] != 0) { SKIP_FOR(); }

                        // search leaf
                        const auto xc = x[im], 
                                   yc = y[im], 
                                   zc = z[im];
                        if(! (      x0_leaf -epsilon*dx_cell <= xc && xc < epsilon*dx_cell + x1_leaf
                                 && y0_leaf -epsilon*dx_cell <= yc && yc < epsilon*dx_cell + y1_leaf
                                 && z0_leaf -epsilon*dx_cell <= zc && zc < epsilon*dx_cell + z1_leaf
                        )) { SKIP_FOR(); }

                        // determine nearest cell; since cell nearest, maybe -0.5 < xx < 0.5
                        const real ii = (xc - x0_leaf) / dx_cell,
                                   jj = (yc - y0_leaf) / dx_cell,
                                   kk = (zc - z0_leaf) / dx_cell;
                        int i = int(ii), j = int(jj), k = int(kk);
                        /*sanity*/ if(i<0||i>=nx_leaf||j<0||j>=nx_leaf||k<0||k>=nx_leaf) { SKIP_FOR(); }
                        const auto xx = ii - i - 0.5, yy = jj - j - 0.5, zz = kk - k - 0.5; // 0.5 == cell-center offset

                        // store result
                        ids[im].i = i;
                        ids[im].j = j;
                        ids[im].k = k;
                        ids[im].l = l;
                        ids[im].lv = lv;
                        ids[im].xx = xx;
                        ids[im].yy = yy;
                        ids[im].zz = zz;
                        calcflg[im] = 1;
                        d_rank[im] = rank;
                    } // [=] __HD__ ()
            }); // foreach
        }
    }

};

#endif
