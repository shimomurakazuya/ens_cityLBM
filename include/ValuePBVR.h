#ifndef VALUEPBVR_H_
#define VALUEPBVR_H_

#ifdef USE_PBVR

#include <iostream>
#include <cstdlib>
#include <vector>

#include "definePrecision.h"
#include "defineMemory.h"
#include "defineAMR.h"
#include "Grid.h"
#include "MeshValue.h"

class Field;

class  ValuePBVR {
public:
    static constexpr int nvariables = 5; // uvw T scalar
    static constexpr int nx = DefAMR::NX_LEAF + 1;
    static constexpr int ny = DefAMR::NX_LEAF + 1;
    static constexpr int nz = DefAMR::NX_LEAF + 1;

private:
    MemType memType_;
    size_t n_leaves_levels_[DefAMR::LV_MAX];
    float* variables_[nvariables]; // [n_variables][n_leaves*nx*ny*nz]
    float* cell_length_; // [n_leaves]
    float* leaf_min_coord_; // [n_leaves]

    struct station { 
        std::string id = "";
        real x = NAN, y = NAN, z = NAN;
        int rank = -99999, l = -99999, i = -99999, j = -99999, k = -99999;
        real ii = -99999, jj = -99999, kk = -99999;
    };
    std::array<station, 12> stations_;

    //// cpu
    //float* variables_cpu_[nvariables]; // [n_variables][n_leaves*nx*ny*nz]
    //float* cell_length_cpu_; // [n_leaves]
    //float* leaf_min_coord_cpu_; // [n_leaves]

public:
    ValuePBVR () {}

    ~ValuePBVR () {
        release();
    }

public:
    void init(const Field* p_field);
    void init_stations(const Field* p_field);

    void update(const Field* p_field);

public:
    MemType  memType() const { return memType_; }

    auto** variables() const noexcept { return const_cast<const float**>(&(variables_[0])); }
    auto*  cell_length()    const noexcept { return cell_length_; }
    auto*  leaf_min_coord() const noexcept { return leaf_min_coord_; }

    auto   n_leaves()  const noexcept { 
        size_t ret = 0;
        for (auto&& n: n_leaves_levels_) {
            ret += n;
        }
        return ret; 
    }

private:
    void  allocate(const size_t n_leaves);
    void  send_gpu_to_cpu();
    void  release();

// call pbvr apis
public:
    void output_pbvr_particles(const Field* p_field, int step);
    void steer_source_posi(Field* p_field, int step);
    void output_scalar_timeseries(Field* p_field, int step);

    void update_pbvr(Field* p_field, int step) {
        output_pbvr_particles(p_field, step);
        steer_source_posi(p_field, step);
        output_scalar_timeseries(p_field, step);
    }
};

#else
class Field;

// wrap by null class
class  ValuePBVR {
public:
    ValuePBVR () {}
    ~ValuePBVR () {}
    void init(const Field*) {}
    void update(const Field*) {};
    void write(const Field*) {};
    MemType  memType() const { return MemType::NullPtr; }

    const float** variables     () const noexcept { return nullptr; }
    const float*  cell_length   () const noexcept { return nullptr; }
    const float*  leaf_min_coord() const noexcept { return nullptr; }
          auto    n_leaves      () const noexcept { return       0; }
};

#endif // ifdef USE_PBVR

#endif // VALUEPBVR_H_

