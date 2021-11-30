#ifndef UTIL_PORT_HPP_
#define UTIL_PORT_HPP_

#include <cstdint>
#include <cassert>
#include <type_traits>
#include "gcd.hpp"

namespace util {

namespace port {

struct backend {};

struct openmp  : backend {};
struct cuda    : backend {};
struct openacc : backend {};

struct range1d {
    int ibegin, iend;

    range1d(int ibegin, int iend): ibegin(ibegin), iend(iend) {}

    range1d(int n=1): ibegin(0), iend(n) {}

    int n() const { return iend - ibegin; }
};

struct range3d {
    range1d x, y, z;

    range3d(range1d x, range1d y=1, range1d z=1): x(x), y(y), z(z) {}

    range3d(int nx, int ny=1, int nz=1): x(nx), y(ny), z(nz) {}
};

struct policy3d {
    range3d range;
    int nx, ny, nz; // only for cuda
    int nu, nv, nw; // only for cuda

    policy3d(range3d range, int max_threads=256): range(range) {
        init_cuda_spec(max_threads);
    }

    policy3d(int nx, int ny=1, int nz=1, int max_threads=256): range(nx, ny, nz) {
        init_cuda_spec(max_threads);
    }

    void init_cuda_spec(int max_threads) {
        nx = std::gcd(range.x.n(), max_threads);
        ny = std::gcd(range.y.n(), max_threads/nx);
        nz = std::gcd(range.z.n(), max_threads/nx/ny);
        nu = (range.x.n() + nx - 1)/nx;
        nv = (range.y.n() + ny - 1)/ny;
        nw = (range.z.n() + nz - 1)/nz;
    }
};

#define PORT_ARGS3D util::port::policy3d policy, int i, int j, int k

template<class Func, class... Args> 
__global__ void pfor3d_helper_cuda_(policy3d policy, Func func, Args... args) {
    int i = threadIdx.x + blockDim.x * blockIdx.x + policy.range.x.ibegin;
    int j = threadIdx.y + blockDim.y * blockIdx.y + policy.range.y.ibegin;
    int k = threadIdx.z + blockDim.z * blockIdx.z + policy.range.z.ibegin;
    if(i < policy.range.x.iend && j < policy.range.y.iend && k < policy.range.z.iend) {
        func(policy, i, j, k, args...);
    }
} // pfor3d_helper_cuda_

template<class Func, class... Args> 
void pfor3d_helper_openmp_(policy3d policy, Func func, Args... args) {
    #pragma omp parallel for collapse(3)
    for(int k=policy.range.z.ibegin; k<policy.range.z.iend; k++) {
    for(int j=policy.range.y.ibegin; j<policy.range.y.iend; j++) {
    for(int i=policy.range.x.ibegin; i<policy.range.x.iend; i++) {
        func(policy, i, j, k, args...);
    }
    }
    }
} // pfor3d_heloper_openmp_

template<class Func, class... Args> 
void pfor3d_helper_openacc_(policy3d policy, Func func, Args... args) {
    #pragma acc hogehoge
    for(int k=policy.range.z.ibegin; k<policy.range.z.iend; k++) {
    for(int j=policy.range.y.ibegin; j<policy.range.y.iend; j++) {
    for(int i=policy.range.x.ibegin; i<policy.range.x.iend; i++) {
        func(policy, i, j, k, args...);
    }
    }
    }
} // pfor3d_heloper_openacc_

template<class ExecutionPolicy, class Func, class... Args>
void pfor3d(
policy3d policy,
Func func,
Args... args
) {
    if(std::is_same<ExecutionPolicy, openmp>::value) {
        pfor3d_helper_openmp_<Func, Args...> (policy, func, args...);
    } else if(std::is_same<ExecutionPolicy, cuda>::value) {
        pfor3d_helper_cuda_<Func, Args...> 
            <<<
                dim3(policy.nu, policy.nv, policy.nw), 
                dim3(policy.nx, policy.ny, policy.nz)
            >>>
            (policy, func, args...);
    } else if(std::is_same<ExecutionPolicy, openacc>::value) {
        pfor3d_helper_openacc_<Func, Args...> (policy, func, args...);
    } else {
        static_assert(std::is_base_of<backend, ExecutionPolicy>::value, "unexpected Execution Policy");
    }
}

} // namespace port

} // namespace util

#endif

