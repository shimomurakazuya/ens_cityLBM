#ifndef MPI_WRAPPER_HPP_
#define MPI_WRAPPER_HPP_

#include <mpi.h>
#include "range.hpp"
#include "mpi_safe_call.hpp"

namespace util {

/// struct mpi: static functions only (just wrapping mpi.h)
struct mpi {
    // basic
    static void init() { MPI_SAFE_CALL(MPI_Init(NULL, NULL)); }
    static void finalize() { MPI_SAFE_CALL(MPI_Finalize()); }
    static void barrier() { MPI_SAFE_CALL(MPI_Barrier(MPI_COMM_WORLD)); }
    static void abort(int exit_code) { MPI_SAFE_CALL(MPI_Abort(MPI_COMM_WORLD, exit_code)); }

    // info
    static int rank() { int ret; MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &ret)); return ret; }
    static int size() { int ret; MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &ret)); return ret; }

    // reduction
    template<typename T> struct MPItypename { static MPI_Datatype name(); }; 
    template<typename T> static T reduce_sum(T t) { T ret = 0; MPI_Allreduce(&t, &ret, 1, MPItypename<T>::name(), MPI_SUM, MPI_COMM_WORLD); return ret; }
    template<typename T> static T reduce_max(T t) { T ret = 0; MPI_Allreduce(&t, &ret, 1, MPItypename<T>::name(), MPI_MAX, MPI_COMM_WORLD); return ret; }
    template<typename T> static T reduce_min(T t) { T ret = 0; MPI_Allreduce(&t, &ret, 1, MPItypename<T>::name(), MPI_MIN, MPI_COMM_WORLD); return ret; }

    // utility
    template<class Func> static void for_each_rank(Func&& func) {
        for(const auto& i: irange(size())) {
            barrier();
            if(rank() == i) { func(); }
        }
        barrier();
    }
    template<class Func> static void on_master(Func&& func) {
        barrier();
        if(rank() == 0) { func(); }
        barrier();
    }
}; // struct mpi

template<> struct mpi::MPItypename<float      > { static MPI_Datatype name() { return MPI_FLOAT      ; } };
template<> struct mpi::MPItypename<double     > { static MPI_Datatype name() { return MPI_DOUBLE     ; } };
template<> struct mpi::MPItypename<long double> { static MPI_Datatype name() { return MPI_LONG_DOUBLE; } };
template<> struct mpi::MPItypename<int        > { static MPI_Datatype name() { return MPI_INT        ; } };
template<> struct mpi::MPItypename<long       > { static MPI_Datatype name() { return MPI_LONG       ; } };
template<> struct mpi::MPItypename<long long  > { static MPI_Datatype name() { return MPI_LONG_LONG  ; } };

/// class mpi_manager: include member var rank, size; init as ctor; finalize as dtor
class mpi_manager {
private:
    int rank_, size_;

public:
    mpi_manager() { 
        mpi::init();
        rank_ = mpi::rank();
        size_ = mpi::size();
    }
    mpi_manager(const mpi_manager&) = delete;
    void operator=(const mpi_manager&) = delete;
    void operator=(mpi_manager&&) = delete;

    ~mpi_manager() { mpi::finalize(); }

    // basic
    void barrier() const { mpi::barrier(); }
    void abort(int exit_code) const { mpi::abort(exit_code); }

    // info
    auto rank() const noexcept  -> decltype(rank_) { return rank_; }
    auto size() const noexcept  -> decltype(size_) { return size_; }

    // reduction
    template<typename T>  T reduce_sum(T t) const { return mpi::reduce_sum(t); }
    template<typename T>  T reduce_max(T t) const { return mpi::reduce_max(t); }
    template<typename T>  T reduce_min(T t) const { return mpi::reduce_min(t); }

    // utility
    template<class Func> void for_each_rank(Func&& func) const { mpi::for_each_rank(func); }
    template<class Func> void on_master(Func&& func) const { mpi::on_master(func); }

}; // struct mpi_manager

}

#endif

