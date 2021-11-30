#pragma once
#ifndef VALUEBUFF_H_
#define VALUEBUFF_H_


#include <iostream>
#include <cstdlib>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"
#include "FuncAllocate.h"


class  ValueBuff {
private:
//    typedef unsigned char   bufftype;
    typedef float           bufftype;

    MemType memType_{MemType::NullPtr};

    int  nQ_;
	bufftype*	sbuff_ = nullptr;
	bufftype*	rbuff_ = nullptr;

	bufftype*	sbuff_host_ = nullptr;
	bufftype*	rbuff_host_ = nullptr;

    const real  ndiv_buff_;

public:
     ValueBuff () : nQ_(27), ndiv_buff_(0.70)
     {
        memType_ = MemType::Host;
#ifdef USE_NVCC
//        memType_ = MemType::Managed;
        memType_ = MemType::Device;
#endif
     }
    ~ValueBuff () { release(); }

public:
    MemType  memType() const { return memType_; }

    template<typename T>       T*  sbuff_cal()       { return  reinterpret_cast<T*>(sbuff_); }
    template<typename T>       T*  rbuff_cal()       { return  reinterpret_cast<T*>(rbuff_); }
    template<typename T> const T*  sbuff_cal() const { return  const_cast<const T*>(reinterpret_cast<T*>(sbuff_)); }
    template<typename T> const T*  rbuff_cal() const { return  const_cast<const T*>(reinterpret_cast<T*>(rbuff_)); }

    template<typename T>       T*  sbuff_host()       { return  reinterpret_cast<T*>(sbuff_host_); }
    template<typename T>       T*  rbuff_host()       { return  reinterpret_cast<T*>(rbuff_host_); }
    template<typename T> const T*  sbuff_host() const { return  const_cast<const T*>(reinterpret_cast<T*>(sbuff_host_)); }
    template<typename T> const T*  rbuff_host() const { return  const_cast<const T*>(reinterpret_cast<T*>(rbuff_host_)); }

public:
    void  init(const int  nn_max)
    {
        allocate(nn_max);
    }

    void  reallocate(const int  nn_max)
    {
        release();
        allocate(nn_max);
    }

private:
    void  allocate(const int nn_max)
    {
        const double factor = ((double)sizeof(real)) / ((double)sizeof(bufftype));
        const long   nn_max_buff = nn_max*nQ_ * ndiv_buff_ * factor;

        FuncAllocate::allocate_value<bufftype>(&sbuff_, nn_max_buff, memType_);
        FuncAllocate::allocate_value<bufftype>(&rbuff_, nn_max_buff, memType_);

        FuncAllocate::allocate_value<bufftype>(&sbuff_host_, nn_max_buff, MemType::Host);
        FuncAllocate::allocate_value<bufftype>(&rbuff_host_, nn_max_buff, MemType::Host);
    }

    void  release()
    {
        FuncAllocate::release_value<bufftype>(sbuff_, memType_);
        FuncAllocate::release_value<bufftype>(rbuff_, memType_);

        FuncAllocate::release_value<bufftype>(sbuff_host_, MemType::Host);
        FuncAllocate::release_value<bufftype>(rbuff_host_, MemType::Host);
    }

};


#endif
