#ifndef __HIP_ALLOCATOR_HPP__
#define __HIP_ALLOCATOR_HPP__
#include <new>
#include <stdexcept>
#include <hip/hip_runtime.h>
#include "hip_safe_call.hpp"

namespace util {

template<class T> 
struct hip_managed_allocator {
  using value_type = T;
  hip_managed_allocator() {}
  template<class U> hip_managed_allocator(const hip_managed_allocator<U>&) {}
  static T* allocate(std::size_t n) {
    T* ret = nullptr;
    if(n > 0) { HIP_SAFE_CALL(hipMallocManaged(&ret, n*sizeof(T))); }
    return ret;
  }
  static void deallocate(T* p, std::size_t n = 0) {
    if(n > 0) { HIP_SAFE_CALL(hipFree(p)); }
  }
};
template<class T, class U> bool operator == (const hip_managed_allocator<T>&, const hip_managed_allocator<U>&) { return true; }
template<class T, class U> bool operator != (const hip_managed_allocator<T>&, const hip_managed_allocator<U>&) { return false; }

template<class T> 
struct hip_device_allocator {
  using value_type = T;
  hip_device_allocator() {}
  template<class U> hip_device_allocator(const hip_device_allocator<U>&) {}
  static T* allocate(std::size_t n) {
    T* ret = nullptr;
    if(n > 0) { HIP_SAFE_CALL(hipMalloc(&ret, n*sizeof(T))); }
    return ret;
  }
  static void deallocate(T* p, std::size_t n = 0) {
    if(n > 0) { HIP_SAFE_CALL(hipFree(p)); }
  }
};
template<class T, class U> bool operator == (const hip_device_allocator<T>&, const hip_device_allocator<U>&) { return true; }
template<class T, class U> bool operator != (const hip_device_allocator<T>&, const hip_device_allocator<U>&) { return false; }

}

#endif
