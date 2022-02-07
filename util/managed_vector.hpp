#ifndef MANAGED_VECTOR
#define MANAGED_VECTOR

#include <vector>

#if defined(USE_NVCC)
  #include "cu_allocator.hpp"
#elif defined(ENABLE_HIP)
  #include "hip_allocator.hpp"
#endif

namespace util {

#if defined(USE_NVCC)
  template<class T> using managed_vector = std::vector<T, cu_managed_allocator<T>>;
#elif defined(ENABLE_HIP)
  template<class T> using managed_vector = std::vector<T, hip_device_allocator<T>>;
#else
  template<class T> using managed_vector = std::vector<T>;
#endif

} // namespace util

#endif
