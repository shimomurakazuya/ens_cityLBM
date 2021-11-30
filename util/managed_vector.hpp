#ifndef MANAGED_VECTOR
#define MANAGED_VECTOR

#include <vector>

#ifdef USE_NVCC
#include "cu_allocator.hpp"
#endif

namespace util {

#ifdef USE_NVCC
template<class T> using managed_vector = std::vector<T, cu_managed_allocator<T>>;
#else
template<class T> using managed_vector = std::vector<T>;
#endif

} // namespace util

#endif
