#pragma once
#ifndef __MEMORYUSAGE_HPP__
#define __MEMORYUSAGE_HPP__

#include <iostream>

inline void print_memory_info() {
  #ifdef USE_NVCC
    // Working on GPUs only
    size_t free_byte, total_byte;
    cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte );
    if(cuda_status != cudaSuccess) {
      std::cout << "Error: cudaMemGetInfo fails, " << cudaGetErrorString(cuda_status) << "\n" << std::endl;
      exit(1);
    }

    auto free_GB  = (double)free_byte/1024./1024./1024.;
    auto total_GB = (double)total_byte/1024./1024./1024.;
    auto used_GB  = total_GB - free_GB;

    std::cout << "free [GB]/total [GB]: " << free_GB << "/" << total_GB << std::endl;
    std::cout << "used [GB]/total [GB]: " << used_GB << "/" << total_GB << std::endl;
  #endif
}

#endif
