#ifndef VA_JOINT_CONFIG_H
#define VA_JOINT_CONFIG_H

#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

using vajoint_uint = unsigned;

inline int get_thread_num() noexcept {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0L;
#endif
}

/// returns the number of non-zero elements of a triangular matrix
template<class T>
constexpr T dim_tri(T const x){
  return (x * (x + 1)) / 2;
}

#endif
