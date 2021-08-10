#ifndef VA_JOINT_CONFIG_H
#define VA_JOINT_CONFIG_H

#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

using vajoint_uint = unsigned;

inline int get_thread_num() noexcept {
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0L;
#endif
}

}

#endif
