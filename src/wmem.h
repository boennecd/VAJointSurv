#ifndef WMEM_H
#define WMEM_H

#include "VA-joint-config.h"
#include <cfaad/AAD.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace wmem {
/// maximum continuous range of memory
constexpr size_t max_ele{131072};

/// sets the working up for some maximum number of threads
void setup_working_memory(const size_t);

/**
 * rewinds the working memory of a given thread. Rewind must be called often
 * to ensure that not too much memory is allocated.
 */
void rewind(const size_t);

/// rewinds the working memory of this thread
inline void rewind(){
  rewind(static_cast<size_t>(get_thread_num()));
}

/// clears the working memory fo all threads
void clear_all();

/// returns a pointer with capacity of some given number of doubles
double * get_double_mem(const size_t);

/// returns a pointer with capacity of some given number of Numbers
cfaad::Number * get_Number_mem(const size_t);
} // namespace wmem

#endif
