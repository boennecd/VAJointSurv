#ifndef WMEM_H
#define WMEM_H

#include "VA-joint-config.h"
#include "cfaad/AAD.h"
#include "simple-mem-stack.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace wmem {
/**
 * sets the working up for some maximum number of threads. This is not thread
 * and possibly invalids all requested memory
 */
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

/// rewind the working memory to the mark for a given thread
void rewind_to_mark(const size_t);

/// rewind the working memory to the mark for this thread
inline void rewind_to_mark(){
  rewind_to_mark(static_cast<size_t>(get_thread_num()));
}

/// sets the mark for a given thread
void set_mark(const size_t);

/// sets the mark for the current thread
inline void set_mark(){
  set_mark(static_cast<size_t>(get_thread_num()));
}

/// rewinds the working memory fo all threads
void rewind_all();

/// clears the working memory fo all threads
void clear_all();

/// returns a pointer with capacity of some given number of doubles
double * get_double_mem(const size_t);

/// returns a pointer with capacity of some given number of Numbers
cfaad::Number * get_Number_mem(const size_t);

/// returns the simple_mem_stack for a given thread
ghqCpp::simple_mem_stack<double> &mem_stack(const size_t);

/// returns the simple_mem_stack for this thread
inline ghqCpp::simple_mem_stack<double> &mem_stack(){
  return mem_stack(get_thread_num());
}
} // namespace wmem

#endif
