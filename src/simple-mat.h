#ifndef SIMPLE_MAT_H
#define SIMPLE_MAT_H

#include <memory>
#include "VA-joint-config.h"
#include <algorithm>

/// a simple matrix container in column-major order
template<class T>
class simple_mat {
  std::unique_ptr<T[]> mem{};
  vajoint_uint n_rows_v{}, n_cols_v{};
  T * external{nullptr};

public:
  /// use external memory
  simple_mat
  (T * external, vajoint_uint const  n_rows, vajoint_uint const n_cols):
  n_rows_v{n_rows}, n_cols_v{n_cols}, external{external} { }

  /// allocates memory for a given size matrix
  simple_mat(vajoint_uint const n_rows, vajoint_uint const n_cols):
  mem{new double[n_cols * n_rows]},
  n_rows_v{n_rows}, n_cols_v{n_cols}, external{nullptr} { }

  /// copy constructor which always copies
  simple_mat(const simple_mat &o):
  mem{new T[o.n_rows_v * o.n_cols_v]},
  n_rows_v{o.n_rows_v},
  n_cols_v{o.n_cols_v}
  {
    std::copy(o.begin(), o.end(), mem.get());
  }

  simple_mat& operator=(const simple_mat &o) {
    mem.reset(new T[o.n_rows_v * o.n_cols_v]);
    n_rows_v = o.n_rows_v;
    n_cols_v = o.n_cols_v;
    std::copy(o.begin(), o.end(), mem.get());
    return *this;
  }

  simple_mat(simple_mat &&o):
  mem{std::move(o.mem)},
  n_rows_v{o.n_rows_v},
  n_cols_v{o.n_cols_v},
  external{o.external} { }

  simple_mat& operator=(simple_mat &&o){
    n_rows_v = o.n_rows_v;
    n_cols_v = o.n_cols_v;
    external = o.external;
    mem = std::move(o.mem);
    return *this;
  }

  vajoint_uint n_cols() const {
    return n_cols_v;
  }
  vajoint_uint n_rows() const {
    return n_rows_v;
  }

  /// returns a pointer to the memory
  double * get(){
    return external ? external : mem.get();
  }
  double const * get() const {
    return external ? external : mem.get();
  }

  /// begin and end end functions
  double * begin() {
    return get();
  }
  double const * begin() const {
    return get();
  }
  double * end() {
    return get() + n_rows() * n_cols();
  }
  double const * end() const {
    return get() + n_rows() * n_cols();
  }

  /// returns the pointer to a given column
  double * col(vajoint_uint const idx){
    return get() + idx * n_rows();
  }
  double const * col(vajoint_uint const idx) const {
    return get() + idx * n_rows();
  }
};

#endif
