#ifndef LP_JOINT_H
#define LP_JOINT_H

#include "VA-joint-config.h"
#include <algorithm>
#include <numeric>
#include <R_ext/RS.h> // for F77_NAME and F77_CALL

extern "C" {
  void F77_NAME(dgeqp3)
    (int const *m, int const *n, double *A, int const *lda,
     int *jpvt, double *tau, double *work, int const *lwork, int *info);

  void F77_NAME(dormqr)
    (char const *side, char const *trans, int const *m, int const *n,
     int const *k, double const *A, int const *lda, double const *tau,
     double *C, int const *ldc, double *work, int *lwork, int *info,
     size_t, size_t);
}

namespace lp_joint {
/**
 * computes x^TL^TLx where L is an upper triangular matrix stored with n x n
 * elements. The last argument must have n elements
 */
inline double quad_form
(double const * __restrict__ x, double const * __restrict__ l,
 vajoint_uint const dim, double * __restrict__ wk_mem) noexcept {
  std::fill(wk_mem, wk_mem + dim, 0);
  for(vajoint_uint j = 0; j < dim; l += dim, ++j)
    for(vajoint_uint i = 0; i <= j; ++i)
      wk_mem[i] += x[j] * l[i];

  double out(0.);
  for(vajoint_uint i = 0; i < dim; ++i)
    out += wk_mem[i] *  wk_mem[i];
  return out;
}

/**
 * computes vec(A)^T vec(X[1:k + l, 1:k + l]) where A is a k x k matrix
 * and X is a n x n matrix with k < n. l is an offset. Both matrices are
 * assumed to be symmetric.
 */
inline double submat_trace
(double const * __restrict__ a, double const * __restrict__ x,
 vajoint_uint const k, vajoint_uint const n, vajoint_uint const l) noexcept {
  double out(0);
  x += l * n + l;
  for(vajoint_uint j = 0; j < k; ++j, a += k, x += n){
    for(vajoint_uint i = 0; i < j; ++i)
      out += 2 * a[i] * x[i];
    out += a[j] * x[j];
  }
  return out;
}

/// computes x <- x + Ab. It is assumed that A is a n x n matrix.
inline void mat_vec_prod
(double const * __restrict__ a, double const * __restrict__ b,
 double * __restrict__ res, vajoint_uint const n) noexcept {
  for(vajoint_uint j = 0; j < n; ++j)
    for(vajoint_uint i = 0; i < n; ++i)
      res[i] += *a++ * b[j];
}

/**
 * computes A[1:k + l, 1:k + l] <- A[1:k + l, 1:k + l] + v * B where A is a
 * n x n matrix and B is a k x k matrix with n > k.
 */
inline void mat_add
(double * __restrict__ A, double const * __restrict__ B,
 vajoint_uint const k, vajoint_uint const n, vajoint_uint const offset,
 double const v)
noexcept {
  A += offset * (n + 1);
  for(vajoint_uint j = 0; j < k; ++j, A += n - k)
    for(vajoint_uint i = 0; i < k; ++i)
      *A++ += v * *B++;
}

/**
 * simple matrix vector product y <- y + Xx where y is a k matrix and X is
 * k x n matrix.
 */
inline void mat_vec
(double * __restrict__ y, double const * __restrict__ X,
 double const * __restrict__ x, vajoint_uint const k, vajoint_uint const n){
  for(vajoint_uint j = 0; j < n; ++j)
    for(vajoint_uint i = 0; i < k; ++i)
      y[i] += *X++ * x[j];
}

/// computes x <- x + v * b.
inline void add(double * __restrict__ x, double const * b,
                vajoint_uint const dim, double const v) noexcept {
  for(vajoint_uint i = 0; i < dim; ++i)
    *x++ += v * *b++;
}

/// computes X <- X + sign bb^T.
template<bool do_add>
inline void rank_one
  (double * __restrict__ x, double const *b, vajoint_uint const n) noexcept {
  for(vajoint_uint j = 0; j < n; ++j, x += n)
    for(vajoint_uint i = 0; i < n; ++i){
      double const val = b[i] * b[j];
      x[i] += do_add ? val : -val;
    }
}

/**
 * Given to upper triangles of matrices A and B, the function copy them into
 * the upper traingle of the matrix
 *
 *   | A 0 |
 *   | 0 B |
 *
 * where the zeros are matrices with zeros.
 */
inline void copy_block_upper_tri
  (double *res, double const *A, double const *B, vajoint_uint const na,
   vajoint_uint const nb){
  // copy A
  for(vajoint_uint j = 0; j < na; ++j, res += j, A += j)
    std::copy(A, A + j + 1, res);

  // copy B
  for(vajoint_uint j = 0; j < nb; ++j, res += j, B += j){
    std::fill(res, res + na, 0);
    res += na;
    std::copy(B, B + j + 1, res);
  }
}

} // namespace lp_joint

#endif
