#ifndef GHQ_LP_UTILS_H
#define GHQ_LP_UTILS_H

#include <R_ext/RS.h> // for F77_NAME and F77_CALL

namespace ghqCpp {

extern "C" {
  void F77_NAME(dtrmm)
    (const char *side, const char *uplo, const char *transa,
     const char *diag, const int *m, const int *n, const double *alpha,
     const double *A, const int *lda, double *B,
     const int *ldb,
     size_t, size_t, size_t, size_t);

  // TODO: possibly use DTPMV for the matrix-vector calls
}

} // namespace ghqCpp

#endif
