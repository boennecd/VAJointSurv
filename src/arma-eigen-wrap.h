#ifndef ARMA_EIGEN_WRAP_H
#define ARMA_EIGEN_WRAP_H

// https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2018-June/010034.html
// TODO: possibly a very bad idea
#ifdef NDEBUG
#define HAS_NDEBUG
#endif

#include "arma-wrap.h"

// https://lists.r-forge.r-project.org/pipermail/rcpp-devel/2018-June/010034.html
// TODO: possibly a very bad idea
#ifdef HAS_NDEBUG
#ifndef NDEBUG
#define NDEBUG
#endif
#undef HAS_NDEBUG
#endif

#include <RcppEigen.h>

#endif
