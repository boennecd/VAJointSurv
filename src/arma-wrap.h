#ifndef ARMA_WRAP_H
#define ARMA_WRAP_H

#define _USE_MATH_DEFINES
#include <cmath>

#ifndef ARMA_DONT_USE_OPENMP
#define ARMA_DONT_USE_OPENMP 1
#endif

#ifndef DO_CHECKS
#define ARMA_NO_DEBUG
#endif

#define ARMA_DONT_PRINT_ERRORS

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include <RcppArmadillo.h>
#pragma GCC diagnostic pop

#endif
