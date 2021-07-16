#ifndef TESTTHAT_WRAPPER_H
#define TESTTHAT_WRAPPER_H

#include <testthat.h>
#include <cmath>
#include <limits.h>

inline bool pass_rel_err
  (double const val, double const truth, double const eps){
  double const thresh =
    std::abs(truth) < eps ? eps : std::abs(truth) * eps;
  return std::abs(val - truth) < thresh;
}

inline bool pass_rel_err(double const val, double const truth){
  return pass_rel_err(val, truth,
                      std::sqrt(std::numeric_limits<double>::epsilon()));
}

#endif
