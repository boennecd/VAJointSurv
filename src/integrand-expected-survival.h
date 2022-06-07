#ifndef INTEGRAND_EXPECTED_SURVIVAL_H
#define INTEGRAND_EXPECTED_SURVIVAL_H

#include "ghq.h"

namespace ghqCpp {

/**
 * used to compute  the expected survival. Thus, it computes
 *
 *   g(U) = exp(-sum w_i[i] * exp(eta[i] + (M.U)[i]))
 *
 * given n weights and offsets w and eta and a matrix M of dimension n x R.
 *
 * The derivatives are computed w.r.t. the vector eta and the matrix M
 */
template<bool comp_grad = false>
class expected_survival_term final : public ghq_problem {
  arma::vec const &eta, &weights;
  arma::mat const &M;

  size_t const v_n_vars = M.n_cols,
               v_n_out
    {comp_grad ? 1  + eta.n_elem + M.n_elem : 1};

public:
  expected_survival_term
  (arma::vec const &eta, arma::vec const &weights, arma::mat const &M);

  size_t n_vars() const { return v_n_vars; }
  size_t n_out() const { return v_n_out; }

  void eval
    (double const *points, size_t const n_points, double * __restrict__ outs,
     simple_mem_stack<double> &mem) const;

  double log_integrand
    (double const *point, simple_mem_stack<double> &mem) const;

  double log_integrand_grad
    (double const *point, double * __restrict__ grad,
     simple_mem_stack<double> &mem) const;

  void log_integrand_hess
    (double const *point, double *hess,
     simple_mem_stack<double> &mem) const;
};

} // namespace ghqCpp

#endif
