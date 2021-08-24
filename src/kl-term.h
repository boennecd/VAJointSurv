#ifndef KL_TERM_H
#define KL_TERM_H

#include "VA-parameter.h"
#include <memory>
#include <limits.h>
#include "arma-wrap.h"
#include "wmem.h"

class kl_term {
  subset_params idx;
  vajoint_uint n_vars;
  vajoint_uint n_wmem_v{2 * n_vars * n_vars};

  /// objects used by setup
  arma::mat vcov_inv,
            vcov_inv_chol,
            vcov_surv_inv,
            vcov_surv_inv_chol;
  double eval_constant = std::numeric_limits<double>::quiet_NaN();

  bool has_vcov = false,
       has_vcov_surv = false;

public:
  kl_term(): idx{}, n_vars{} { }
  kl_term(subset_params const &idx);

  /**
   * sets up objects to evaluate the divergence. Has to be called prior to
   * calling eval and grad */
  void setup(double const *param, double *wk_mem);

  /// Evaluates the lower bound term
  double eval(double const *param, double *wk_mem) const;

  /**
   * Evaluates the gradient of the lower bound term and adds it to the result.
   * The lower bound term is returned
   */
  double grad(double *g, double const *param, double *wk_mem) const;

  size_t n_wmem() const {
    return n_wmem_v;
  }
};

#endif
