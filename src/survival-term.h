#ifndef SURVIVAL_TERM_H
#define SURVIVAL_TERM_H

#include "bases.h"
#include <cfaad/AAD.h>
#include <array>

namespace survival {

using joint_bases::basisMixin;
using joint_bases::bases_vector;

/**
 * computes the approximate expected cumulative hazard times minus one. That is
 *
 *   E[int_0^t exp(z^T.delta + b(s)^T.omega + alpha^T.M(s)u + v)ds]
 *     ~ int_0^t exp(z^T.delta + b(s)^T.omega + (alpha^T, 1).hat(M)(s).zeta
 *                + (alpha^T, 1).hat(M)(s).Psi.hat(M)(s)^T.(alpha^T, 1)^T / 2)ds
 *
 * where hat(M)(s) is given by
 *
 *    (( M(s)    0 ),
 *     (    0    1 ))
 */
class expected_cum_hazzard {
  /// the base for the time-varying fixed effect
  basisMixin const &b;
  /// the basis for the time-varying random effects
  bases_vector const &bases_rng;
  /// the number of fixed effects
  vajoint_uint const n_fixef;

  /// the number of time-varying random effect basis function plus one
  vajoint_uint const n_basis_rng_p1
  { 1 + std::accumulate(
      bases_rng.begin(), bases_rng.end(), vajoint_uint{},
      [](vajoint_uint x, joint_bases::bases_vector::value_type const &r){
        return x + r->get_n_basis();
      })
  };

  /// the largest basis dimension
  vajoint_uint const max_base_dim
  {
    ([&]{
      vajoint_uint out{b.get_n_basis()};
      for(auto &b : bases_rng)
        out = std::max(out, b->get_n_basis());
      return out;
    })()
  };

public:
  expected_cum_hazzard
  (basisMixin const &b, bases_vector const &bases_rng,
   vajoint_uint const n_fixef):
  b(b), bases_rng(bases_rng), n_fixef(n_fixef) { }

  /// struct to hold nodes, weights, and the number of nodes
  struct node_weight {
    /// the nodes should be between zero and one
    double const * ns, * ws;
    vajoint_uint n_nodes;
  };

  /**
   * evaluates the approximate expected cumulative hazard between
   * lower_bound and upper_bound times minus 1. The last two arguments are
   * for working memory.
   */
  template<class T>
  T eval(node_weight const &nws, double const lower_bound,
         double const upper_bound, double const *z, T const * delta,
         T const * omega, T const * alpha, T const *zeta, T const * Psi,
         T * wk_mem, double * dwk_mem) const {
    T out{0};

    for(vajoint_uint i = 0; i < nws.n_nodes; ++i){
      if(nws.ws[i] == 0)
        continue;
      T log_term{0};
      double const node_val
        {(upper_bound - lower_bound) * nws.ns[i] + lower_bound};
      b(dwk_mem, node_val);

      // compute the term from the time-varying fixed effects
      log_term += cfaad::dotProd(dwk_mem, dwk_mem + b.get_n_basis(),
                                 omega);

      // construct the (alpha^T, 1).hat(M)(s) vector
      T * const alpha_M = wk_mem;
      wk_mem += n_basis_rng_p1;
      {
        vajoint_uint idx{};
        for(vajoint_uint i = 0; i < bases_rng.size(); ++i){
          (*bases_rng[i])(dwk_mem, node_val);
          for(vajoint_uint j = 0; j < bases_rng[i]->get_n_basis(); ++j, ++idx)
            alpha_M[idx] = alpha[i] * dwk_mem[j];
        }
        alpha_M[idx] = 1;
      }

      // add the mean log hazard term
      log_term += cfaad::dotProd(alpha_M, alpha_M + n_basis_rng_p1, zeta);

      // add the quadratic term
      cfaad::matVecProd(Psi, Psi + n_basis_rng_p1 * n_basis_rng_p1,
                        alpha_M, alpha_M + n_basis_rng_p1, wk_mem, false);
      log_term += cfaad::dotProd(alpha_M, alpha_M + n_basis_rng_p1, wk_mem) / 2;

      // add the hazard term multiplied by the weight
      out += nws.ws[i] * exp(log_term);
    }

    // compute the fixed effect on the log hazard scale, add it, and return
    return out * exp(cfaad::dotProd(z, z + n_fixef, delta));
  }

  /**
   * returns the working memory of eval. The first elements is the required T
   * memory and the second element is the required double memory. */
  std::array<vajoint_uint, 2> get_wkmem(){
    return { 2 * n_basis_rng_p1, max_base_dim};
  }
};

} // namespace survival

#endif
