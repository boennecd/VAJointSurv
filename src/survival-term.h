#ifndef SURVIVAL_TERM_H
#define SURVIVAL_TERM_H

#include "bases.h"
#include <cfaad/AAD.h>
#include <array>
#include "simple-mat.h"
#include "VA-parameter.h"
#include <stdexcept>

namespace survival {

using joint_bases::basisMixin;
using joint_bases::bases_vector;

/**
 * struct to hold nodes, weights, and the number of nodes. In this namespace,
 * the quadrature that are used are some quadrature rule on the interval (0, 1)
 */
struct node_weight {
  /// the nodes should be between zero and one
  double const * ns, * ws;
  vajoint_uint n_nodes;
};

/**
 * computes the approximate expected cumulative hazard times minus one. That is
 *
 *   E[int_0^t exp(z^T.fixef + b(s)^T.fixef_time + association^T.M(s)u + v)ds]
 *     ~ int_0^t
 *       exp(z^T.fixef + b(s)^T.fixef_time + (association^T, 1).hat(M)(s).VA_mean
 *         + (association^T, 1).hat(M)(s).VA_vcov.hat(M)(s)^T.(association^T, 1)^T / 2)ds
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
        return x + r->n_basis();
      })
  };

  /// the largest basis dimension
  vajoint_uint const max_base_dim
  {
    ([&]{
      vajoint_uint out{b.n_basis()};
      for(auto &b : bases_rng)
        out = std::max(out, b->n_basis());
      return out;
    })()
  };

public:
  expected_cum_hazzard
  (basisMixin const &b, bases_vector const &bases_rng,
   vajoint_uint const n_fixef):
  b(b), bases_rng(bases_rng), n_fixef(n_fixef) { }

  /**
   * evaluates the approximate expected cumulative hazard between
   * lower and upper times minus 1. The last two arguments are
   * for working memory.
   */
  template<class T>
  T operator()(node_weight const &nws, double const lower,
               double const upper, double const *design, T const * fixef,
               T const * fixef_time, T const * association, T const *VA_mean,
               T const * VA_vcov, T * wk_mem, double * dwk_mem) const {
    T out{0};
    double const delta_bound{upper - lower};

    for(vajoint_uint i = 0; i < nws.n_nodes; ++i){
      if(nws.ws[i] == 0)
        continue;

      // compute the term from the time-varying fixed effects
      double const node_val{delta_bound * nws.ns[i] + lower};
      b(dwk_mem, node_val);
      T log_hazard
        {cfaad::dotProd(dwk_mem, dwk_mem + b.n_basis(), fixef_time)};

      // construct the (association^T, 1).hat(M)(s) vector
      T * const association_M = wk_mem;
      wk_mem += n_basis_rng_p1;
      {
        vajoint_uint idx{};
        for(vajoint_uint i = 0; i < bases_rng.size(); ++i){
          (*bases_rng[i])(dwk_mem, node_val);
          for(vajoint_uint j = 0; j < bases_rng[i]->n_basis(); ++j, ++idx)
            association_M[idx] = association[i] * dwk_mem[j];
        }
        association_M[idx] = 1;
      }

      // add the mean log hazard term
      log_hazard += cfaad::dotProd
        (association_M, association_M + n_basis_rng_p1, VA_mean);

      // add the quadratic term
      cfaad::matVecProd(VA_vcov, VA_vcov + n_basis_rng_p1 * n_basis_rng_p1,
                        association_M, association_M + n_basis_rng_p1, wk_mem,
                        false);
      log_hazard += cfaad::dotProd
        (association_M, association_M + n_basis_rng_p1, wk_mem) / 2;

      // add the hazard term multiplied by the weight
      out += nws.ws[i] * exp(log_hazard);
    }

    // compute the fixed effect on the log hazard scale, add it, and return
    return delta_bound * out *
      exp(cfaad::dotProd(design, design + n_fixef, fixef));
  }

  /**
   * returns the needed working memory of eval. The first elements is the
   * required T memory and the second element is the required double memory.
   */
  std::array<vajoint_uint, 2> get_wkmem(){
    return { 2 * n_basis_rng_p1, max_base_dim };
  }
};

/**
 * holds input for the observations of each type of outcome. It is used to
 * ease parsing of data in the constructor of survival_dat.
 */
struct obs_input {
  /// the number of observations with this type of outcome
  vajoint_uint n_obs;
  /// pointer to the lower bounds for possible right-truncation
  double const *lbs;
  /// pointer to the event time or censoring time
  double const *ubs;
  /// pointer to the event indicators (1: an event; otherwise no event)
  double const *event;
};

/// class that holds survival terms in the lower bound
class survival_dat {
  /**
   * the bases for the time-varying fixed effects (one for each type of
   * outcome)
   */
  bases_vector bases_fix;
  /// the bases for the time-varying random effects
  bases_vector bases_rng;
  /// design matrices for the fixed effects (one for each type of outcome)
  std::vector<simple_mat<double> > design_mats;
  /**
   * functor to compute the approximation of the expected cumulative hazard (one
   * for each type of outcome)
   */
  std::vector<expected_cum_hazzard> cum_hazs;
  /// vector of vectors with lower bounds, upper bounds, and outcome indicators
  struct obs_info_obj {
    double lb, ub;
    bool event;
  };
  std::vector<std::vector<obs_info_obj> > obs_info;

public:
  /// the indices of the parameters
  subset_params const par_idx;
  /// the number of type of outcomes
  vajoint_uint const n_outcomes = bases_fix.size();
  /// the required working memory
  std::array<vajoint_uint, 2> wkmem_val
  {
    ([&]{
      std::array<vajoint_uint, 2> out{0,0};
      for(auto &ch : cum_hazs){
        auto ch_mem = ch.get_wkmem();
        out[0] = std::max(out[0], ch_mem[0]);
        out[1] = std::max(out[1], ch_mem[1]);
      }

      vajoint_uint const n_shared_p1{par_idx.n_shared() + 1};
      out[0] += n_shared_p1 * (n_shared_p1 + 1);

      vajoint_uint d_xtra{n_shared_p1};
      for(auto &b : bases_fix)
        d_xtra = std::max(b->n_basis(), d_xtra);
      out[1] += d_xtra;

      return out;
    })()
  };

  survival_dat
    (bases_vector const &bases_fix, bases_vector const &bases_rng,
     std::vector<simple_mat<double> > &design_mats,
     subset_params const &par_idx, std::vector<obs_input> const &input):
    bases_fix{joint_bases::clone_bases(bases_fix)},
    bases_rng{joint_bases::clone_bases(bases_rng)},
    design_mats{design_mats},
    par_idx{par_idx}
  {
    if(par_idx.surv_info().size() != n_outcomes)
      throw std::invalid_argument("surv_info().size() != n_outcomes");
    if(design_mats.size() != n_outcomes)
      throw std::invalid_argument("design_mats.size() != n_outcomes");
    if(bases_fix.size() != n_outcomes)
      throw std::invalid_argument("bases_fix.size() != n_outcomes");
    if(input.size() != n_outcomes)
      throw std::invalid_argument("input.size() != n_outcomes");
    if(bases_rng.size() != par_idx.marker_info().size())
      throw std::invalid_argument("bases_rng.size() != marker_info().size()");

    cum_hazs.reserve(n_outcomes);
    for(vajoint_uint i = 0; i < n_outcomes; ++i)
      cum_hazs.emplace_back(*bases_fix[i], bases_rng, design_mats[i].n_rows());

    obs_info.resize(n_outcomes);
    for(vajoint_uint i = 0; i < n_outcomes; ++i){
      obs_info[i].reserve(input[i].n_obs);
      if(input[i].n_obs != design_mats[i].n_cols())
        throw std::invalid_argument
          ("input[i].n_obs != design_mats[i].n_cols()");

      for(vajoint_uint j = 0; j < input[i].n_obs; ++j)
        obs_info[i].push_back
          (obs_info_obj
            {input[i].lbs[j], input[i].ubs[j], input[i].event[j] == 1});
    }
  }

  /// returns the number of lower bound terms of a given type
  vajoint_uint get_n_terms(const vajoint_uint type){
    return obs_info[type].size();
  }

  /// evaluates the lower bound of observation idx for the type of outcome
  template<class T>
  T operator()
    (T const *param, T *wk_mem, const vajoint_uint idx, const vajoint_uint type,
     double * dwk_mem, node_weight const &nws) const {
    // get the information for the outcome and event type
    obs_info_obj const &info{obs_info[type][idx]};
    expected_cum_hazzard const &haz{cum_hazs[type]};
    auto const &surv_info{par_idx.surv_info()[type]};

    // compute the approximate expected log hazard if needed
    T out{0};
    double const * const design = design_mats[type].col(idx);
    vajoint_uint const n_shared{par_idx.n_shared()},
                    n_shared_p1{n_shared + 1};
    if(info.event){
      out -= cfaad::dotProd(design, design + surv_info.n_fix,
                            param + surv_info.idx_fix);

      // TODO: we can avoid re-computing these bases
      (*bases_fix[type])(dwk_mem, info.ub);
      out -= cfaad::dotProd(dwk_mem, dwk_mem + surv_info.n_variying,
                            param + surv_info.idx_varying);

      vajoint_uint offset{};
      for(size_t i = 0; i < bases_rng.size(); ++i){
        (*bases_rng[i])(dwk_mem, info.ub);
        auto M_VA_mean = cfaad::dotProd
          (dwk_mem, dwk_mem + bases_rng[i]->n_basis(),
           param + par_idx.va_mean() + offset);
        out -= param[surv_info.idx_association + i] * M_VA_mean;

        offset += bases_rng[i]->n_basis();
      }

      // the frailty term
      out -= param[par_idx.va_mean() + n_shared + type];
    }

    // add the term from the approximate expected cumulative hazard
    T * const VA_mean{wk_mem};
    for(vajoint_uint i = 0; i < n_shared; ++i)
      VA_mean[i] = param[par_idx.va_mean() + i];
    VA_mean[n_shared] = param[par_idx.va_mean() + n_shared + type];
    wk_mem += n_shared_p1;

    T * const VA_vcov{wk_mem};
    {
      vajoint_uint const rng_dim{n_shared + par_idx.n_shared_surv()};
      T const * const vcov_full{param + par_idx.va_vcov()};
      for(vajoint_uint j = 0; j < n_shared; ++j){
        for(vajoint_uint i = 0; i < n_shared; ++i)
          VA_vcov[i + j * n_shared_p1] = vcov_full[i + j * rng_dim];

        VA_vcov[n_shared + j * n_shared_p1] =
          vcov_full[n_shared + type + j * rng_dim];
      }

      vajoint_uint const offset{n_shared + type};
      for(vajoint_uint i = 0; i < n_shared; ++i)
        VA_vcov[i + n_shared * n_shared_p1] = vcov_full[i + offset * rng_dim];

      VA_vcov[n_shared + n_shared * n_shared_p1] =
        vcov_full[offset + offset * rng_dim];
    }
    wk_mem += n_shared_p1 * n_shared_p1;

    out += haz
      (nws, info.lb, info.ub, design, param + surv_info.idx_fix,
       param + surv_info.idx_varying, param + surv_info.idx_association,
       VA_mean, VA_vcov, wk_mem, dwk_mem);

    return out;
  }

  /**
   * returns the needed working memory of eval. The first elements is the
   * required T memory and the second element is the required double memory.
   */
  std::array<vajoint_uint, 2> get_wkmem(){
    return wkmem_val;
  }
};

} // namespace survival

#endif
