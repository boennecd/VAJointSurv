#ifndef SURVIVAL_TERM_H
#define SURVIVAL_TERM_H

#include <numeric>
#include "bases.h"
#include "cfaad/AAD.h"
#include <array>
#include "simple-mat.h"
#include "VA-parameter.h"
#include <stdexcept>
#include "JointSurv-misc.h"

namespace survival {

using joint_bases::basisMixin;
using joint_bases::bases_vector;

/**
 * computes the approximate expected cumulative hazard times minus one. That is
 *
 *   E[int_0^t exp(z^T.fixef + b(s)^T.fixef_vary + association^T.M(s)u + v)ds]
 *     ~ int_0^t
 *       exp(z^T.fixef + b(s)^T.fixef_vary + (association^T, 1).hat(M)(s).VA_mean
 *         + (association^T, 1).hat(M)(s).VA_vcov.hat(M)(s)^T.(association^T, 1)^T / 2)ds
 *
 * where hat(M)(s) is given by
 *
 *    (( M(s)    0 ),
 *     (    0    1 ))
 *
 * M(s) = hat(M)(s) if there is no frailty
 */
class expected_cum_hazzard {
  /// the base for the time-varying fixed effect
  std::unique_ptr<basisMixin> b;
  /// the number of basis function of b
  vajoint_uint b_n_basis_v{b->n_basis()};
  /// the basis for the time-varying random effects
  bases_vector bases_rng;
  /// the number of basis functions of each rng base expansion
  std::vector<vajoint_uint> rng_n_basis_v{
    ([&]{
      std::vector<vajoint_uint> out;
      out.reserve(bases_rng.size());
      for(auto &b : bases_rng)
        out.emplace_back(b->n_basis());
      return out;
    })()
  };
  /// the number weights of each rng base expansion
  std::vector<vajoint_uint> rng_n_weights_v{
    ([&]{
      std::vector<vajoint_uint> out;
      out.reserve(bases_rng.size());
      for(auto &b : bases_rng)
        out.emplace_back(b->n_weights());
      return out;
    })()
  };
  /// the number of fixed effects
  vajoint_uint n_fixef;
  /// the derivative/integral argument to pass to bases_rng
  std::vector<std::vector<int> > ders_v;
  /// is there a frailty term?
  bool with_frailty_v;

  /// the number of time-varying random effect basis function plus one
  vajoint_uint n_basis_rng_p1
  { 1 + std::accumulate(
      bases_rng.begin(), bases_rng.end(), vajoint_uint{},
      [](vajoint_uint x, joint_bases::bases_vector::value_type const &r){
        return x + r->n_basis();
      })
  };

  /// the largest basis dimension
  vajoint_uint max_base_dim
  {
    ([&]{
      vajoint_uint out{b->n_basis()};
      for(auto &b : bases_rng)
        out = std::max(out, b->n_basis());
      return out;
    })()
  };

  /// the required working memory
  std::array<size_t, 2> n_wmem_v
  {
    ([&]{
      std::array<size_t, 2> out {0, 0};

      out[1] = b->n_wmem();
      for(auto &b : bases_rng)
        out[1] = std::max(out[1], b->n_wmem());

      out[0] += n_basis_rng_p1;
      out[1] += max_base_dim;
      return out;
    })()
  };

  double scale_node_val
    (double const lower, double const upper, double const node) const {
    return (upper - lower) * node + lower;
  }

public:
  expected_cum_hazzard
  (basisMixin const &b_in, bases_vector const &bases_rng,
   vajoint_uint const n_fixef, std::vector<std::vector<int> > const &ders,
   bool const with_frailty):
  b(b_in.clone()), bases_rng(joint_bases::clone_bases(bases_rng)),
  n_fixef(n_fixef), ders_v(ders), with_frailty_v{with_frailty} {
    if(ders.size() != bases_rng.size())
      throw std::invalid_argument(
          "ders size does not match the number of basis functions (" +
            std::to_string(ders.size()) + ", " +
            std::to_string(bases_rng.size()) + ")");
  }

  vajoint_uint b_n_basis() const { return b_n_basis_v; }
  vajoint_uint rng_n_basis(size_t const idx) const {
    return rng_n_basis_v[idx];
  }
  vajoint_uint rng_n_weights(size_t const idx) const {
    return rng_n_weights_v[idx];
  }

  bool with_frailty() const { return with_frailty_v; }

  /**
   * evaluates the approximate expected cumulative hazard between
   * lower and upper times minus 1. The wk_mem and dwk_mem arguments are
   * for working memory. The last argument is possibly cached basis expansions.
   * Use a null pointer if there is no caching
   */
  template<class T>
  T operator()
    (node_weight const &nws, double const lower, double const upper,
     double const *design, double const * const fixef_design_varying,
     double const * const rng_design_varying, T const * fixef, T const * fixef_vary,
     T const * association, T const *VA_mean, T const * VA_vcov,
     T * wk_mem, double * dwk_mem, double const * cached_expansions) const {
    T out{0};
    bool const use_cache = cached_expansions;

    T * const association_M = wk_mem;
    wk_mem += n_basis_rng_p1;
    double * const dwk_mem_basis{dwk_mem + max_base_dim};

    vajoint_uint const n_vars
      {with_frailty() ? n_basis_rng_p1 : n_basis_rng_p1 - 1};

    for(vajoint_uint i = 0; i < nws.n_nodes; ++i){
      T fixef_term;

      if(use_cache){
        // compute the term from the time-varying fixed effects
        fixef_term =
          cfaad::dotProd(cached_expansions,
                         cached_expansions + b_n_basis(), fixef_vary);
        cached_expansions += b_n_basis();

        // construct the (association^T, 1).hat(M)(s) vector
        vajoint_uint idx{}, idx_association{};
        for(vajoint_uint j = 0; j < bases_rng.size(); ++j){
          for(vajoint_uint k = 0; k < rng_n_basis(j); ++k)
            association_M[idx + k] = 0;

          for(size_t l = 0; l < ders()[j].size(); ++l){
            for(vajoint_uint k = 0; k < rng_n_basis(j); ++k)
              association_M[idx + k] +=
                association[idx_association] * *cached_expansions++;
            ++idx_association;
          }
          idx += rng_n_basis(j);
        }
        association_M[idx] = 1;

      } else {
        // compute the term from the time-varying fixed effects
        double const node_val{scale_node_val(lower, upper, nws.ns[i])};
        (*b)(dwk_mem, dwk_mem_basis, node_val, fixef_design_varying);
        fixef_term =
          cfaad::dotProd(dwk_mem, dwk_mem + b_n_basis(), fixef_vary);

        // construct the (association^T, 1).hat(M)(s) vector
        vajoint_uint idx{}, idx_association{};
        double const * rng_design_varying_j{rng_design_varying};
        for(vajoint_uint j = 0; j < bases_rng.size(); ++j){
          for(vajoint_uint k = 0; k < rng_n_basis(j); ++k)
            association_M[idx + k] = 0;

          for(int der : ders()[j]){
            (*bases_rng[j])
              (dwk_mem, dwk_mem_basis, node_val, rng_design_varying_j, der);
            for(vajoint_uint k = 0; k < rng_n_basis(j); ++k)
              association_M[idx + k] +=
                association[idx_association] * dwk_mem[k];
            ++idx_association;
          }
          idx += rng_n_basis(j);
          rng_design_varying_j += rng_n_weights(j);
        }
        association_M[idx] = 1;
      }

      // add the mean log hazard term
      auto mean_term = cfaad::dotProd
        (association_M, association_M + n_vars, VA_mean);

      // add the quadratic term
      auto quad_term = cfaad::quadFormSym
        (VA_vcov, association_M, association_M + n_vars) / 2;

      // add the hazard term multiplied by the weight
      out += nws.ws[i] * exp(quad_term + mean_term + fixef_term);
    }

    // compute the fixed effect on the log hazard scale, add it, and return
    return (upper - lower) * out *
      exp(cfaad::dotProd(design, design + n_fixef, fixef));
  }

  /**
   * returns the needed working memory of eval. The first elements is the
   * required T memory and the second element is the required double memory.
   */
  std::array<size_t, 2> const& n_wmem() const {
    return n_wmem_v;
  }

  std::vector<std::vector<int> > const & ders() const {
    return ders_v;
  }

  /// sets the cached basis expansions
  void cache_expansions
    (double const lower, double const upper, double * cache_mem,
     double * wk_mem, node_weight const &nws,
     double const * const fixef_design_varying,
     double const * const rng_design_varying){
    // evaluates the cached expansions
    for(vajoint_uint i = 0; i < nws.n_nodes; ++i)
      cache_mem = cache_expansion_at
        (scale_node_val(lower, upper, nws.ns[i]), cache_mem, wk_mem,
         fixef_design_varying, rng_design_varying);
  }

  /***
   * sets the cached expansion at given node value. Returns a pointer at the end
   * of the used memory
   */
  double * cache_expansion_at
    (double const at, double * cache_mem, double * wk_mem,
     double const *fixef_design_varying, double const * rng_design_varying){
    (*b)(cache_mem, wk_mem, at, fixef_design_varying);
    cache_mem += b_n_basis();

    for(vajoint_uint base = 0; base < bases_rng.size(); ++base){
      for(int der : ders()[base]){
        (*bases_rng[base])(cache_mem, wk_mem, at, rng_design_varying, der);
        cache_mem += rng_n_basis(base);
      }
      rng_design_varying += rng_n_weights(base);
    }

    return cache_mem;
  }

  vajoint_uint cache_mem_per_node() const {
    vajoint_uint out{b_n_basis()};
    for(vajoint_uint i = 0; i < bases_rng.size(); ++i)
      out += ders()[i].size() * rng_n_basis(i);
    return out;
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
   * design matrices for the time-varying fixed effects (one for each type of
   * outcome)
   */
  std::vector<simple_mat<double> > fixef_design_varying_mats;
  /**
   * design matrices for the time-varying random effects (one for each type of
   * outcome)
   */
  std::vector<simple_mat<double> > rng_design_varying_mats;
  /**
   * functors to compute the approximation of the expected cumulative hazard
   * (one for each type of outcome)
   */
  std::vector<expected_cum_hazzard> cum_hazs;
  /// vector of vectors with lower bounds, upper bounds, and outcome indicators
  struct obs_info_obj {
    double lb, ub;
    bool event;
  };
  std::vector<std::vector<obs_info_obj> > obs_info;
  /// the largest basis dimension
  vajoint_uint max_basis_dim{};

  /// the indices of the parameters
  subset_params par_idx;
  /// the number of type of outcomes
  vajoint_uint n_outcomes_v = bases_fix.size();
  /// the required working memory
  std::array<size_t, 2> wmem_w;

  /// holds memory for the cached expansions
  std::vector<simple_mat<double> > cached_expansions;

  /// the cached quadrature nodes and weights
  std::vector<double> cached_nodes, cached_weights;

  bool has_cached_expansions() const {
    return cached_expansions.size() > 0;
  }

public:
  survival_dat() = default;

  survival_dat
    (bases_vector const &bases_fix_in, bases_vector const &bases_rng_in,
     std::vector<simple_mat<double> > &design_mats,
     std::vector<simple_mat<double> > &fixef_design_varying_mats,
     std::vector<simple_mat<double> > &rng_design_varying_mats,
     subset_params const &par_idx, std::vector<obs_input> const &input,
     std::vector<std::vector<std::vector<int> > > &ders):
    bases_fix{joint_bases::clone_bases(bases_fix_in)},
    bases_rng{joint_bases::clone_bases(bases_rng_in)},
    design_mats{design_mats},
    fixef_design_varying_mats{fixef_design_varying_mats},
    rng_design_varying_mats{rng_design_varying_mats},
    par_idx{par_idx}
  {
    if(par_idx.surv_info().size() != n_outcomes_v)
      throw std::invalid_argument("surv_info().size() != n_outcomes");
    if(design_mats.size() != n_outcomes_v)
      throw std::invalid_argument("design_mats.size() != n_outcomes");
    if(fixef_design_varying_mats.size() != n_outcomes_v)
      throw std::invalid_argument("fixef_design_varying_mats.size() != n_outcomes");
    if(rng_design_varying_mats.size() != n_outcomes_v)
      throw std::invalid_argument("rng_design_varying_mats.size() != n_outcomes");
    if(bases_fix.size() != n_outcomes_v)
      throw std::invalid_argument("bases_fix.size() != n_outcomes");
    if(input.size() != n_outcomes_v)
      throw std::invalid_argument("input.size() != n_outcomes");
    if(bases_rng.size() != par_idx.marker_info().size())
      throw std::invalid_argument("bases_rng.size() != marker_info().size()");
    if(ders.size() != n_outcomes_v)
      throw std::invalid_argument(
          "ders size does not match the number of survival types (" +
            std::to_string(ders.size()) + ", " +
            std::to_string(n_outcomes_v) + ")");

    // create the cum_hazs
    cum_hazs.reserve(n_outcomes_v);
    for(vajoint_uint i = 0; i < n_outcomes_v; ++i)
      cum_hazs.emplace_back
        (*bases_fix[i], bases_rng, design_mats[i].n_rows(), ders[i],
         par_idx.surv_info()[i].with_frailty);

    // set the required working memory
    wmem_w = {0, 0};
    for(auto &ch : cum_hazs){
      auto ch_mem = ch.n_wmem();
      wmem_w[0] = std::max(wmem_w[0], ch_mem[0]);
      wmem_w[1] = std::max(wmem_w[1], ch_mem[1]);
    }

    vajoint_uint const n_shared_p1{par_idx.n_shared() + 1};
    wmem_w[0] += n_shared_p1 * (n_shared_p1 + 1);

    max_basis_dim = n_shared_p1;
    for(auto &b : bases_fix)
      max_basis_dim = std::max<vajoint_uint>(b->n_basis(), max_basis_dim);
    wmem_w[1] += max_basis_dim;

    // add the observations
    obs_info.resize(n_outcomes_v);
    for(vajoint_uint type = 0; type < n_outcomes_v; ++type){
      // check the par_info matches
      auto &surv_info_type = par_idx.surv_info()[type];
      if(design_mats[type].n_rows() != surv_info_type.n_fix)
        throw std::invalid_argument("esign_mats[type].n_rows() != surv_info_type.n_fix");
      if(bases_fix[type]->n_basis() != surv_info_type.n_variying)
        throw std::invalid_argument("bases_fix[type]->n_basis() != surv_info_type.n_variying");
      if(ders[type].size() != par_idx.marker_info().size())
        throw std::invalid_argument("ders[type].size() != par_idx.marker_info().size()");
      for(unsigned j = 0; j < par_idx.marker_info().size(); ++j)
        if(ders[type][j].size() != surv_info_type.n_associations[j])
          throw std::invalid_argument("ders[type][j].size() != surv_info_type.n_associations[j]");

      if(fixef_design_varying_mats[type].n_rows() != bases_fix[type]->n_weights())
        throw std::invalid_argument("fixef_design_varying_mats[type].n_rows() != bases_fix[type]->n_weights()");

      unsigned n_weights_rng{};
      for(unsigned j = 0; j < par_idx.marker_info().size(); ++j)
        n_weights_rng += cum_hazs[type].rng_n_weights(j);
      if(rng_design_varying_mats[type].n_rows() != n_weights_rng)
        throw std::invalid_argument("rng_design_varying_mats[type].n_rows() != n_weights_rng");

      // add the outcomes
      obs_info[type].reserve(input[type].n_obs);
      if(input[type].n_obs != design_mats[type].n_cols())
        throw std::invalid_argument
          ("input[type].n_obs != design_mats[type].n_cols()");

      for(vajoint_uint obs = 0; obs < input[type].n_obs; ++obs)
        obs_info[type].push_back
          (obs_info_obj
            {input[type].lbs[obs], input[type].ubs[obs], input[type].event[obs] == 1});
    }
  }

  void set_cached_expansions(node_weight const &nws){
    if(has_cached_expansions()){
      /// check if we already use the same quadrature rule
      bool is_same_rule
        {nws.n_nodes == cached_nodes.size() &&
          nws.n_nodes == cached_weights.size()};
      for(vajoint_uint i = 0; i < nws.n_nodes && is_same_rule; ++i)
        is_same_rule &= nws.ns[i] == cached_nodes[i] &&
          nws.ws[i] == cached_weights[i];

      if(is_same_rule)
        return;
    }

    vajoint_uint const n_nodes{nws.n_nodes};
    cached_weights.resize(n_nodes);
    std::copy(nws.ws, nws.ws + n_nodes, cached_weights.begin());
    cached_nodes.resize(n_nodes);
    std::copy(nws.ns, nws.ns + n_nodes, cached_nodes.begin());

    cached_expansions.clear();
    cached_expansions.reserve(obs_info.size());
    for(size_t type = 0; type < obs_info.size(); ++type){
      auto &info_objs = obs_info[type];
      auto &haz_type = cum_hazs[type];
      auto const &fixef_design_varying_mat = fixef_design_varying_mats[type];
      auto const &rng_design_varying_mat = rng_design_varying_mats[type];

      size_t const n_basis_cols{(n_nodes + 1) * info_objs.size()};
      std::vector<double> wk_mem(haz_type.n_wmem()[1]);

      cached_expansions.emplace_back
        (haz_type.cache_mem_per_node(), n_basis_cols);
      for(size_t obs = 0; obs < info_objs.size(); ++obs){
        double * cache_mem{cached_expansions.back().col(obs * (n_nodes + 1))};

        // store the event time as the first column
        if(info_objs[obs].event)
          cache_mem = haz_type.cache_expansion_at
            (info_objs[obs].ub, cache_mem, wk_mem.data(),
             fixef_design_varying_mat.col(obs),
             rng_design_varying_mat.col(obs));

        // use the rest of the columns for the terms from the cumulative hazard
        haz_type.cache_expansions
          (info_objs[obs].lb, info_objs[obs].ub, cache_mem, wk_mem.data(), nws,
           fixef_design_varying_mat.col(obs),
           rng_design_varying_mat.col(obs));
      }
    }
  }

  /// clears the cached expansions
  void clear_cached_expansions(){
    cached_expansions.clear();
    cached_expansions.shrink_to_fit();

    cached_nodes.clear();
    cached_nodes.shrink_to_fit();

    cached_weights.clear();
    cached_weights.shrink_to_fit();
  }

  /// returns the number of lower bound terms of a given type
  vajoint_uint n_terms(const vajoint_uint type) const {
    return obs_info[type].size();
  }

  /// returns the number of types of survival processes
  vajoint_uint n_outcomes() const {
    return n_outcomes_v;
  }

  /// evaluates the lower bound of observation idx for the type of outcome
  template<class T>
  T operator()
    (T const *param, T *wk_mem, const vajoint_uint idx, const vajoint_uint type,
     double * dwk_mem, node_weight nws) const {
    // get the information for the outcome and event type
    obs_info_obj const &info{obs_info[type][idx]};
    expected_cum_hazzard const &haz{cum_hazs[type]};
    auto const &surv_info{par_idx.surv_info()[type]};

    double const * cached_expansions_pass{nullptr};
    bool const use_cache{has_cached_expansions()};
    if(use_cache){
      nws = { cached_nodes.data(), cached_weights.data(),
              static_cast<vajoint_uint>(cached_nodes.size()) };

      cached_expansions_pass = cached_expansions[type].col
        ((nws.n_nodes + 1) * idx);
    }

    // compute the approximate expected log hazard if needed
    T out{0};
    double const * const design{design_mats[type].col(idx)},
                 * const fixef_design_varying
                  {fixef_design_varying_mats[type].col(idx)},
                 * const rng_design_varying
                 {rng_design_varying_mats[type].col(idx)};
    vajoint_uint const n_shared{par_idx.n_shared()},
                    n_shared_p1{n_shared + 1};
    if(info.event){
      out -= cfaad::dotProd(design, design + surv_info.n_fix,
                            param + surv_info.idx_fix);

      double * const basis_wmem{dwk_mem + max_basis_dim};

      if(use_cache){
        out -= cfaad::dotProd
          (cached_expansions_pass, cached_expansions_pass + haz.b_n_basis(),
           param + surv_info.idx_varying);
        cached_expansions_pass += haz.b_n_basis();

        vajoint_uint offset{}, idx_association{surv_info.idx_association};
        for(size_t i = 0; i < bases_rng.size(); ++i){
          for(size_t j = 0; j  < haz.ders()[i].size(); ++j){
            auto M_VA_mean = cfaad::dotProd
              (cached_expansions_pass,
               cached_expansions_pass + haz.rng_n_basis(i),
               param + par_idx.va_mean() + offset);
            cached_expansions_pass += haz.rng_n_basis(i);
            out -= param[idx_association++] * M_VA_mean;
          }

          offset += haz.rng_n_basis(i);

        }

      } else {
        (*bases_fix[type])(dwk_mem, basis_wmem, info.ub, fixef_design_varying);
        out -= cfaad::dotProd(dwk_mem, dwk_mem + haz.b_n_basis(),
                              param + surv_info.idx_varying);

        double const * rng_design_varying_j{rng_design_varying};
        vajoint_uint offset{}, idx_association{surv_info.idx_association};
        for(size_t i = 0; i < bases_rng.size(); ++i){
          for(int der : haz.ders()[i]){
            (*bases_rng[i])
              (dwk_mem, basis_wmem, info.ub, rng_design_varying_j, der);
            auto M_VA_mean = cfaad::dotProd
              (dwk_mem, dwk_mem + haz.rng_n_basis(i),
               param + par_idx.va_mean() + offset);
            out -= param[idx_association++] * M_VA_mean;
          }

          offset += haz.rng_n_basis(i);
          rng_design_varying_j += haz.rng_n_weights(i);

        }
      }

      // the frailty term
      if(haz.with_frailty())
        out -=
          param[par_idx.va_mean() + n_shared + par_idx.frailty_offset(type)];
    }

    // add the term from the approximate expected cumulative hazard
    T * const VA_mean{wk_mem};
    for(vajoint_uint i = 0; i < n_shared; ++i)
      VA_mean[i] = param[par_idx.va_mean() + i];

    if(haz.with_frailty())
      VA_mean[n_shared] =
        param[par_idx.va_mean() + n_shared + par_idx.frailty_offset(type)];
    wk_mem += n_shared_p1;

    T * const VA_vcov{wk_mem};
    vajoint_uint const rng_dim{n_shared + par_idx.n_shared_surv()};

    if(haz.with_frailty()){
      T const * const vcov_full{param + par_idx.va_vcov()};
      for(vajoint_uint j = 0; j < n_shared; ++j){
        for(vajoint_uint i = 0; i < n_shared; ++i)
          VA_vcov[i + j * n_shared_p1] = vcov_full[i + j * rng_dim];

        VA_vcov[n_shared + j * n_shared_p1] =
          vcov_full[n_shared + par_idx.frailty_offset(type) + j * rng_dim];
      }

      vajoint_uint const offset{n_shared + par_idx.frailty_offset(type)};
      for(vajoint_uint i = 0; i < n_shared; ++i)
        VA_vcov[i + n_shared * n_shared_p1] = vcov_full[i + offset * rng_dim];

      VA_vcov[n_shared + n_shared * n_shared_p1] =
        vcov_full[offset + offset * rng_dim];

    } else {
      T const * const vcov_full{param + par_idx.va_vcov()};
      for(vajoint_uint j = 0; j < n_shared; ++j)
        for(vajoint_uint i = 0; i < n_shared; ++i)
          VA_vcov[i + j * n_shared] = vcov_full[i + j * rng_dim];

    }
    wk_mem += n_shared_p1 * n_shared_p1;

    out += haz
      (nws, info.lb, info.ub, design, fixef_design_varying, rng_design_varying,
       param + surv_info.idx_fix, param + surv_info.idx_varying,
       param + surv_info.idx_association, VA_mean, VA_vcov, wk_mem, dwk_mem,
       cached_expansions_pass);

    return out;
  }

  /**
   * returns the needed working memory of eval. The first elements is the
   * required T memory and the second element is the required double memory.
   */
  std::array<size_t, 2> const & n_wmem() const {
    return wmem_w;
  }
};

} // namespace survival

#endif
