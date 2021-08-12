#ifndef MARKER_TERM_H
#define MARKER_TERM_H

#include "VA-joint-config.h"
#include "VA-parameter.h"
#include "bases.h"
#include <iterator>
#include <type_traits>
#include <memory>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <vector>
#include <unordered_map>
#include <cfaad/AAD.h>
#include "simple-mat.h"
#include <limits>

namespace marker {
/// used for pre-computed data for eval
struct comp_dat {
  /// the indices of the non-missing outcomes
  const std::vector<vajoint_uint> indices;
  /// the factorization of the covariance matrix
  const cfaad::CholFactorization vcov_factorization;
  /// dimension of the random effect
  const vajoint_uint n_rngs;

  /**
   * current parameters, working memory, the indices of the model, and
   * a missingness flag
   */
  comp_dat(double const *param, double *wk_mem, subset_params const &par_idx,
           unsigned const missingness_flag);
};

/// holds data for marker observations
class marker_dat {
public:
  /// the indices of the parameters
  subset_params const par_idx;
  /// maximum number of possible markers per observation
  static constexpr unsigned max_markers{31};
  /// the maximum number of markers per observations
  const vajoint_uint n_markers
    { static_cast<vajoint_uint>(par_idx.get_marker_info().size()) };
  /// the number of observations
  const vajoint_uint n_obs;

private:
  /// the bases for the time-varying fixed effects
  joint_bases::bases_vector bases_fix;
  /// the bases for the time-varying random effects
  joint_bases::bases_vector bases_rng;
  /// the number of fixed effects
  const vajoint_uint n_fixed_effects
    {
      std::accumulate(
        par_idx.get_marker_info().begin(), par_idx.get_marker_info().end(),
        vajoint_uint{},
        [](const vajoint_uint x, subset_params::marker const &r){
          return x + r.n_fix;
        })
    };
  /// the number of time-varying fixed effect basis function
  const vajoint_uint n_basis_fix
    { std::accumulate(
        bases_fix.begin(), bases_fix.end(), vajoint_uint{},
        [](vajoint_uint x, joint_bases::bases_vector::value_type const &r){
          return x + r->get_n_basis();
        })
    };
  /// the number of time-varying random effect basis function
  const vajoint_uint n_basis_rng
    { std::accumulate(
        bases_rng.begin(), bases_rng.end(), vajoint_uint{},
        [](vajoint_uint x, joint_bases::bases_vector::value_type const &r){
          return x + r->get_n_basis();
        })
    };
  /// the first dimension of the combined design matrix
  const vajoint_uint dim_design
    {n_fixed_effects + n_basis_fix + n_basis_rng};

  /**
   * contains the stacked columns of the fixed effect design matrix,
   * time-varying fixed design matrix, and time-varying random effect
   * design matrix in that order.
   */
  simple_mat<double> design_mats{dim_design, n_obs};
  /// contains the matrix with the outcomes
  simple_mat<double> outcomes{n_markers, n_obs};

  /// bit flags for missing variables
  std::vector<unsigned> missingness{std::vector<unsigned>(n_obs)};

  /// the offsets between the time-varying random effects for each marker
  std::vector<vajoint_uint> offsets_rng;

  /// the needed working memory
  size_t n_wmem{
    n_markers * (n_markers + 1) + n_basis_rng * n_basis_rng +
    par_idx.get_n_shared() * (par_idx.get_n_shared() + n_markers)};

  /**
   * the pre-computed data for eval. We do not want to re-compute the matrix
   * factorization of the covariance matrix for each outcome. Thus, we do it
   * once for all observed missingness patterns. See setup */
  std::unordered_map<unsigned, comp_dat> pre_comp_dat;

public:
  /// creates the object and allocates the memory that is needed
  marker_dat
  (subset_params const &par_idx, const vajoint_uint n_obs,
   joint_bases::bases_vector const &bases_fix,
   joint_bases::bases_vector const &bases_rng):
  par_idx{par_idx}, n_obs{n_obs},
  bases_fix{joint_bases::clone_bases(bases_fix)},
  bases_rng{joint_bases::clone_bases(bases_rng)}
  {
    if(n_markers != bases_fix.size())
      throw std::runtime_error
        ("number of markers and time-varying fixed effect bases differs");
    if(n_markers != bases_rng.size())
      throw std::runtime_error
        ("number of markers and time-varying random effect bases differs");

    if(n_markers > max_markers)
      throw std::runtime_error("too many markers");

    // set all outcomes to be missing
    unsigned def_miss{};
    for(unsigned i = 0; i < n_markers; ++i)
      def_miss |= (1u << i);
    std::fill(missingness.begin(), missingness.end(), def_miss);

    // sets the indices of the time-varying random effects
    offsets_rng.reserve(n_markers);
    vajoint_uint next_idx{};
    for(vajoint_uint i = 0; i < n_markers; ++i){
      offsets_rng.emplace_back(next_idx);
      next_idx += bases_rng[i]->get_n_basis();
    }
  }

  /**
   * sets design matrices for the varying effects using an iterator over the
   * observation times.
   */
  template<class I>
  void set_design_mats(I obs_time){
    static_assert(std::is_same<typename std::iterator_traits<I>::value_type,
                               double>::value, "iterator is not to doubles");

    for(vajoint_uint i = 0; i < n_obs; ++i, ++obs_time){
      double *mem = design_mats.col(i) + n_fixed_effects;
      for(auto &x : bases_fix){
        x->operator()(mem, *obs_time);
        mem += x->get_n_basis();
      }
      for(auto &x : bases_rng){
        x->operator()(mem, *obs_time);
        mem += x->get_n_basis();
      }
    }
  }

  /**
   * sets the outcome for a given index for a given individual. The outcome
   * of this individual is then set as observed */
  void set_outcome(const vajoint_uint idx, const unsigned obs_type,
                   const double value){
    outcomes.col(idx)[obs_type] = value;
    // mark as observed
    missingness[idx] &= ~(1u << obs_type);
  }

  /// sets the fixed effect design matrix of a given type of outcome
  template<class I>
  void set_fixef_design(const vajoint_uint idx, const vajoint_uint obs_type,
                        I values){
    static_assert(std::is_same<typename std::iterator_traits<I>::value_type,
                               double>::value, "iterator is not to doubles");

    std::copy(values, values + par_idx.get_marker_info()[obs_type].n_fix,
              design_mats.col(idx) + par_idx.get_fixef_idx_marker(obs_type));
  }

  /**
   * Sets up objects to evaluate the expected log conditional density. Has to be
   * called prior to calling eval.
   */
  void setup(double const *param, double *wk_mem){
    pre_comp_dat.clear();
    for(unsigned missingness_flag : missingness)
      if(pre_comp_dat.find(missingness_flag) == pre_comp_dat.end())
        pre_comp_dat.emplace(
          missingness_flag, comp_dat{param, wk_mem, par_idx, missingness_flag});
  }

  /// computes the expected log conditional density of a given observation
  template<class T>
  T eval(T const *param, T *wk_mem, const vajoint_uint idx){
    unsigned const missingness_flag = missingness[idx];
    comp_dat const &c_dat = pre_comp_dat.at(missingness_flag);
    const std::vector<vajoint_uint> &indices = c_dat.indices;
    const cfaad::CholFactorization &vcov_factorization =
      c_dat.vcov_factorization;
    vajoint_uint const n_indices = indices.size();

    // TODO: could handle can optimize for the scenarios where one or more
    //       markers are observed

    // construct the copy of the covariance matrix
    T * const Sig = wk_mem;
    {
      T const * const s{param + par_idx.get_idx_error_term()};
      for(vajoint_uint j = 0; j < n_indices; ++j)
        for(vajoint_uint i = 0; i < n_indices; ++i)
          wk_mem[i + j * n_indices] = s[indices[i] + n_markers * indices[j]];
    }
    wk_mem += n_indices * n_indices;

    // construct the copy of random effect covariance matrix
    T * const Psi = wk_mem;
    {
      T const * const s{param + par_idx.get_idx_va_vcov()};
      vajoint_uint const n_rngs{n_basis_rng + par_idx.get_n_shared_surv()};
      vajoint_uint j_idx{};
      for(vajoint_uint j = 0; j < n_indices; ++j){
        vajoint_uint const offset_j{offsets_rng[j]},
                        n_effects_j
                        {par_idx.get_marker_info()[indices[j]].n_rng};
        for(vajoint_uint jj = 0; jj < n_effects_j; ++jj, ++j_idx){
          vajoint_uint i_idx{};
          for(vajoint_uint i = 0; i < n_indices; ++i){
            vajoint_uint const offset_i{offsets_rng[i]},
                            n_effects_i
                            {par_idx.get_marker_info()[indices[i]].n_rng};
            for(vajoint_uint ii = 0; ii < n_effects_i; ++ii, ++i_idx)
              Psi[i_idx + j_idx * c_dat.n_rngs] =
                s[i_idx + offset_i + (j_idx + offset_j) * n_rngs];
          }
        }
      }
    }
    wk_mem += c_dat.n_rngs * c_dat.n_rngs;

    // compute the difference between the observed and expected outcome
    T out{0}; // the output
    double const * const design{design_mats.col(idx)};
    {
      T * const delta = wk_mem;
      for(vajoint_uint i = 0; i < n_indices; ++i){
        delta[i] = outcomes.col(idx)[indices[i]];

        { // fixed effect
          vajoint_uint const offset{par_idx.get_fixef_idx_marker(indices[i])};
          double const * d{design + offset};
          delta[i] -= cfaad::dotProd
            (d, d + par_idx.get_marker_info()[indices[i]].n_fix,
             param + offset);
        }
        { // time-varying fixed effect
          vajoint_uint const offset{par_idx.get_varying_idx_marker(indices[i])};
          double const * d{design + offset};
          delta[i] -= cfaad::dotProd
            (d, d + par_idx.get_marker_info()[indices[i]].n_variying,
             param + offset);
        }
        { // time-varying random effects
          unsigned const offset{offsets_rng[indices[i]]};
          double const * d{design + offset + n_fixed_effects + n_basis_fix};
          delta[i] -= cfaad::dotProd
            (d, d + par_idx.get_marker_info()[indices[i]].n_rng,
             param + par_idx.get_idx_va_mean() + offset);
        }
      }
      out += cfaad::quadFormInv(delta, Sig, vcov_factorization);
    }

    // add the determinant term and the constant
    constexpr double log_2_pi{1.83787706640935};
    out += static_cast<double>(n_indices) * log_2_pi;
    out += cfaad::logDeter(Sig, vcov_factorization);

    // add the trace term. Start by computing Psi.M
    T * const Psi_M = wk_mem;
    wk_mem += c_dat.n_rngs * n_indices;
    for(vajoint_uint i = 0; i < c_dat.n_rngs * n_indices; ++i)
      Psi_M[i] = 0;

    for(vajoint_uint i = 0; i < n_indices; ++i){
      unsigned const offset{offsets_rng[indices[i]]},
                    n_rng_i{par_idx.get_marker_info()[indices[i]].n_rng};
      double const * d{design + offset + n_fixed_effects + n_basis_fix};
      cfaad::matVecProd
        (Psi + offset * c_dat.n_rngs, Psi + (offset + n_rng_i) * c_dat.n_rngs,
         d, d + n_rng_i, Psi_M + offset * c_dat.n_rngs, false);
    }

    // compute M^T.Psi.M
    T * const MT_Psi_M = wk_mem;
    for(vajoint_uint i = 0; i < n_indices * n_indices; ++i)
      MT_Psi_M[i] = 0;

    for(vajoint_uint i = 0; i < n_indices; ++i){
      unsigned const offset{offsets_rng[indices[i]]},
                    n_rng_i{par_idx.get_marker_info()[indices[i]].n_rng};
      double const * d{design + offset + n_fixed_effects + n_basis_fix};
      // TODO: the code is wrong. We need a leading dimension
      cfaad::matVecProd
        (Psi_M + offset * c_dat.n_rngs,
         Psi_M + (offset + i + 1) * c_dat.n_rngs, d, d + n_rng_i,
         MT_Psi_M + offset * n_indices, true);
    }

    out += cfaad::trInvMatMat(Sig, MT_Psi_M, vcov_factorization);

    return out / 2;
  }

  size_t get_n_wmem() const {
    return n_wmem;
  }
};

/// this is a helper class to setup a marker_dat object
struct setup_marker_dat_helper {
  /// the fixed effect design matrix
  simple_mat<double> fixef_design;
  /// pointer to id of each observation (sorted)
  int const * ids;
  /// pointer to the observation time of each observation (sorted within id)
  double const *obs_time;
  /// pointer to the observed outcomes
  double const *obs;

  /// returns the number of observations
  vajoint_uint n_obs() const {
    return fixef_design.n_cols();
  }

  /// constructs the helper class and validates the data
  setup_marker_dat_helper
    (double * fixef, vajoint_uint const n_fixef, vajoint_uint const arg_n_obs,
     int const *ids, double const *obs_time, double const *obs):
    fixef_design{fixef, n_fixef, arg_n_obs}, ids{ids}, obs_time{obs_time},
    obs{obs} {
      // checks that the ids and observations times are sorted
      int const * ids_i{ids};
      double const * obs_time_i{obs_time};
      for(vajoint_uint i = 1; i < n_obs(); ++i, ++ids_i, ++obs_time_i){
        if(ids_i[0] > ids_i[1])
          throw std::invalid_argument("id > next id" );
        if(ids_i[0] == ids_i[1] && obs_time_i[0] >= obs_time_i[1])
          throw std::invalid_argument("obs_time >= next obs_time");
      }
    }
};

/// creates a comp_dat object after checking the input arguments
marker_dat get_comp_dat
  (std::vector<setup_marker_dat_helper> &input_dat,
   subset_params const &par_idx,
   joint_bases::bases_vector const &bases_fix,
   joint_bases::bases_vector const &bases_rng);

} // namespace marker

#endif
