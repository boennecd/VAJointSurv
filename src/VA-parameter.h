#ifndef VA_PARAMETER_H
#define VA_PARAMETER_H
#include "VA-joint-config.h"
#include <vector>
#include <string>
#include <stdexcept>

/**
 * handles finding parameters vector elements. Zero-based indices are used for
 * all class members. The template argument show whether it is with the
 * parameterization where the covariance matrices are triangular matrices. */
class subset_params {
  static constexpr bool is_traingular_default = false;

public:
  struct marker {
    /**
     * number of fixed effects, varying fixed effects, and number of random
     * effects
     */
    vajoint_uint n_fix, n_variying, n_rng;
    vajoint_uint idx_fix = 0,
                 idx_varying = 0;

    marker(vajoint_uint const n_fix, vajoint_uint const n_variying,
           vajoint_uint const n_rng):
      n_fix(n_fix), n_variying(n_variying), n_rng(n_rng) { }
  };
  struct surv {
    /// number of fixed effects and varying fixed effects
    vajoint_uint n_fix, n_variying;
    /// number of associations parameters for each marker
    std::vector<unsigned> n_associations;
    /// is there a frailty term
    bool with_frailty;
    vajoint_uint idx_fix = 0,
                 idx_varying = 0,
                 idx_association = 0,
                 frailty_offset = 0;

    surv(vajoint_uint const n_fix, vajoint_uint const n_variying,
         const std::vector<unsigned> &n_associations, bool const with_frailty):
      n_fix(n_fix), n_variying(n_variying), n_associations{n_associations},
      with_frailty{with_frailty} { }
  };

private:
  // contains data for each type of outcome
  std::vector<marker> marker_info_v;
  std::vector<surv> surv_info_v;

  vajoint_uint idx_error_term = 0,
               idx_shared_effect = 0,
               idx_shared_surv = 0,
               idx_va_mean = 0,
               idx_va_vcov = 0,

               n_params_v = 0,
               n_params_w_va_v = 0,
               n_shared_effect = 0,
               // using the triangular parameterization
               idx_error_term_triangular = 0,
               idx_shared_effect_triangular = 0,
               idx_shared_surv_triangular = 0,
               idx_va_mean_triangular = 0,
               idx_va_vcov_triangular = 0,

               n_params_triangular_v = 0,
               n_parms_w_va_triangular_v = 0,
               n_shared_surv_v = 0;

  void re_compute_indices(){
    /// fill in the indices from the markers and the shared effect
    vajoint_uint idx{0L};
    n_shared_effect = 0;
    for(auto &info : marker_info_v){
      info.idx_fix = idx;
      idx += info.n_fix;
      n_shared_effect += info.n_rng;
    }

    for(auto &info : marker_info_v){
      info.idx_varying = idx;
      idx += info.n_variying;
    }

    /// fill in the indices from the survival outcomes
    n_shared_surv_v = 0;
    for(auto &info : surv_info_v){
      info.idx_fix = idx;
      idx += info.n_fix;
      info.idx_varying = idx;
      idx += info.n_variying;
      info.idx_association = idx;
      for(unsigned n : info.n_associations)
        idx += n;
      info.frailty_offset = n_shared_surv_v;
      if(info.with_frailty)
        ++n_shared_surv_v;
    }

    /// fill in the indices for the covariance matrices
    ([&,idx]() mutable {
      idx_error_term = idx;
      idx += static_cast<vajoint_uint>(
        marker_info_v.size() * marker_info_v.size()); // -Wconversion
      idx_shared_effect = idx;
      idx += n_shared_effect * n_shared_effect;
      idx_shared_surv = idx;
      idx += n_shared_surv_v * n_shared_surv_v;
      n_params_v = idx;

      idx_va_mean = idx;
      vajoint_uint const va_dim{n_shared_effect + n_shared_surv_v};
      idx += va_dim;
      idx_va_vcov = idx;
      n_params_w_va_v = idx + va_dim * va_dim;
    })();

    idx_error_term_triangular = idx;
    idx += dim_tri(marker_info_v.size());
    idx_shared_effect_triangular = idx;
    idx += dim_tri(n_shared_effect);
    idx_shared_surv_triangular = idx;
    idx += dim_tri(n_shared_surv_v);
    n_params_triangular_v = idx;

    idx_va_mean_triangular = idx;
    vajoint_uint const va_dim{n_shared_effect + n_shared_surv_v};
    idx += va_dim;
    idx_va_vcov_triangular = idx;
    n_parms_w_va_triangular_v = idx + dim_tri(va_dim);
  }

public:
  std::vector<marker> const & marker_info() const {
    return marker_info_v;
  }

  std::vector<surv> const & surv_info() const {
    return surv_info_v;
  }

  /// adds a marker to the model
  void add_marker(marker const &info){
    if(surv_info_v.size() > 0)
      throw std::runtime_error("marker added after survival terms");
    marker_info_v.push_back(info);
    re_compute_indices();
  }

  /**
   * adds a survival outcome to the model. This must be done after the all the
   * markers has been added because of the number of association parameter.
   */
  void add_surv(surv const &info){
    if(info.n_associations.size() != marker_info().size())
      throw std::invalid_argument
        ("new surv objects nubmer of associations parameters do not match the number of markers");
    surv_info_v.push_back(info);
    re_compute_indices();
  }

  /// returns the index of the fixed effect coefficients for a given marker
  template<bool is_traingular = is_traingular_default>
  vajoint_uint fixef_marker(vajoint_uint const idx) const {
    return marker_info_v[idx].idx_fix;
  }

  /**
   * returns the index of the fixed varying effect's coefficients for a given
   * marker
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint fixef_vary_marker(vajoint_uint const idx) const {
    return marker_info_v[idx].idx_varying;
  }

  /**
   * the covariance matrices are last. Thus, it is useful to know where they
   * start which this function returns
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint vcov_start() const {
    return vcov_marker<is_traingular>();
  }

  /// end of the covariance matrices
  template<bool is_traingular = is_traingular_default>
  vajoint_uint vcov_end() const {
    return va_mean<is_traingular>();
  }

  /// returns the index of the covariance matrix for the markers' residual
  template<bool is_traingular = is_traingular_default>
  vajoint_uint vcov_marker() const {
    return is_traingular ? idx_error_term_triangular : idx_error_term;
  }

  /// returns the index of the covariance matrix for the shared effect
  template<bool is_traingular = is_traingular_default>
  vajoint_uint vcov_vary() const{
    return is_traingular ? idx_shared_effect_triangular : idx_shared_effect;
  }

  /**
   * returns the index of the fixed effect coefficients for a given survival
   * outcome
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint fixef_surv(vajoint_uint const idx) const {
    return surv_info_v[idx].idx_fix;
  }

  /**
   * returns the index of the fixed varying effect's coefficients for a given
   * survival outcome
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint fixef_vary_surv(vajoint_uint const idx) const {
    return surv_info_v[idx].idx_varying;
  }

  /**
   * returns the index of the association parameter for a given survival outcome
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint association(vajoint_uint const idx)
  const {
    return surv_info_v[idx].idx_association;
  }

  /// returns the offset for the frailty variable for the survival type
  vajoint_uint frailty_offset(vajoint_uint const idx) const {
    return surv_info_v[idx].frailty_offset;
  }

  /// returns the index of the covarinace matrix for the frailties
  template<bool is_traingular = is_traingular_default>
  vajoint_uint vcov_surv() const {
    return is_traingular ? idx_shared_surv_triangular : idx_shared_surv;
  }

  /// returns the number of parameters
  template<bool is_traingular = is_traingular_default>
  vajoint_uint n_params() const {
    return is_traingular ? n_params_triangular_v : n_params_v;
  }

  /// returns the location of the mean of the VA distribution
  template<bool is_traingular = is_traingular_default>
  vajoint_uint va_mean() const {
    return is_traingular ? idx_va_mean_triangular : idx_va_mean;
  }

  /// returns the end of the mean of the VA distribution
  template<bool is_traingular = is_traingular_default>
  vajoint_uint va_mean_end() const {
    return va_vcov<is_traingular>();
  }

  /// returns the location of the covariance of the VA distribution
  template<bool is_traingular = is_traingular_default>
  vajoint_uint va_vcov() const {
    return is_traingular ? idx_va_vcov_triangular : idx_va_vcov;
  }

  /// returns the end of va_vcov
  template<bool is_traingular = is_traingular_default>
  vajoint_uint va_vcov_end() const {
    return n_params_w_va<is_traingular>();
  }

  /**
   * returns the start of the variational parameters. The variational parameters
   * lie in continuous memory.
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint va_par_start() const {
    return va_mean<is_traingular>();
  }

  /**
   * returns the end of the variational parameters. The variational parameters
   * lie in continuous memory.
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint va_par_end() const {
    return va_vcov_end<is_traingular>();
  }

  /// returns the number of shared random effects
  vajoint_uint n_shared() const {
    return n_shared_effect;
  }

  /// returns the number of shared random effects for the survival outcomes
  vajoint_uint n_shared_surv() const {
    return n_shared_surv_v;
  }

  /**
   * returns the number parameters including the variational parameters for one
   * cluster.
   */
  template<bool is_traingular = is_traingular_default>
  vajoint_uint n_params_w_va() const {
    return is_traingular ? n_parms_w_va_triangular_v : n_params_w_va_v;
  }

  /// returns the number of VA parameters
  template<bool is_traingular = is_traingular_default>
  vajoint_uint n_va_params() const {
    return n_params_w_va<is_traingular>() - n_params<is_traingular>();
  }

  /// gives the names of the parameters
  std::vector<std::string> param_names
    (bool const is_traingular = is_traingular_default) const;

  /// gives the names of the variational parameters
  std::vector<std::string> va_param_names
    (bool const is_traingular = is_traingular_default) const;
};

#endif
