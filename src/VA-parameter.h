#ifndef VA_PARAMETER_H
#define VA_PARAMETER_H
#include "VA-joint-config.h"
#include "vector"

/// handles finding parameters vector elements. Zero-based indices are used for
/// all class members
class subset_params {
public:
  struct marker {
    vajoint_uint n_fix, n_variying, n_rng;
    vajoint_uint idx_fix = 0,
                 idx_varying = 0;

    marker(vajoint_uint const n_fix, vajoint_uint const n_variying,
           vajoint_uint const n_rng):
      n_fix(n_fix), n_variying(n_variying), n_rng(n_rng) { }
  };
  struct survival {
    vajoint_uint n_fix, n_variying;
    vajoint_uint idx_fix = 0,
                 idx_varying = 0,
                 idx_association = 0;

    survival(vajoint_uint const n_fix, vajoint_uint const n_variying):
      n_fix(n_fix), n_variying(n_variying) { }
  };

private:
  // contains data for each type of outcome
  std::vector<marker> marker_info;
  std::vector<survival> survival_info;

  vajoint_uint idx_error_term = 0,
              idx_shared_effect = 0,
              idx_shared_survival = 0,
              n_params = 0,
              n_shared_effect = 0;

  inline void re_compute_indices(){
    /// fill in the indices from the markers and the shared effect
    vajoint_uint idx { 0L };
    n_shared_effect = 0;
    for(auto &info : marker_info){
      info.idx_fix = idx;
      idx += info.n_fix;
      info.idx_varying = idx;
      idx += info.n_variying;

      n_shared_effect += info.n_rng;
    }

    idx_error_term = idx;
    idx += marker_info.size() * marker_info.size();
    idx_shared_effect = idx;
    idx += n_shared_effect * n_shared_effect;

    /// fill in the indices from the survival outcomes
    for(auto &info : survival_info){
      info.idx_fix = idx;
      idx += info.n_fix;
      info.idx_varying = idx;
      idx += info.n_variying;
      info.idx_association = idx;
      idx += marker_info.size();
    }

    idx_shared_survival = idx;
    n_params = idx + survival_info.size() * survival_info.size();
  }

public:
  inline std::vector<marker> const & get_marker_info() const {
    return marker_info;
  }

  inline std::vector<survival> const & get_survival_info() const {
    return survival_info;
  }

  /// adds a marker to the model
  inline void add_marker(marker const &info){
    marker_info.push_back(info);
    re_compute_indices();
  }

  /// adds a survival outcome to the model
  inline void add_survival(survival const &info){
    survival_info.push_back(info);
    re_compute_indices();
  }

  /// returns the index of the fixed effect coefficients for a given marker
  inline vajoint_uint get_fixef_idx_marker(vajoint_uint const idx) const {
    return marker_info[idx].idx_fix;
  }

  /**
   * returns the index of the fixed varying effect's coefficients for a given
   * marker
   */
  inline vajoint_uint get_varying_idx_marker(vajoint_uint const idx) const {
    return marker_info[idx].idx_varying;
  }

  /// returns the index of the covarinace matrix for the markers' residual
  inline vajoint_uint get_idx_error_term() const {
    return idx_error_term;
  }

  /// returns the index of the covarinace matrix for the shared effect
  inline vajoint_uint get_idx_shared_effect() const{
    return idx_shared_effect;
  }

  /**
   * returns the index of the fixed effect coefficients for a given survival
   * outcome
   */
  inline vajoint_uint get_fixef_idx_survival(vajoint_uint const idx) const {
    return survival_info[idx].idx_fix;
  }

  /**
   * returns the index of the fixed varying effect's coefficients for a given
   * survival outcome
   */
  inline vajoint_uint get_varying_idx_survival(vajoint_uint const idx) const {
    return survival_info[idx].idx_varying;
  }

  /**
   * returns the index of the association parameter for a given survival outcome
   */
  inline vajoint_uint get_idx_association_parameter(vajoint_uint const idx)
  const {
    return survival_info[idx].idx_association;
  }

  /// returns the index of the covarinace matrix for the frailties
  inline vajoint_uint get_idx_shared_survival() const {
    return idx_shared_survival;
  }

  /// returns the number of parameters
  inline vajoint_uint get_n_params() const {
    return n_params;
  }
};

#endif
