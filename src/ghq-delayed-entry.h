#include "ghq.h"
#include "simple-mat.h"
#include "bases.h"
#include "VA-parameter.h"
#include "JointSurv-misc.h"
#include <numeric>

namespace survival {

struct delayed_dat {
  struct cluster_obs {
    /// the type of survival outcome and the index among those types
    vajoint_uint type, index;
    /// the entry time
    double entry_time;
  };
  /// defines all the delayed entries for a cluster
  using cluster_info = std::vector<cluster_obs>;

private:
  /**
   * the bases for the time-varying fixed effects (one for each type of
   * outcome)
   */
  joint_bases::bases_vector bases_fix;
  /// the bases for the time-varying random effects
  joint_bases::bases_vector bases_rng;
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
  /// the number of basis functions of each rng base expansion
  std::vector<vajoint_uint> rng_n_basis_v {
    ([&]{
      std::vector<vajoint_uint> out;
      out.reserve(bases_rng.size());
      for(auto &b : bases_rng)
        out.emplace_back(b->n_basis());
      return out;
    })()
  };
  /// the cumulative sum of the number weights starting at zero
  std::vector<vajoint_uint> rng_n_weights_cumsum_v{
    ([&]{
      std::vector<vajoint_uint> out;
      out.reserve(bases_rng.size() + 1L);
      out.emplace_back(0);
      for(auto &b : bases_rng)
        out.emplace_back(out.back() + b->n_weights());
      return out;
    })()
  };

  /// like calling bases_rng[idx]->n_basis()
  vajoint_uint rng_n_basis(size_t const idx) const {
    return rng_n_basis_v[idx];
  }

  vajoint_uint rng_n_weights_cumsum(size_t const idx) const {
    return rng_n_weights_cumsum_v[idx];
  }

  /// the number of fixed effects
  vajoint_uint n_fixef;
  /// the derivative/integral argument to pass to bases_rng
  std::vector<std::vector<std::vector<int> > > ders_v;

  /// the indices of the parameters
  subset_params par_idx;

  /// the info needed to compute each cluster
  std::vector<cluster_info> v_cluster_infos;

  /**
   * vector with index from the 1,...,E survival type to the index of the
   * frailty. It can be random accessed at index 0,...,E-1
   */
  std::vector<vajoint_uint> frailty_map;

  /// helper class to evaluate the expected hazard and to compute gradient of it
  struct eval_data;
  friend class eval_data;
  struct eval_data {
    /// the scaled quadrature weights for each survival outcome in the cluster
    std::vector<double> quad_weights;

    /**
     * suppose that there are l survival outcomes in the cluster. This contains
     * the evaluated basis expansion at each quadrature node in the format
     *
     *     basis dim
     *     ~~~~~~~~~
     *   o|---------|
     *   u|---------|
     *   t|---------|
     *   1|---------|
     *     ~~~~~~~~~
     *   .|----|
     *   .|----|
     *   .|----|
     *     ~~~~~~~
     *   o|-------|
     *   u|-------|
     *   t|-------|
     *   l|-------|
     *     ~~~~~~~
     *
     * The dimension of the basis functions may differ if the survival outcomes
     * are of different types. Thus, we store the blocks as a vector of
     * matrices.
     */
    std::vector<simple_mat<double> > fixef_vary_basis;

    /**
     * input to construct the random effect design matrix. The random effect
     * design matrix, disregarding the frailty, is of the form
     *
     *    marker 1   ...    marker k
     *    ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * o |--------|--------|--------|
     * u |--------|--------|--------|
     * t |--------|--------|--------|
     * 1 |--------|--------|--------|
     *    ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * . |--------|--------|--------|
     * . |--------|--------|--------|
     * . |--------|--------|--------|
     * . |--------|--------|--------|
     *    ~~~~~~~~~~~~~~~~~~~~~~~~~~
     * o |--------|--------|--------|
     * u |--------|--------|--------|
     * t |--------|--------|--------|
     * l |--------|--------|--------|
     *    ~~~~~~~~~~~~~~~~~~~~~~~~~~
     *
     * Each of the l x k blocks consist of the sum of rows of basis expansions
     * scaled by the association parameters. Thus, we have to store each block
     * which as
     *
     * vector <-- marker 1 to k
     *   vector <-- outcome 1 to l
     *     vector <-- the associations for outcome l and marker k
     *       matrix <-- blocks for each type of association. Dimensions are
     *                  <n quadrature nodes> <dimension of the basis>
     *
     * TODO: a multidimensional arrays will likely be more efficient. All but
     *       the association dimension are equal and that dimension could be
     *       dealt with using some zero padding.
     */
    std::vector<std::vector<std::vector<simple_mat<double> > > > rng_basis;

    eval_data
      (delayed_dat const &dat, node_weight const &nws,
       cluster_info const &info, ghqCpp::simple_mem_stack<double> &mem);

    /**
     * if there are E survival types then this is a vector with E elements that
     * can be random accessed with the the indices of "active" frailties
     * variables in the cluster for each type (a subset of (1, ..., E)).
     */
    std::vector<vajoint_uint> idx_active_frailty;

    /// the inverse of idx_active_frailty. Thus, the size may be less than E
    std::vector<vajoint_uint> idx_inv_active_fraitly;

    size_t n_active_frailties() const { return idx_inv_active_fraitly.size(); }
  };

  vajoint_uint n_markers() const { return bases_rng.size(); }

  /// struct to hide parts of the implementation
  struct impl;
  friend class impl;

  /// holds memory for the cached data
  std::vector<eval_data> cached_expansions;

  /// the cached quadrature nodes and weights
  std::vector<double> cached_nodes, cached_weights;

  bool has_cached_expansions() const {
    return cached_expansions.size() > 0;
  }

public:
  delayed_dat() = default;

  delayed_dat(joint_bases::bases_vector const &bases_fix_in,
              joint_bases::bases_vector const &bases_rng_in,
              std::vector<simple_mat<double> > &design_mats,
              std::vector<simple_mat<double> > &fixef_design_varying,
              std::vector<simple_mat<double> > &rng_design_varying,
              subset_params const &par_idx,
              std::vector<cluster_info> const &cluster_infos,
              std::vector<std::vector<std::vector<int> > > &ders);

  void set_cached_expansions
    (node_weight const &nws, ghqCpp::simple_mem_stack<double> &mem);

  /// clears the cached expansions
  void clear_cached_expansions();

  /// returns information about each cluster
  std::vector<cluster_info> const & cluster_infos() const {
    return v_cluster_infos;
  }

  /// evaluates the delayed entry term for a given cluster
  double operator()
    (double const *param, ghqCpp::simple_mem_stack<double> &mem,
     const vajoint_uint cluster_index, node_weight const &nws,
     ghqCpp::ghq_data const &ghq_dat) const;

  /**
   * evaluates the delayed entry term for a given cluster and adds the gradient
   * to the passed pointer (i.e. the pointer is not overwritten)
   */
  double grad
    (double const *param, double *gr, ghqCpp::simple_mem_stack<double> &mem,
     const vajoint_uint cluster_index, node_weight const &nws,
     ghqCpp::ghq_data const &ghq_dat) const;
};

} // namespace survival
