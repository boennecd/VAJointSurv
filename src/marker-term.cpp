#include "marker-term.h"
#include <algorithm>

namespace marker {

comp_dat::comp_dat
(double const *param, double *wk_mem, subset_params const &par_idx,
 std::uint32_t const missingness_flag):
indices{([&]{
    vajoint_uint const n_markers = par_idx.marker_info().size();
    // find the non-missing variables
    std::vector<vajoint_uint> out;
    if(missingness_flag == 0){
      // all present
      out.resize(n_markers);
      std::iota(out.begin(), out.end(), 0);

    } else {
      // some are missing
      out.reserve(n_markers);
      unsigned v{1};
      for(unsigned i = 0; i < n_markers; ++i, v *= 2)
        if(!(missingness_flag & v))
          out.emplace_back(i);
    }

    return out;
})()},
vcov_factorization{([&]() -> cfaad::CholFactorization {
  /// create the subset of the covariance matrix
  vajoint_uint const n_markers = par_idx.marker_info().size();
  double const * const sig = param + par_idx.vcov_marker();
  size_t const n_indices{indices.size()};
  for(vajoint_uint j = 0; j < n_indices; ++j)
    for(vajoint_uint i = 0; i < n_indices; ++i)
      wk_mem[i + j * n_indices] = sig[indices[i] + n_markers * indices[j]];

  return { wk_mem, static_cast<int>(n_indices), true };
})()},
n_rngs{([&]{
  vajoint_uint out{};
  for(vajoint_uint idx : indices)
    out += par_idx.marker_info()[idx].n_rng;
  return out;
})()}{ }

comp_dat_return get_comp_dat
(std::vector<setup_marker_dat_helper> &input_dat,
 subset_params const &par_idx,
 joint_bases::bases_vector const &bases_fix,
 joint_bases::bases_vector const &bases_rng){
  // basic checks
  if(input_dat.size() != bases_fix.size())
    throw std::invalid_argument("input_dat.size() != bases_fix.size()");
  if(input_dat.size() != bases_rng.size())
    throw std::invalid_argument("input_dat.size() != bases_rng.size()");

  // compute the number of observations
  vajoint_uint const n_markers(input_dat.size());

  // compute the maximum number of unique observations
  vajoint_uint max_obs{};
  for(auto &x : input_dat)
    max_obs += x.n_obs();

  // the id for each unique observation
  std::vector<int> unique_ids;
  unique_ids.reserve(max_obs);
  // the observation time for each unique observation
  std::vector<double> unique_obs_time;
  unique_obs_time.reserve(max_obs);

  // pointer to each markers current id, observation time, and current index
  std::vector<int const*> ids_ptr(n_markers);
  std::vector<double const*> obs_time_ptr(n_markers);
  std::vector<vajoint_uint> idx(n_markers);

  // fills the pointers and indices
  auto fill_ptrs = [&]{
    for(vajoint_uint i = 0; i < n_markers; ++i){
      ids_ptr[i] = input_dat[i].ids;
      obs_time_ptr[i] = input_dat[i].obs_time;
      idx[i] = 0;
    }
  };
  fill_ptrs(); // fill already

  // increments the i'th set of pointers
  auto inc_ptrs = [&](vajoint_uint const i){
    ++ids_ptr[i];
    ++obs_time_ptr[i];
    ++idx[i];
  };

  // pointer to end of ids_ptr
  std::vector<int const*> ids_ptr_end(n_markers);
  for(vajoint_uint i = 0; i < n_markers; ++i)
    ids_ptr_end[i] = input_dat[i].ids + input_dat[i].n_obs();

  // returns true if all pointers are at the end
  auto all_finished = [&]{
    for(vajoint_uint i = 0; i < n_markers; ++i)
      if(ids_ptr[i] != ids_ptr_end[i])
        return false;
    return true;
  };

  // compute the number of observations
  for(; !all_finished();){
    int min_id{std::numeric_limits<int>::max()};
    double min_time{std::numeric_limits<double>::max()};

    for(vajoint_uint i = 0; i < n_markers; ++i){
      if(ids_ptr[i] == ids_ptr_end[i])
        continue;
      if(*ids_ptr[i] < min_id ||
          (*obs_time_ptr[i] < min_time && *ids_ptr[i] == min_id)){
        // update the minimums
        min_id = *ids_ptr[i];
        min_time = *obs_time_ptr[i];
      }
    }

    unique_ids.emplace_back(min_id);
    unique_obs_time.emplace_back(min_time);

    // increment the pointers at the minimums
    for(vajoint_uint i = 0; i < n_markers; ++i){
      if(ids_ptr[i] == ids_ptr_end[i] ||
         *ids_ptr[i] > min_id ||
         *obs_time_ptr[i] > min_time)
         continue;

      inc_ptrs(i);
    }
  }

  // create the object and initialize
  vajoint_uint const n_obs(unique_ids.size());
  marker_dat out{par_idx, n_obs, bases_fix, bases_rng};

  // fill in the design matrices for the time-varying effects
  out.set_design_mats(unique_obs_time.begin());

  // fill in the outcomes and fixed effect design matrix
  fill_ptrs();
  for(vajoint_uint i = 0; i < n_obs; ++i){
    for(vajoint_uint j = 0; j < n_markers; ++j){
      if(ids_ptr[j] == ids_ptr_end[j] ||
         *ids_ptr[j] > unique_ids[i] ||
         *obs_time_ptr[j] > unique_obs_time[i])
         continue;

      out.set_outcome(i, j, input_dat[j].obs[idx[j]]);
      out.set_fixef_design(i, j, input_dat[j].fixef_design.col(idx[j]));
      inc_ptrs(j);
    }
  }

  return { std::move(out), std::move(unique_ids) };
}
} // namespace marker
