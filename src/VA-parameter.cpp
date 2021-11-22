#include "VA-parameter.h"
#include <stdexcept>

std::vector<std::string> subset_params::param_names
  (bool const is_traingular) const {
  if(!is_traingular)
    throw std::runtime_error("not implemented");

  std::vector<std::string> out(n_params<true>());

  // handle the fixed effects
  {
    vajoint_uint m_idx{};
    for(auto &x : marker_info_v){
      ++m_idx;
      std::string const prefix{"mark" + std::to_string(m_idx)};

      for(vajoint_uint i = 0; i < x.n_fix; ++i)
        out[i + x.idx_fix] = prefix + "_fixef" + std::to_string(i + 1);

      for(vajoint_uint i = 0; i < x.n_variying; ++i)
        out[i + x.idx_varying] = prefix + "_fixef_vary" + std::to_string(i + 1);
    }
  }

  {
    vajoint_uint s_idx{};
    for(auto &x : surv_info_v){
      ++s_idx;
      std::string const prefix{"surv" + std::to_string(s_idx)};

      for(vajoint_uint i = 0; i < x.n_fix; ++i)
        out[i + x.idx_fix] = prefix + "_fixef" + std::to_string(i + 1);

      for(vajoint_uint i = 0; i < x.n_variying; ++i)
        out[i + x.idx_varying] = prefix + "_fixef_vary" + std::to_string(i + 1);

      unsigned idx_assoc{x.idx_association};
      for(vajoint_uint i = 0; i < x.n_associations.size(); ++i)
        for(vajoint_uint j = 0; j < x.n_associations[i]; ++j)
          out[idx_assoc++] = prefix + "_assoc" +
            std::to_string(i + 1) + "_" + std::to_string(j + 1);
    }
  }

  // handle the covariance matrices
  {
    vajoint_uint const dim_vcov_marker{
      static_cast<vajoint_uint>(dim_tri(marker_info_v.size()))};
    for(vajoint_uint i = 0; i < dim_vcov_marker; ++i)
      out[i + vcov_marker<true>()] = "vcov_marker" + std::to_string(i + 1);
  }
  {
    vajoint_uint const dim_vcov_surv{dim_tri(n_shared_surv_v)};
    for(vajoint_uint i = 0; i < dim_vcov_surv; ++i)
      out[i + vcov_surv<true>()] = "vcov_surv" + std::to_string(i + 1);
  }

  vajoint_uint const dim_vcov_vary{dim_tri(n_shared())};
  for(vajoint_uint i = 0; i < dim_vcov_vary; ++i)
    out[i + vcov_vary<true>()] = "vcov_vary" + std::to_string(i + 1);

  return out;

}

/// gives the names of the variational parameters
std::vector<std::string> subset_params::va_param_names
  (bool const is_traingular) const {
  if(!is_traingular)
    throw std::runtime_error("not implemented");

  std::vector<std::string> out(n_va_params<true>());

  vajoint_uint idx{}, m_idx{};
  for(auto &x : marker_info()){
    ++m_idx;
    for(vajoint_uint i = 0; i < x.n_rng; ++i)
      out[idx++] = "mark" + std::to_string(m_idx) + "_rng" + std::to_string(i + 1);
  }

  for(vajoint_uint i = 0; i < surv_info().size(); ++i)
    if(surv_info()[i].with_frailty)
      out[idx++] = "frailty" + std::to_string(i + 1);

  vajoint_uint vcov_dim(n_shared() + n_shared_surv_v);
  vcov_dim = dim_tri(vcov_dim);
  for(vajoint_uint i = 0; i < vcov_dim; ++i)
    out[idx++] = "VA_vcov" + std::to_string(i + 1);

  return out;
}
