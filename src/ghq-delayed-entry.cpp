#include "ghq-delayed-entry.h"
#include "integrand-expected-survival.h"
#include <numeric>
#include <set>

namespace survival {

delayed_dat::delayed_dat
  (joint_bases::bases_vector const &bases_fix_in,
   joint_bases::bases_vector const &bases_rng_in,
   std::vector<simple_mat<double> > &design_mats, subset_params const &par_idx,
   std::vector<cluster_info> const &cluster_infos,
   std::vector<std::vector<std::vector<int> > > &ders):
  bases_fix{joint_bases::clone_bases(bases_fix_in)},
  bases_rng{joint_bases::clone_bases(bases_rng_in)},
  design_mats{design_mats},
  ders_v{ders},
  par_idx{par_idx},
  v_cluster_infos{cluster_infos},
  frailty_map{
    ([&]{
      vajoint_uint const n_types = par_idx.surv_info().size();
      std::vector<vajoint_uint> out(n_types);
      vajoint_uint idx{};
      for(size_t i = 0; i < n_types; ++i)
        if(par_idx.surv_info()[i].with_frailty)
          out[i] = idx++;
        else
          out[i] = n_types;
      return out;
    })()
  } {
    if(bases_fix_in.size() != par_idx.surv_info().size())
      throw std::invalid_argument
        ("bases_fix_in.size() != par_idx.surv_info().size()");
    else if(bases_rng_in.size() != par_idx.marker_info().size())
      throw std::invalid_argument
        ("bases_rng_in.size() != par_idx.marker_info().size()");
    else if(design_mats.size() != par_idx.surv_info().size())
      throw std::invalid_argument
        ("design_mats.size() != par_idx.surv_info().size()");
    for(size_t i = 0; i < par_idx.surv_info().size(); ++i)
      if(design_mats[i].n_rows() != par_idx.surv_info()[i].n_fix)
        throw std::invalid_argument
          ("design_mats[i].n_rows() != par_idx.surv_info()[i].n_fix");

    for(auto &c_info : cluster_infos)
      for(auto &obs : c_info)
        if(obs.type >= par_idx.surv_info().size())
          throw std::invalid_argument
            ("obs.type >= par_idx.surv_info().size()");
        else if(obs.index >= design_mats[obs.type].n_cols())
          throw std::invalid_argument
            ("obs.index >= design_mats[obs.type].n_cols()");

    for(size_t i = 0; i < ders.size(); ++i){
      if(ders[i].size() != bases_rng_in.size())
        throw std::invalid_argument("ders_i.size() != bases_rng_in.size()");

      for(size_t j = 0; j < ders[i].size(); ++j)
        if(ders[i][j].size() != par_idx.surv_info()[i].n_associations[j])
          throw std::invalid_argument
            ("ders[i][j].size() != par_idx.surv_info()[i].n_associations[j]");
    }
  }

delayed_dat::eval_data::eval_data
  (delayed_dat const &dat, node_weight const &nws,
   delayed_dat::cluster_info const &info,
   ghqCpp::simple_mem_stack<double> &mem){

  vajoint_uint const n_outcomes = info.size(),
                           n_gl = nws.n_nodes,
                  n_gl_outcomes{n_gl * n_outcomes},
                      n_markers = dat.n_markers();

  // set the quadrature weights for the cumulative hazards and the time points
  // at which we evaluate the bases
  double * const time_points(mem.get(n_gl_outcomes));
  quad_weights.reserve(n_gl_outcomes);

  {
    double *time_points_i{time_points};
    for(auto &obs : info)
      for(vajoint_uint i = 0; i < n_gl; ++i){
        *time_points_i++ = obs.entry_time * nws.ns[i];
        quad_weights.emplace_back(obs.entry_time * nws.ws[i]);
      }
  }

  // fill in the design matrix for the time-varying fixed effects
  fixef_vary_basis.reserve(n_outcomes);
  {
    double *time_points_i{time_points};
    for(auto &obs : info){
      auto const &fixef_base = dat.bases_fix[obs.type];
      auto const n_basis = fixef_base->n_basis();
      double * const basis_mem{mem.get(fixef_base->n_wmem() + n_basis)},
             * const basis_mem_wk{basis_mem + n_basis};

      fixef_vary_basis.emplace_back(n_gl, n_basis);
      auto &mat = fixef_vary_basis.back();
      for(vajoint_uint i = 0; i < n_gl; ++i, ++time_points_i){
        (*fixef_base)(basis_mem, basis_mem_wk, *time_points_i);

        // copy the transpose
        for(vajoint_uint j = 0; j < n_basis; ++j)
          mat.col(j)[i] = basis_mem[j];
      }
    }
  }

  // fill in the design matrix for the random effects
  rng_basis.reserve(n_markers);

  for(vajoint_uint k = 0; k < n_markers; ++k){
    rng_basis.emplace_back();
    std::vector<std::vector<simple_mat<double> > > &marker_k{rng_basis.back()};
    marker_k.reserve(n_outcomes);

    auto const &fixef_base = dat.bases_rng[k];
    auto const n_basis = fixef_base->n_basis();
    double * const basis_mem{mem.get(fixef_base->n_wmem() + n_basis)},
           * const basis_mem_wk{basis_mem + n_basis};

    double *time_points_i{time_points};
    for(vajoint_uint l = 0; l < n_outcomes; ++l, time_points_i += n_gl){
      auto &obs = info[l];

      auto &ders_kl = dat.ders_v[obs.type][k];
      marker_k.emplace_back(ders_kl.size(), simple_mat<double>{n_gl, n_basis});
      std::vector<simple_mat<double> > &marker_kl{marker_k.back()};

      for(size_t i = 0; i < ders_kl.size(); ++i){
        auto &simple_mat_i = marker_kl[i];
        for(size_t h = 0; h < n_gl; ++h){
          (*fixef_base)(basis_mem, basis_mem_wk, time_points_i[h], ders_kl[i]);

          // copy the transpose
          for(vajoint_uint j = 0; j < n_basis; ++j)
            simple_mat_i.col(j)[h] = basis_mem[j];
        }
      }
    }
  }

  // fill in the maps to active frailties
  std::set<vajoint_uint> type_active;
  for(auto &obs : info)
    if(dat.par_idx.surv_info()[obs.type].with_frailty)
      type_active.emplace(obs.type);

  idx_active_frailty.resize(dat.par_idx.surv_info().size());
  std::fill(idx_active_frailty.begin(),
            idx_active_frailty.end() + dat.par_idx.surv_info().size(), 0);

  idx_inv_active_fraitly.resize(type_active.size());

  auto type = type_active.begin();
  for(size_t i = 0; i < type_active.size(); ++i, ++type){
    idx_active_frailty[*type] = i;
    idx_inv_active_fraitly[i] = dat.frailty_map[*type];
  }
}

double delayed_dat::operator()
  (double const *param, ghqCpp::simple_mem_stack<double> &mem,
   const vajoint_uint cluster_index, node_weight const &nws,
   ghqCpp::ghq_data const &ghq_dat) const {
  auto const &info{cluster_infos()[cluster_index]};
  vajoint_uint const n_outcomes = info.size(),
                     n_gl = nws.n_nodes,
                     n_gl_outcomes{n_gl * n_outcomes},
                     n_shared{par_idx.n_shared()};

  eval_data e_dat{*this, nws, info, mem}; // TODO: can cache this
  vajoint_uint const n_rng = n_shared + e_dat.n_active_frailties();

  // get the memory we need
  double * const etas
    {mem.get(n_gl_outcomes + n_gl_outcomes * n_rng + n_rng * n_rng)},
         * const rng_design{etas + n_gl_outcomes},
         * const vcov{rng_design + n_gl_outcomes * n_rng};
  auto mem_marker = mem.set_mark_raii();

  // setup the eta offsets. These consists of time-varying and fixed part
  // We start with the fixed part
  for(vajoint_uint i = 0; i < n_outcomes; ++i){
    auto &obs = info[i];
    auto &design_mat_i = design_mats[obs.type];
    double const *x{design_mat_i.col(obs.index)};
    double const offset
      {std::inner_product(x, x + design_mat_i.n_rows(),
                          param + par_idx.fixef_surv(obs.type), 0.)};
    std::fill(etas + i * n_gl, etas + (i + 1) * n_gl, offset);
  }

  // then the time-varying part
  {
    double *etas_i{etas};
    for(size_t l = 0; l < info.size(); ++l, etas_i += n_gl){
      auto &obs = info[l];
      double const * const fixef_vary_par
        {param + par_idx.fixef_vary_surv(obs.type)};

      auto &fixef_vary_basis_l = e_dat.fixef_vary_basis[l];
      for(vajoint_uint j = 0; j < fixef_vary_basis_l.n_cols(); ++j)
        for(vajoint_uint i = 0; i < fixef_vary_basis_l.n_rows(); ++i)
          etas_i[i] += fixef_vary_basis_l.col(j)[i] * fixef_vary_par[j];
    }
  }

  // fill the rng_design matrix
  std::fill(rng_design, rng_design + n_gl_outcomes * n_rng, 0);
  {
    double * rng_design_k{rng_design};
    for(vajoint_uint k = 0; k < n_markers();
        rng_design_k += rng_n_basis(k) * n_gl_outcomes, ++k){
      std::vector<std::vector<simple_mat<double> > > &marker_k =
        e_dat.rng_basis[k];

      size_t const n_basis = rng_n_basis(k);

      double * rng_design_kl{rng_design_k};
      for(size_t l = 0; l < info.size(); ++l, rng_design_kl += n_gl){
        auto &obs = info[l];
        auto const &rng_basis_kl = marker_k[l];
        double const * association{param + par_idx.association(obs.type)};

        // find the right association parameters
        auto &ders_l = ders_v[obs.type];
        for(size_t kk = 0; kk < k; ++kk)
          association += ders_l[kk].size();

        for(size_t v = 0; v < rng_basis_kl.size(); ++v){
          auto const &rng_basis_klv = rng_basis_kl[v];
          for(size_t j = 0; j < n_basis; ++j)
            for(size_t i = 0; i < n_gl; ++i)
              rng_design_kl[i + j * n_gl_outcomes] +=
                association[v] * rng_basis_klv.col(j)[i];
        }
      }
    }
  }

  // fill in the frailty columns
  if(e_dat.n_active_frailties()){
    double * rng_design_frailty_part{rng_design + n_gl_outcomes * n_shared};

    for(size_t l = 0; l < info.size(); ++l, rng_design_frailty_part += n_gl){
      auto &obs = info[l];
      if(!par_idx.surv_info()[obs.type].with_frailty)
        continue;

      auto idx = e_dat.idx_active_frailty[obs.type];
      double *cp{rng_design_frailty_part + idx * n_gl_outcomes};
      std::fill(cp, cp + n_gl, 1);
    }
  }

  // fill the covariance matrix
  {
    double * vcov_j{vcov};
    double const * vcov_vary{param + par_idx.vcov_vary()};
    for(size_t j = 0; j < n_shared; ++j, vcov_j += n_rng,
        vcov_vary += n_shared){
      std::copy(vcov_vary, vcov_vary + n_shared, vcov_j);
      std::fill(vcov_j + n_shared, vcov_j + n_rng, 0);
    }

    double const * const vcov_surv{param + par_idx.vcov_surv()};
    auto const n_shared_surv = par_idx.n_shared_surv();
    size_t const n_left{e_dat.n_active_frailties()};

    for(size_t j = 0; j < n_left; ++j, vcov_j += n_rng){
      std::fill(vcov_j, vcov_j + n_shared, 0);

      size_t const offset{e_dat.idx_inv_active_fraitly[j] * n_shared_surv};
      for(size_t i = 0; i < e_dat.idx_inv_active_fraitly.size(); ++i)
        vcov_j[i + n_shared] =
          vcov_surv[e_dat.idx_inv_active_fraitly[i] + offset];
    }
  }

  arma::vec ws_vec(&e_dat.quad_weights[0], n_gl_outcomes, false),
          etas_vec(etas, n_gl_outcomes, false);
  arma::mat rng_design_mat(rng_design, n_gl_outcomes, n_rng, false),
                  vcov_mat(vcov, n_rng, n_rng, false);
  ghqCpp::expected_survival_term<false> surv_term
    (etas_vec, ws_vec, rng_design_mat, vcov_mat);
  ghqCpp::adaptive_problem prob(surv_term, mem, 1e-6);

  return ghqCpp::ghq(ghq_dat, prob, mem, 1000)[0];
}

}
