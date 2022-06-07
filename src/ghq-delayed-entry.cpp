#include "ghq-delayed-entry.h"
#include "integrand-expected-survival.h"
#include <numeric>
#include <set>

namespace survival {

delayed_dat::delayed_dat
  (joint_bases::bases_vector const &bases_fix_in,
   joint_bases::bases_vector const &bases_rng_in,
   std::vector<simple_mat<double> > &design_mats,
   std::vector<simple_mat<double> > &fixef_design_varying_mats,
   std::vector<simple_mat<double> > &rng_design_varying_mats,
   subset_params const &par_idx,
   std::vector<cluster_info> const &cluster_infos,
   std::vector<std::vector<std::vector<int> > > &ders):
  bases_fix{joint_bases::clone_bases(bases_fix_in)},
  bases_rng{joint_bases::clone_bases(bases_rng_in)},
  design_mats{design_mats},
  fixef_design_varying_mats{fixef_design_varying_mats},
  rng_design_varying_mats{rng_design_varying_mats},
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
    else if(fixef_design_varying_mats.size() != par_idx.surv_info().size())
      throw std::invalid_argument
      ("fixef_design_varying_mats.size() != par_idx.surv_info().size()");
    else if(rng_design_varying_mats.size() != par_idx.surv_info().size())
      throw std::invalid_argument
      ("rng_design_varying_mats.size() != par_idx.surv_info().size()");

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
        else if(obs.index >= fixef_design_varying_mats[obs.type].n_cols())
          throw std::invalid_argument
          ("obs.index >= fixef_design_varying_mats[obs.type].n_cols()");
        else if(obs.index >= rng_design_varying_mats[obs.type].n_cols())
          throw std::invalid_argument
          ("obs.index >= rng_design_varying_mats[obs.type].n_cols()");

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
    for(size_t obs_idx = 0; obs_idx < info.size(); ++obs_idx){
      auto &obs = info[obs_idx];
      auto const &fixef_base = dat.bases_fix[obs.type];
      auto const n_basis = fixef_base->n_basis();
      double * const basis_mem{mem.get(fixef_base->n_wmem() + n_basis)},
             * const basis_mem_wk{basis_mem + n_basis};

      fixef_vary_basis.emplace_back(n_gl, n_basis);
      auto &mat = fixef_vary_basis.back();
      for(vajoint_uint i = 0; i < n_gl; ++i, ++time_points_i){
        (*fixef_base)
          (basis_mem, basis_mem_wk, *time_points_i,
           dat.fixef_design_varying_mats[obs.type].col(obs.index));

        // copy the transpose
        for(vajoint_uint j = 0; j < n_basis; ++j)
          mat.col(j)[i] = basis_mem[j];
      }
    }
  }

  // fill in the design matrix for the random effects
  rng_basis.reserve(n_markers);

  for(vajoint_uint mark = 0; mark < n_markers; ++mark){
    rng_basis.emplace_back();
    std::vector<std::vector<simple_mat<double> > > &marker_k{rng_basis.back()};
    marker_k.reserve(n_outcomes);

    auto const &rng_base = dat.bases_rng[mark];
    auto const n_basis = rng_base->n_basis();
    double * const basis_mem{mem.get(rng_base->n_wmem() + n_basis)},
           * const basis_mem_wk{basis_mem + n_basis};

    double *time_points_i{time_points};
    for(vajoint_uint obs_idx = 0; obs_idx < n_outcomes;
        ++obs_idx, time_points_i += n_gl){
      auto &obs = info[obs_idx];

      auto &ders_kl = dat.ders_v[obs.type][mark];
      marker_k.emplace_back(ders_kl.size(), simple_mat<double>{n_gl, n_basis});
      std::vector<simple_mat<double> > &marker_kl{marker_k.back()};

      for(size_t der = 0; der < ders_kl.size(); ++der){
        auto &simple_mat_i = marker_kl[der];
        for(size_t h = 0; h < n_gl; ++h){
          (*rng_base)
            (basis_mem, basis_mem_wk, time_points_i[h],
             dat.rng_design_varying_mats[obs.type].col(obs.index) +
               dat.rng_n_weights_cumsum(mark),
             ders_kl[der]);

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
  std::fill(idx_active_frailty.begin(), idx_active_frailty.end(), 0);

  idx_inv_active_fraitly.resize(type_active.size());

  auto type = type_active.begin();
  for(size_t i = 0; i < type_active.size(); ++i, ++type){
    idx_active_frailty[*type] = i;
    idx_inv_active_fraitly[i] = dat.frailty_map[*type];
  }
}

void delayed_dat::set_cached_expansions
  (node_weight const &nws, ghqCpp::simple_mem_stack<double> &mem){
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
  cached_expansions.reserve(cluster_info().size());

  for(auto &info : cluster_infos()){
    mem.reset_to_mark();
    cached_expansions.emplace_back(*this, nws, info, mem);
  }
}

void delayed_dat::clear_cached_expansions(){
  cached_expansions.clear();
  cached_expansions.shrink_to_fit();

  cached_nodes.clear();
  cached_nodes.shrink_to_fit();

  cached_weights.clear();
  cached_weights.shrink_to_fit();
}

struct delayed_dat::impl {
  delayed_dat const &dat;
  delayed_dat::cluster_info const &info;
  node_weight const &nws;
  eval_data const &e_dat;
  ghqCpp::simple_mem_stack<double> &mem;

  vajoint_uint const n_outcomes = info.size(),
                     n_gl = nws.n_nodes,
                     n_gl_outcomes{n_gl * n_outcomes},
                     n_shared{dat.par_idx.n_shared()},
                     n_rng = n_shared + e_dat.n_active_frailties();

  double * const etas
    {mem.get(n_gl_outcomes + n_gl_outcomes * n_rng + n_rng * n_rng)},
         * const rng_design{etas + n_gl_outcomes},
         * const vcov{rng_design + n_gl_outcomes * n_rng};

  ghqCpp::simple_mem_stack<double>::return_memory_handler mem_mark
    {mem.set_mark_raii()};

  impl(delayed_dat const &dat, delayed_dat::cluster_info const &info,
       node_weight const &nws, eval_data const &e_dat,
       ghqCpp::simple_mem_stack<double> &mem, double const *param):
    dat{dat}, info{info}, nws{nws}, e_dat{e_dat}, mem{mem} {
    	// setup the eta offsets. These consists of time-varying and fixed part
      for(vajoint_uint i = 0; i < n_outcomes; ++i){
        // start with the fixed part
        auto &obs = info[i];
        double * const etas_i{etas + i * n_gl};
        {
          auto &design_mat_i = dat.design_mats[obs.type];
          double const *x{design_mat_i.col(obs.index)};
          double const offset
            {std::inner_product(x, x + design_mat_i.n_rows(),
                                param + dat.par_idx.fixef_surv(obs.type), 0.)};
          std::fill(etas_i, etas_i + n_gl, offset);
        }

        double const * const fixef_vary_par
          {param + dat.par_idx.fixef_vary_surv(obs.type)};

        // then the time-varying part
        auto &fixef_vary_basis_i = e_dat.fixef_vary_basis[i];
        for(vajoint_uint j = 0; j < fixef_vary_basis_i.n_cols(); ++j)
          for(vajoint_uint l = 0; l < fixef_vary_basis_i.n_rows(); ++l)
            etas_i[l] += fixef_vary_basis_i.col(j)[l] * fixef_vary_par[j];
      }

      // fill the rng_design matrix
      std::fill(rng_design, rng_design + n_gl_outcomes * n_rng, 0);
      {
        double * rng_design_k{rng_design};
        for(vajoint_uint k = 0; k < dat.n_markers();
            rng_design_k += dat.rng_n_basis(k) * n_gl_outcomes, ++k){
          std::vector<std::vector<simple_mat<double> > > const &marker_k =
            e_dat.rng_basis[k];

          size_t const n_basis = dat.rng_n_basis(k);

          double * rng_design_kl{rng_design_k};
          for(size_t l = 0; l < info.size(); ++l, rng_design_kl += n_gl){
            auto &obs = info[l];
            auto const &rng_basis_kl = marker_k[l];
            double const * association{param + dat.par_idx.association(obs.type)};

            // find the right association parameters
            auto &ders_l = dat.ders_v[obs.type];
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
          if(!dat.par_idx.surv_info()[obs.type].with_frailty)
            continue;

          auto idx = e_dat.idx_active_frailty[obs.type];
          double *cp{rng_design_frailty_part + idx * n_gl_outcomes};
          std::fill(cp, cp + n_gl, 1);
        }
      }

      // fill the covariance matrix
      double * vcov_j{vcov};
      double const * vcov_vary{param + dat.par_idx.vcov_vary()};
      for(size_t j = 0; j < n_shared; ++j, vcov_j += n_rng,
          vcov_vary += n_shared){
        std::copy(vcov_vary, vcov_vary + n_shared, vcov_j);
        std::fill(vcov_j + n_shared, vcov_j + n_rng, 0);
      }

      double const * const vcov_surv{param + dat.par_idx.vcov_surv()};
      auto const n_shared_surv = dat.par_idx.n_shared_surv();
      size_t const n_left{e_dat.n_active_frailties()};

      for(size_t j = 0; j < n_left; ++j, vcov_j += n_rng){
        std::fill(vcov_j, vcov_j + n_shared, 0);

        size_t const offset{e_dat.idx_inv_active_fraitly[j] * n_shared_surv};
        for(size_t i = 0; i < e_dat.idx_inv_active_fraitly.size(); ++i)
          vcov_j[i + n_shared] =
            vcov_surv[e_dat.idx_inv_active_fraitly[i] + offset];
      }
    }
};

double delayed_dat::operator()
  (double const *param, ghqCpp::simple_mem_stack<double> &mem,
   const vajoint_uint cluster_index, node_weight const &nws,
   ghqCpp::ghq_data const &ghq_dat) const {
  // use the helper to set up the objects need for the quadrature
  auto const &info{cluster_infos()[cluster_index]};

  std::unique_ptr<eval_data> local_eval_data;
  if(!has_cached_expansions())
    local_eval_data.reset(new eval_data{*this, nws, info, mem});
  eval_data const &e_dat
    {has_cached_expansions()
      ? cached_expansions[cluster_index]
      : *local_eval_data.get()};

  impl im{*this, info, nws, e_dat, mem, param};

  // apply the quadrature
  double * const etas{im.etas},
         * const rng_design{im.rng_design},
         * const vcov{im.vcov};

  vajoint_uint const n_gl_outcomes{im.n_gl_outcomes},
                     n_rng{im.n_rng};

  arma::vec ws_vec
      (const_cast<double*>(&e_dat.quad_weights[0]), n_gl_outcomes, false),
          etas_vec(etas, n_gl_outcomes, false);
  arma::mat rng_design_mat(rng_design, n_gl_outcomes, n_rng, false),
                  vcov_mat(vcov, n_rng, n_rng, false);
  ghqCpp::expected_survival_term<false> surv_term_inner
    (etas_vec, ws_vec, rng_design_mat);
  ghqCpp::rescale_problem<false> surv_term(vcov_mat, surv_term_inner);
  ghqCpp::adaptive_problem prob(surv_term, mem, 1e-6);

  double res{};
  ghqCpp::ghq(&res, ghq_dat, prob, mem, 200);
  return std::log(res);
}

double delayed_dat::grad
  (double const *param, double *gr, ghqCpp::simple_mem_stack<double> &mem,
   const vajoint_uint cluster_index, node_weight const &nws,
   ghqCpp::ghq_data const &ghq_dat) const {
  // use the helper to set up the objects need for the quadrature
  auto const &info{cluster_infos()[cluster_index]};

  std::unique_ptr<eval_data> local_eval_data;
  if(!has_cached_expansions())
    local_eval_data.reset(new eval_data{*this, nws, info, mem});
  eval_data const &e_dat
    {has_cached_expansions()
      ? cached_expansions[cluster_index]
      : *local_eval_data.get()};

  impl im{*this, info, nws, e_dat, mem, param};

  // apply the quadrature
  double * const etas{im.etas},
         * const rng_design{im.rng_design},
         * const vcov{im.vcov};

  vajoint_uint const n_gl_outcomes{im.n_gl_outcomes},
                             n_rng{im.n_rng},
                        n_outcomes{im.n_outcomes},
                              n_gl{im.n_gl},
                          n_shared{im.n_shared};

  arma::vec ws_vec
      (const_cast<double*>(&e_dat.quad_weights[0]), n_gl_outcomes, false),
  etas_vec(etas, n_gl_outcomes, false);
  arma::mat rng_design_mat(rng_design, n_gl_outcomes, n_rng, false),
  vcov_mat(vcov, n_rng, n_rng, false);
  ghqCpp::expected_survival_term<true> surv_term_inner
    (etas_vec, ws_vec, rng_design_mat);
  ghqCpp::rescale_problem<true> surv_term(vcov_mat, surv_term_inner);
  ghqCpp::adaptive_problem prob(surv_term, mem, 1e-6);

  size_t const n_res{prob.n_out()};
  double * __restrict__ res{mem.get(n_res)};
  auto mem_mark = mem.set_mark_raii();
  ghqCpp::ghq(res, ghq_dat, prob, mem, 200);
  double const fn_exp{res[0]},
               fn{std::log(fn_exp)};

  // to get the derivative w.r.t. the logarithm
  for(size_t i = 1; i < n_res; ++i)
    res[i] /= fn_exp;

  double const * const d_eta{res + 1},
               * const d_rng_design{d_eta + n_gl_outcomes};

  // handle the derivatives w.r.t. eta
  {
    double const * d_eta_i{d_eta};
    for(vajoint_uint i = 0; i < n_outcomes; ++i, d_eta_i += n_gl){
      // start with the fixed effect part
      auto &obs = info[i];
      {
        auto &design_mat_i = design_mats[obs.type];
        double const *x{design_mat_i.col(obs.index)};
        double const sum_d_eta = std::accumulate(d_eta_i, d_eta_i + n_gl, 0.);
        double * const outcome{gr + par_idx.fixef_surv(obs.type)};
        for(vajoint_uint j = 0; j < design_mat_i.n_rows(); ++j)
          outcome[j] += sum_d_eta * x[j];
      }

      // then the time-varying part
      double * const outcome{gr + par_idx.fixef_vary_surv(obs.type)};

      auto &fixef_vary_basis_i = e_dat.fixef_vary_basis[i];
      for(vajoint_uint j = 0; j < fixef_vary_basis_i.n_cols(); ++j)
        for(vajoint_uint i = 0; i < fixef_vary_basis_i.n_rows(); ++i)
          outcome[j] += fixef_vary_basis_i.col(j)[i] * d_eta_i[i];
    }
  }

  // the derivatives for the association parameters
  {
    double const * d_rng_design_k{d_rng_design};
    for(vajoint_uint k = 0; k < n_markers();
        d_rng_design_k += rng_n_basis(k) * n_gl_outcomes, ++k){
      std::vector<std::vector<simple_mat<double> > > const &marker_k =
        e_dat.rng_basis[k];

      size_t const n_basis = rng_n_basis(k);

      double const * d_rng_design_kl{d_rng_design_k};
      for(size_t l = 0; l < info.size(); ++l, d_rng_design_kl += n_gl){
        auto &obs = info[l];
        auto const &rng_basis_kl = marker_k[l];
        double * outcome{gr + par_idx.association(obs.type)};

        // find the right association parameters
        auto &ders_l = ders_v[obs.type];
        for(size_t kk = 0; kk < k; ++kk)
          outcome += ders_l[kk].size();

        for(size_t v = 0; v < rng_basis_kl.size(); ++v){
          auto const &rng_basis_klv = rng_basis_kl[v];
          for(size_t j = 0; j < n_basis; ++j)
            for(size_t i = 0; i < n_gl; ++i)
              outcome[v] += d_rng_design_kl[i + j * n_gl_outcomes] *
                rng_basis_klv.col(j)[i];
        }
      }
    }
  }

  // the derivatives w.r.t. the covariance matrix
  double const * d_vcov_j{d_rng_design + n_rng * n_gl_outcomes};
  double * gr_vcov_vary{gr + par_idx.vcov_vary()};
  for(size_t j = 0; j < n_shared; ++j, d_vcov_j += n_rng,
      gr_vcov_vary += n_shared){
    for(vajoint_uint i = 0; i < n_shared; ++i)
      gr_vcov_vary[i] += d_vcov_j[i];
  }

  double * const gr_vcov_surv{gr + par_idx.vcov_surv()};
  auto const n_shared_surv = par_idx.n_shared_surv();
  size_t const n_left{e_dat.n_active_frailties()};

  for(size_t j = 0; j < n_left; ++j, d_vcov_j += n_rng){
    size_t const offset{e_dat.idx_inv_active_fraitly[j] * n_shared_surv};
    for(size_t i = 0; i < e_dat.idx_inv_active_fraitly.size(); ++i)
      gr_vcov_surv[e_dat.idx_inv_active_fraitly[i] + offset] +=
        d_vcov_j[i + n_shared];
  }

  return fn;
}
} // namespace survival
