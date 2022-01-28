#include "kl-term.h"
#include "lp-joint.h"
#include <algorithm>

inline void copy_sub_mat
  (double *wk_mem, double const *mat, vajoint_uint const offset,
   vajoint_uint const dim, vajoint_uint const n_vars){
  double const * ele{mat + offset * (offset + dim + 1)};
  for(vajoint_uint j = 0; j < dim; ++j, ele += n_vars, wk_mem += dim)
    std::copy(ele, ele + dim, wk_mem);
}

kl_term::kl_term(subset_params const &idx):
  idx(idx),
  n_vars(idx.n_shared() + idx.n_shared_surv())
  { }

void kl_term::setup(double const *param, double *wk_mem,
                    lb_terms which){
  eval_constant = 0;
  which_terms = which;

  has_vcov = idx.n_shared() &&
    (which_terms == lb_terms::all || which_terms == lb_terms::markers);
  if(has_vcov){
    vajoint_uint const dim{idx.n_shared()};
    vcov_fac.reset(new cfaad::CholFactorization
                     (param + idx.vcov_vary(), dim, true));

    eval_constant += log(vcov_fac->determinant()) - static_cast<double>(dim);
  }

  has_vcov_surv = idx.n_shared_surv() &&
    (which_terms == lb_terms::all || which_terms == lb_terms::surv);
  if(has_vcov_surv){
    vajoint_uint const dim{idx.n_shared_surv()};
    vcov_surv_fac.reset(new cfaad::CholFactorization(
      param + idx.vcov_surv(), dim, true));

    eval_constant +=
      log(vcov_surv_fac->determinant()) - static_cast<double>(dim);
  }
}

double kl_term::eval(double const *param, double *wk_mem) const {
  double out(eval_constant);
  if(!has_vcov and !has_vcov_surv)
    return out;

  vajoint_uint const n_shared{idx.n_shared()},
                n_shared_surv{idx.n_shared_surv()};

  double const * const va_mean{param + idx.va_mean()},
               * const va_vcov{param + idx.va_vcov()};

    if(which_terms == lb_terms::markers && n_shared_surv){
      copy_sub_mat(wk_mem, va_vcov, 0, n_shared, n_vars);

      double log_det;
      arma::mat vcov_va_mat(const_cast<double*>(wk_mem), n_shared, n_shared,
                            false);
      if(!arma::log_det_sympd(log_det, vcov_va_mat))
        throw std::runtime_error("kl_term: log_det_sympd(vcov_va_mat) failed");
      out -= log_det;

    }
    else if(which_terms == lb_terms::surv && n_shared){
      copy_sub_mat(wk_mem, va_vcov, n_shared, n_shared_surv, n_vars);

      double log_det;
      arma::mat vcov_va_mat(const_cast<double*>(wk_mem), n_shared_surv,
                            n_shared_surv, false);
      if(!arma::log_det_sympd(log_det, vcov_va_mat))
        throw std::runtime_error("kl_term: log_det_sympd(vcov_va_mat) failed");
      out -= log_det;

    } else {
      double log_det;
      arma::mat vcov_va_mat(const_cast<double*>(va_vcov), n_vars, n_vars,
                            false);
      if(!arma::log_det_sympd(log_det, vcov_va_mat))
        throw std::runtime_error("kl_term: log_det_sympd(vcov_va_mat) failed");
      out -= log_det;

    }

  // handle the non-determinant term
  auto handle_terms = [&]
  (vajoint_uint const offset, vajoint_uint const dim,
   cfaad::CholFactorization const &fact, double const *org){
    // copy parts of the VA covariance matrix
    double * const sub_VA_vcov{wk_mem};
    copy_sub_mat(sub_VA_vcov, va_vcov, offset, dim, n_vars);

    // handle the two terms
    double term{};
    term += cfaad::quadFormInv(va_mean + offset, sub_VA_vcov, fact);
    term += cfaad::trInvMatMat(org, sub_VA_vcov, fact);

    return term;
  };

  if(has_vcov)
    out += handle_terms(0, n_shared, *vcov_fac, param + idx.vcov_vary());
  if(has_vcov_surv)
    out += handle_terms
    (n_shared, n_shared_surv, *vcov_surv_fac, param + idx.vcov_surv());

  return out / 2;
}

double kl_term::grad(double *g, double const *param, double *wk_mem) const {
  if(!has_vcov and !has_vcov_surv)
    return eval(param, wk_mem);

  double const * const va_mean = param + idx.va_mean(),
               * const va_vcov = param + idx.va_vcov();

  // handle some of the terms from the VA covariance matrix
  vajoint_uint const n_shared = idx.n_shared(),
                n_shared_surv = idx.n_shared_surv();

  auto handle_sub_part = [&](vajoint_uint const offset, vajoint_uint const dim){
    copy_sub_mat(wk_mem, va_vcov, offset, dim, n_vars);
    arma::mat va_cov_mat(const_cast<double*>(wk_mem),
                         dim, dim, false),
              inv_mat   (wk_mem + dim * dim, dim, dim, false);

    if(!arma::inv_sympd(inv_mat, va_cov_mat))
      throw std::runtime_error("inv(va_cov_mat) failed");

    double *der_term{g + idx.va_vcov() + offset * (offset + dim + 1)};
    double const *term{inv_mat.begin()};

    for(vajoint_uint i = 0; i < dim; ++i, der_term += n_vars, term += dim)
      for(vajoint_uint j = 0; j < dim; ++j)
        der_term[j] -= .5 * term[j];
  };

  if(which_terms == lb_terms::markers && n_shared_surv)
    handle_sub_part(0, n_shared);
  else if(which_terms == lb_terms::surv && n_shared)
    handle_sub_part(n_shared, n_shared_surv);
  else {
    arma::mat va_cov_mat(const_cast<double*>(va_vcov),
                         n_vars, n_vars, false),
              inv_mat   (wk_mem, n_vars, n_vars, false);

    if(!arma::inv_sympd(inv_mat, va_cov_mat))
      throw std::runtime_error("inv(va_cov_mat) failed");

    lp_joint::add(g + idx.va_vcov(), inv_mat.begin(), n_vars * n_vars,
                  -0.5);
 }

  auto add_vcov_term =
    [&]
    (vajoint_uint const dim, vajoint_uint const offset,
     vajoint_uint const idx_par, cfaad::CholFactorization const &fact){
    arma::mat i1(wk_mem, dim, dim, false),
              i2(i1.end(), dim, dim, false);

    {
      double const *  vcov_ele = param + idx_par;
      double const * va_vcov_ele = va_vcov + offset * (offset + dim + 1);
      double * __restrict__ to = i1.memptr();
      for(vajoint_uint j = 0; j < dim; ++j, va_vcov_ele += n_vars - dim)
        for(vajoint_uint i = 0; i < dim; ++i)
          *to++ = *vcov_ele++ - *va_vcov_ele++;
    }

    double const * const va_mean_sub = va_mean + offset;
    lp_joint::rank_one<false>(i1.memptr(), va_mean_sub, dim);

    for(vajoint_uint i = 0; i < dim; ++i)
      fact.solve(i1.memptr() + i * dim);
    i2 = i1.t();
    for(vajoint_uint i = 0; i < dim; ++i)
      fact.solve(i2.memptr() + i * dim);

    lp_joint::add(g + idx_par, i2.begin(), dim * dim, .5);

    // copy the upper triangular matrix to the full matrix
    double * const fact_inv{wk_mem};
    {
      double const * v_inv_ele(fact.get_inv());
      for(vajoint_uint j = 0; j < dim; ++j){
        for(vajoint_uint i = 0; i < j; ++i, ++v_inv_ele){
          fact_inv[i + j * dim] = *v_inv_ele;
          fact_inv[j + i * dim] = *v_inv_ele;
        }
        fact_inv[j + j * dim] = *v_inv_ele++;
      }
    }

    // handle the terms from the VA covariance matrix
    lp_joint::mat_add(g + idx.va_vcov(), fact_inv,
                      dim, n_vars, offset, .5);

    // handle the terms from the VA mean
    lp_joint::mat_vec_prod
      (fact_inv, va_mean_sub, g + idx.va_mean() + offset, dim);
  };

  if(has_vcov)
    add_vcov_term(n_shared, 0, idx.vcov_vary(), *vcov_fac);
  if(has_vcov_surv)
    add_vcov_term(n_shared_surv, n_shared, idx.vcov_surv(), *vcov_surv_fac);

  return eval(param, wk_mem);
}
