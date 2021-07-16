#include "kl-term.h"
#include "lp-joint.h"

kl_term::kl_term(subset_params const &idx):
  idx(idx),
  n_vars(idx.get_n_shared() + idx.get_n_shared_surv()),
  vcov_inv     (idx.get_n_shared(), idx.get_n_shared()),
  vcov_inv_chol(idx.get_n_shared(), idx.get_n_shared()),
  vcov_surv_inv
    (idx.get_n_shared_surv(), idx.get_n_shared_surv()),
  vcov_surv_inv_chol
    (idx.get_n_shared_surv(), idx.get_n_shared_surv())
  { }

void kl_term::setup(double const *param, double *wk_mem){
  double &local_eval_constant = eval_constant;
  local_eval_constant = 0;

  auto setup_vcov_decomps =
    [&local_eval_constant]
    (double const *m, vajoint_uint const dim, arma::mat &vcov_inv,
     arma::mat &vcov_inv_chol){
      if(dim == 0)
        return false;

      arma::mat vcov(const_cast<double*>(m), dim, dim, false);
      bool has_vcov = false;
      for(double x : vcov)
        if(x != 0){
          has_vcov = true;
          break;
        }
      if(!has_vcov)
        return false;

      if(!arma::inv(vcov_inv, vcov))
        throw std::runtime_error("kl_term: inv() failed");
      if(!arma::chol(vcov_inv_chol, vcov_inv))
        throw std::runtime_error("kl_term: chol(inv()) failed");

      double log_det;
      if(!arma::log_det_sympd(log_det, vcov))
        throw std::runtime_error("kl_term: log_det_sympd(vcov) failed");

      local_eval_constant += (log_det - static_cast<double>(dim)) / 2.;

      return true;
    };

  has_vcov = setup_vcov_decomps
    (param + idx.get_idx_shared_effect(), idx.get_n_shared(),
     vcov_inv, vcov_inv_chol);

  has_vcov_surv = setup_vcov_decomps
    (param + idx.get_idx_shared_surv(), idx.get_n_shared_surv(),
     vcov_surv_inv, vcov_surv_inv_chol);
}

double kl_term::eval(double const *param, double *wk_mem) {
  double out(eval_constant);
  if(!has_vcov or !has_vcov_surv)
    return out;

  double const * const va_mean = param + idx.get_idx_va_mean(),
               * const va_vcov = param + idx.get_idx_va_vcov();

   {
      double log_det;
      arma::mat vcov_va_mat(const_cast<double*>(va_vcov), n_vars, n_vars,
                            false);
      if(!arma::log_det_sympd(log_det, vcov_va_mat))
        throw std::runtime_error("kl_term: log_det_sympd(vcov_va_mat) failed");
      out -= log_det / 2;
   }

  vajoint_uint const n_shared = idx.get_n_shared();
  if(has_vcov){
    double term(0.);
    term += lp_joint::quad_form(va_mean, vcov_inv_chol.memptr(),
                                n_shared, wk_mem);

    term += lp_joint::submat_trace(vcov_inv.memptr(), va_vcov,
                                   n_shared, n_vars, 0);

    out += term / 2;
  }

  vajoint_uint const n_shared_surv = idx.get_n_shared_surv();
  if(has_vcov_surv){
    double term(0.);
    term += lp_joint::quad_form
      (va_mean + n_shared, vcov_surv_inv_chol.memptr(),
       n_shared_surv, wk_mem);

    term += lp_joint::submat_trace(vcov_surv_inv.memptr(), va_vcov,
                                   n_shared_surv, n_vars, n_shared);

    out += term / 2;
  }

  return out;
}

double kl_term::grad(double *g, double const *param, double *wk_mem) {
  if(!has_vcov or !has_vcov_surv)
    return eval(param, wk_mem);

  double const * const va_mean = param + idx.get_idx_va_mean(),
               * const va_vcov = param + idx.get_idx_va_vcov();

  // handle some of the terms from the VA covariance matrix
  vajoint_uint const n_shared = idx.get_n_shared(),
                n_shared_surv = idx.get_n_shared_surv(),
                       n_vars = n_shared + n_shared_surv;

 {
    arma::mat va_cov_mat(const_cast<double*>(va_vcov),
                         n_vars, n_vars, false),
              inv_mat   (wk_mem, n_vars, n_vars, false);

    if(!arma::inv_sympd(inv_mat, va_cov_mat))
      throw std::runtime_error("inv(va_cov_mat) failed");

    lp_joint::add(g + idx.get_idx_va_vcov(), inv_mat.begin(), n_vars * n_vars,
                  -0.5);
 }

  auto add_vcov_term =
    [&]
    (vajoint_uint const dim, vajoint_uint const offset,
     vajoint_uint const idx_par,
     arma::mat const &inv_mat){
    arma::mat i1(wk_mem, dim, dim, false),
              i2(i1.end(), dim, dim, false);

    {
      double const *  vcov_surv = param + idx_par;
      double const * va_vcov_ele = va_vcov + offset * (offset + dim + 1);
      double * __restrict__ to = i1.memptr();
      for(vajoint_uint j = 0; j < dim; ++j, va_vcov_ele += n_vars - dim)
        for(vajoint_uint i = 0; i < dim; ++i)
          *to++ = *vcov_surv++ - *va_vcov_ele++;
    }

    double const * const va_mean_sub = va_mean + offset;
    lp_joint::rank_one<false>(i1.memptr(), va_mean_sub, dim);

    i2 = inv_mat * i1;
    i1 = i2 * inv_mat;

    lp_joint::add(g + idx_par, i1.begin(), dim * dim, .5);

    // handle the terms from the VA covariance matrix
    lp_joint::mat_add(g + idx.get_idx_va_vcov(), inv_mat.memptr(),
                      dim, n_vars, offset, .5);

    // handle the terms from the VA mean
    lp_joint::mat_vec_prod(
      inv_mat.memptr(), va_mean_sub, g + idx.get_idx_va_mean() + offset, dim);
  };

  if(has_vcov)
    add_vcov_term(n_shared, 0, idx.get_idx_shared_effect(), vcov_inv);
  if(has_vcov_surv)
    add_vcov_term(n_shared_surv, n_shared, idx.get_idx_shared_surv(),
                  vcov_surv_inv);

  return eval(param, wk_mem);
}
