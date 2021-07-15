#include "kl-term.h"
#include "lp-joint.h"

kl_term::kl_term(subset_params const &idx):
  idx(idx),
  n_vars(idx.get_n_shared() + idx.get_n_shared_survival()),
  vcov_inv     (idx.get_n_shared(), idx.get_n_shared()),
  vcov_inv_chol(idx.get_n_shared(), idx.get_n_shared()),
  vcov_survival_inv
    (idx.get_n_shared_survival(), idx.get_n_shared_survival()),
  vcov_survival_inv_chol
    (idx.get_n_shared_survival(), idx.get_n_shared_survival())
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

      local_eval_constant += (log_det - dim) / 2.;

      return true;
    };

  has_vcov = setup_vcov_decomps
    (param + idx.get_idx_shared_effect(), idx.get_n_shared(),
     vcov_inv, vcov_inv_chol);

  has_vcov_survival = setup_vcov_decomps
    (param + idx.get_idx_shared_survival(), idx.get_n_shared_survival(),
     vcov_survival_inv, vcov_survival_inv_chol);
}

double kl_term::eval(double const *param, double *wk_mem) {
  double out(eval_constant);
  if(!has_vcov or !has_vcov_survival)
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

  vajoint_uint const n_shared_survival = idx.get_n_shared_survival();
  if(has_vcov_survival){
    double term(0.);
    term += lp_joint::quad_form
      (va_mean + n_shared, vcov_survival_inv_chol.memptr(),
       n_shared_survival, wk_mem);

    term += lp_joint::submat_trace(vcov_survival_inv.memptr(), va_vcov,
                                   n_shared_survival, n_vars, n_shared);

    out += term / 2;
  }

  return out;
}

double kl_term::grad(double *g, double const *param, double *wk_mem) {
  if(!has_vcov or !has_vcov_survival)
    return eval(param, wk_mem);

  double const * const va_mean = param + idx.get_n_params(),
               * const va_vcov = va_mean + n_vars;



  return eval(param, wk_mem);
}
