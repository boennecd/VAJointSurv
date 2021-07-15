#include <testthat.h>
#include "kl-term.h"
#include <memory.h>

context("testing kl-terms") {
  test_that("eval gives the right result") {
    /*
     set.seed(1)
     n_shared <- 2L
     n_shared_survival <- 3L
     n_vars <- n_shared + n_shared_survival

     Omega <- drop(rWishart(1, n_vars, diag(n_vars)))
     Xi <- drop(rWishart(1, n_shared_survival, diag(n_shared_survival)))
     Psi <- drop(rWishart(1, n_shared, diag(n_shared)))
     zeta <- rnorm(n_vars)

     val <- -determinant(Omega)$modulus + determinant(Xi)$modulus +
     determinant(Psi)$modulus +
     drop(zeta[1:n_shared] %*% solve(Psi, zeta[1:n_shared])) +
     drop(zeta[-(1:n_shared)] %*% solve(Xi, zeta[-(1:n_shared)])) +
     sum(diag(solve(Psi, Omega[1:n_shared, 1:n_shared]))) +
     sum(diag(solve(Xi, Omega[-(1:n_shared), -(1:n_shared)]))) -
     n_shared - n_shared_survival
     dput(val / 2)
     dput(Xi)
     dput(Psi)
     dput(Omega)
     dput(zeta)
     */
    constexpr vajoint_uint n_shared = 2,
                  n_shared_survival = 3,
                             n_vars = n_shared + n_shared_survival;
    constexpr double Xi[n_shared_survival * n_shared_survival] = { 1.96775053611171, -1.73597597741474, 0.529397523176239, -1.73597597741474, 3.24256054526995, -0.292627703276501, 0.529397523176239, -0.292627703276501, 0.634396281932773 },
                    Psi[n_shared * n_shared] = { 2.4606560951913, 0.789983565757713, 0.789983565757713, 0.892097273439034},
                  Omega[n_vars * n_vars] = { 2.42434323779257, 1.9812109601339, -2.3977488177111, 0.896508989006271, -0.967290384087283, 1.9812109601339, 8.7605890723572, -4.44094380859342, -0.0834669056878007, -6.70896207863171, -2.3977488177111, -4.44094380859342, 6.14892949801278, 1.97812834810877, 4.9338943130402, 0.896508989006271, -0.0834669056878007, 1.97812834810877, 3.33690095112284, 1.98372476564407, -0.967290384087283, -6.70896207863171, 4.9338943130402, 1.98372476564407, 7.74887957345459 },
                  zeta[n_vars] = { 1.08576936214569, -0.69095383969683, -1.28459935387219, 0.046726172188352, -0.235706556439501 };

    subset_params params;
    params.add_marker({ 2, 2, 1 });
    params.add_marker({ 1, 4, 1 });
    params.add_survival({ 5, 2 });
    params.add_survival({ 1, 4 });
    params.add_survival({ 5, 2 });

    // create and fill parameter vector
    vajoint_uint const n_params_w_va = params.get_n_parms_w_va();
    std::unique_ptr<double[]> par(new double[n_params_w_va]);
    std::fill(par.get(), par.get() + n_params_w_va, 0.);

    auto fill_par = [](double const *value, vajoint_uint const n_ele,
                       double *out){
      std::copy(value, value + n_ele, out);
    };
    fill_par(Xi, n_shared_survival * n_shared_survival,
             par.get() + params.get_idx_shared_survival());
    fill_par(Psi, n_shared * n_shared,
             par.get() + params.get_idx_shared_effect());
    fill_par(Omega, n_vars * n_vars,
             par.get() + params.get_idx_va_vcov());
    fill_par(zeta, n_vars,
             par.get() + params.get_idx_va_mean());

    // compute the kl term
    {
      kl_term term(params);
      std::unique_ptr<double[]> mem(new double[term.get_n_dmen()]);
      term.setup(par.get(), mem.get());
      constexpr double truth = 14.58945197638;
      expect_true(
        std::abs(term.eval(par.get()) - truth) < std::abs(truth) * 1e-8);

      // also works when provided with the working memory
      expect_true(std::abs(term.eval(par.get(), mem.get()) - truth) <
        std::abs(truth) * 1e-8);
    }
  }
}
