#include "testthat-wrapper.h"
#include "kl-term.h"
#include <memory>
#include "log-cholesky.h"

namespace {
/*
 set.seed(1)
 n_shared <- 2L
 n_shared_surv <- 3L
 n_vars <- n_shared + n_shared_surv

 Omega <- drop(rWishart(1, n_vars, diag(n_vars)))
 Xi <- drop(rWishart(1, n_shared_surv, diag(n_shared_surv)))
 Psi <- drop(rWishart(1, n_shared, diag(n_shared)))
 zeta <- rnorm(n_vars)

 f <- function(Omega, Xi, Psi, zeta)
 (-determinant(Omega)$modulus + determinant(Xi)$modulus +
 determinant(Psi)$modulus +
 drop(zeta[1:n_shared] %*% solve(Psi, zeta[1:n_shared])) +
 drop(zeta[-(1:n_shared)] %*% solve(Xi, zeta[-(1:n_shared)])) +
 sum(diag(solve(Psi, Omega[1:n_shared, 1:n_shared]))) +
 sum(diag(solve(Xi, Omega[-(1:n_shared), -(1:n_shared)]))) -
 n_shared - n_shared_surv)/2

 dput(f(Omega, Xi, Psi, zeta))
 dput(Xi)
 dput(Psi)
 dput(Omega)
 dput(zeta)

 log_chol <- function(x){
 x <- chol(x)
 diag(x) <- log(diag(x))
 x[upper.tri(x, TRUE)]
 }
 log_chol_inv <- function(x){
 n <- (sqrt(8 * length(x) + 1) - 1) / 2
 out <- matrix(0, n, n)
 out[upper.tri(out, TRUE)] <- x
 diag(out) <- exp(diag(out))
 crossprod(out)
 }
 dput(Xi_chol <- log_chol(Xi))
 dput(Psi_chol <- log_chol(Psi))
 dput(Omega_chol <- log_chol(Omega))

 g <- function(x){
 Xi_chol <- x[seq_along(Xi_chol)]
 Psi_chol <- x[seq_along(Psi_chol) + length(Xi_chol)]
 Omega_chol <- x[seq_along(Omega_chol) + length(Xi_chol) + length(Psi_chol)]
 zeta <- tail(x, length(zeta))
 f(Omega = log_chol_inv(Omega_chol), Psi = log_chol_inv(Psi_chol),
 Xi = log_chol_inv(Xi_chol), zeta = zeta)
 }

 g(c(Xi_chol, Psi_chol, Omega_chol, zeta))
 deriv <- numDeriv::grad(g, c(Xi_chol, Psi_chol, Omega_chol, zeta))
 dput(Xi_deriv <- deriv[seq_along(Xi_chol)])
 dput(Psi_deriv <- deriv[seq_along(Psi_chol) + length(Xi_chol)])
 dput(Omega_deriv <- deriv[
 seq_along(Omega_chol) + length(Xi_chol) + length(Psi_chol)])
 dput(zeta_deriv <- tail(deriv, length(zeta)))

# only with the markers
 f <- function(Omega, Xi, Psi, zeta)
 (-determinant(Omega[1:n_shared, 1:n_shared])$modulus + # determinant(Xi)$modulus +
 determinant(Psi)$modulus +
 drop(zeta[1:n_shared] %*% solve(Psi, zeta[1:n_shared])) +
# drop(zeta[-(1:n_shared)] %*% solve(Xi, zeta[-(1:n_shared)])) +
 sum(diag(solve(Psi, Omega[1:n_shared, 1:n_shared]))) +
# sum(diag(solve(Xi, Omega[-(1:n_shared), -(1:n_shared)])))
 -n_shared)/2

 dput(f(Omega, Xi, Psi, zeta))
 deriv <- numDeriv::grad(g, c(Xi_chol, Psi_chol, Omega_chol, zeta))

 dput(Xi_deriv <- deriv[seq_along(Xi_chol)])
 dput(Psi_deriv <- deriv[seq_along(Psi_chol) + length(Xi_chol)])
 dput(Omega_deriv <- deriv[
 seq_along(Omega_chol) + length(Xi_chol) + length(Psi_chol)])
 dput(zeta_deriv <- tail(deriv, length(zeta)))

# only survival
 f <- function(Omega, Xi, Psi, zeta)
 (-determinant(Omega[-(1:n_shared), -(1:n_shared)])$modulus +
 determinant(Xi)$modulus +
# determinant(Psi)$modulus +
# drop(zeta[1:n_shared] %*% solve(Psi, zeta[1:n_shared])) +
 drop(zeta[-(1:n_shared)] %*% solve(Xi, zeta[-(1:n_shared)])) +
# sum(diag(solve(Psi, Omega[1:n_shared, 1:n_shared]))) +
 sum(diag(solve(Xi, Omega[-(1:n_shared), -(1:n_shared)]))) -
 n_shared_surv)/2

 dput(f(Omega, Xi, Psi, zeta))
 deriv <- numDeriv::grad(g, c(Xi_chol, Psi_chol, Omega_chol, zeta))

 dput(Xi_deriv <- deriv[seq_along(Xi_chol)])
 dput(Psi_deriv <- deriv[seq_along(Psi_chol) + length(Xi_chol)])
 dput(Omega_deriv <- deriv[
 seq_along(Omega_chol) + length(Xi_chol) + length(Psi_chol)])
 dput(zeta_deriv <- tail(deriv, length(zeta)))
 */


constexpr vajoint_uint n_shared = 2,
              n_shared_surv = 3,
                         n_vars = n_shared + n_shared_surv;
constexpr double Xi[n_shared_surv * n_shared_surv] { 1.96775053611171, -1.73597597741474, 0.529397523176239, -1.73597597741474, 3.24256054526995, -0.292627703276501, 0.529397523176239, -0.292627703276501, 0.634396281932773 },
                Psi[n_shared * n_shared] { 2.4606560951913, 0.789983565757713, 0.789983565757713, 0.892097273439034},
              Omega[n_vars * n_vars] { 2.42434323779257, 1.9812109601339, -2.3977488177111, 0.896508989006271, -0.967290384087283, 1.9812109601339, 8.7605890723572, -4.44094380859342, -0.0834669056878007, -6.70896207863171, -2.3977488177111, -4.44094380859342, 6.14892949801278, 1.97812834810877, 4.9338943130402, 0.896508989006271, -0.0834669056878007, 1.97812834810877, 3.33690095112284, 1.98372476564407, -0.967290384087283, -6.70896207863171, 4.9338943130402, 1.98372476564407, 7.74887957345459 },
              zeta[n_vars] { 1.08576936214569, -0.69095383969683, -1.28459935387219, 0.046726172188352, -0.235706556439501 },
           Xi_chol[dim_tri(n_shared_surv)] { 0.338445515244742, -1.23753842192996, 0.268556296839291, 0.377395645981701, 0.133336360814841, -0.373073361500074 },
          Psi_chol[dim_tri(n_shared)] { 0.450214009873517, 0.503607972233726, -0.224335373954299 },
        Omega_chol[dim_tri(n_vars)] { 0.442780328966089, 1.2724293214294, 0.982962307931023, -1.53995004190371, -0.928567034713538, 0.534977211455984, 0.575781351653492, -0.305388387156356, 1.51178116845085, -0.23369758278885, -0.621240580541804, -2.2146998871775, 1.12493091814311, -0.0449336090152294, 0.0872100127375849 };

constexpr double true_kl_term = 14.58945197638;

constexpr double Xi_deriv[dim_tri(n_shared_surv)] { -5.80311449958224, -3.32225269517676, -5.84896643091192, -3.36986740874183,
                                                    -4.70975985678572, -9.17044789534625 },
                Psi_deriv[dim_tri(n_shared)] { -0.427035663737666, -0.0740751249053278, -12.8125448527808 },
              Omega_deriv[dim_tri(n_vars)] { -0.61960594919438, 1.20999189744859, 10.1852434869524, -1.08994359451185,
                                             0.348923601855757, 2.67579154568366, -0.430153063324311, -0.259508809540019,
                                             1.68692480260345, -0.612398324838771, -0.268131842823628, -3.90191228270156,
                                             0.754827608496878, -0.264923488296771, 1.51071293143399 },
              zeta_deriv[n_vars] { 0.963963111992338, -1.62815076284324, -1.38009039095917, -0.682457259304228,
                                   0.465330561405237 };

// without the survival part
constexpr double true_kl_term_no_surv{5.43857866944945},
                 Xi_deriv_no_s[] = {0., 0., 0., 0., 0., 0.},
                 Psi_deriv_no_s[] = {-0.427035663832965, -0.0740751251107456, -12.8125448527058},
                 Omega_deriv_no_s[] = {-0.619605949523746, 1.20999189732096, 10.1852434869334, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                 zeta_deriv_no_s[] = {0.963963112174406, -1.62815076286685, 0, 0, 0};

// without markers
constexpr double true_kl_term_no_marker{7.47606508560383},
                 Xi_deriv_no_m[] = {-5.803114499632, -3.32225269517634, -5.84896643129221, -3.36986740815176,
                                    -4.70975985768593, -9.17044789529971},
                 Psi_deriv_no_m[] = {0, 0, 0},
                 Omega_deriv_no_m[] = {0, 0, 0, -0.619384955968862, 0.201941440223768, 3.26394555054402,
                                       -0.808825953840981, -0.314820162313017, 1.31531055565136, 0.15054243930918,
                                       -0.39063512977232, -3.50835643881455, 0.858372808603325, -0.221739103733836,
                                       2.19167623243605},
                 zeta_deriv_no_m[] = {0, 0, -1.3800903909731, -0.682457259664746, 0.46533056191001};

} // namespace

context("testing kl-terms") {
  test_that("eval gives the right result") {
    double const eps = std::sqrt(std::numeric_limits<double>::epsilon());

    subset_params params;
    params.add_marker({ 2, 2, 1 });
    params.add_marker({ 1, 4, 1 });
    params.add_surv({ 5, 2, {1, 1}, true });
    params.add_surv({ 1, 4, {1, 1}, true });
    params.add_surv({ 5, 2, {1, 1}, true });

    // create and fill parameter vector
    vajoint_uint const n_params_w_va = params.n_params_w_va();
    std::unique_ptr<double[]> par(new double[n_params_w_va]);
    std::fill(par.get(), par.get() + n_params_w_va, 0.);

    auto fill_par = [](double const *value, vajoint_uint const n_ele,
                       double *out){
      std::copy(value, value + n_ele, out);
    };
    fill_par(Xi, n_shared_surv * n_shared_surv,
             par.get() + params.vcov_surv());
    fill_par(Psi, n_shared * n_shared,
             par.get() + params.vcov_vary());
    fill_par(Omega, n_vars * n_vars,
             par.get() + params.va_vcov());
    fill_par(zeta, n_vars,
             par.get() + params.va_mean());

    // compute the kl term
    kl_term term(params);
    std::unique_ptr<double[]> mem(new double[term.n_wmem()]);
    term.setup(par.get(), mem.get());
    expect_true(pass_rel_err(term.eval(par.get(), mem.get()), true_kl_term));

    // check the gradient
    std::unique_ptr<double[]> gr(new double[n_params_w_va]),
                          gr_res(new double[params.n_params_w_va<true>()]);
    std::fill(gr.get(), gr.get() + n_params_w_va, 0.);
    std::fill(gr_res.get(), gr_res.get() + params.n_params_w_va<true>(), 0.);

    double val = term.grad(gr.get(), par.get(), mem.get());
    expect_true(pass_rel_err(val, true_kl_term));

    auto test_grad = [&](double const *Xi_deriv, double const *Psi_deriv,
                         double const *Omega_deriv, double const *zeta_deriv){
      {
        double *g_out = gr_res.get() + params.vcov_surv<true>();
        double const *g_in = gr.get() + params.vcov_surv();
        log_chol::dpd_mat::get(Xi_chol, n_shared_surv, g_out, g_in);

        for(vajoint_uint i = 0; i < dim_tri(n_shared_surv); ++i)
          expect_true(pass_rel_err(g_out[i], Xi_deriv[i]));
      }

      {
        double *g_out = gr_res.get() + params.vcov_vary<true>();
        double const *g_in = gr.get() + params.vcov_vary();
        log_chol::dpd_mat::get(Psi_chol, n_shared, g_out, g_in);

        for(vajoint_uint i = 0; i < dim_tri(n_shared); ++i)
          expect_true(pass_rel_err(g_out[i], Psi_deriv[i]));
      }

      {
        double *g_out = gr_res.get() + params.va_vcov<true>();
        double const *g_in = gr.get() + params.va_vcov();
        log_chol::dpd_mat::get(Omega_chol, n_vars, g_out, g_in);

        for(vajoint_uint i = 0; i < dim_tri(n_vars); ++i)
          expect_true(pass_rel_err(g_out[i], Omega_deriv[i]));
      }

      {
        double const *g_out = gr.get() + params.va_mean();
        for(vajoint_uint i = 0; i < n_vars; ++i)
          expect_true(pass_rel_err(g_out[i], zeta_deriv[i]));
      }
    };

    test_grad(Xi_deriv, Psi_deriv, Omega_deriv, zeta_deriv);

    term.setup(par.get(), mem.get(), lb_terms::markers);
    expect_true(
      pass_rel_err(term.eval(par.get(), mem.get()), true_kl_term_no_surv));

    std::fill(gr.get(), gr.get() + n_params_w_va, 0.);
    std::fill(gr_res.get(), gr_res.get() + params.n_params_w_va<true>(), 0.);

    val = term.grad(gr.get(), par.get(), mem.get());
    expect_true(pass_rel_err(val, true_kl_term_no_surv));

    test_grad(Xi_deriv_no_s, Psi_deriv_no_s, Omega_deriv_no_s, zeta_deriv_no_s);

    term.setup(par.get(), mem.get(), lb_terms::surv);
    expect_true(
      pass_rel_err(term.eval(par.get(), mem.get()), true_kl_term_no_marker));


    // clean up
    wmem::clear_all();
  }

  test_that("eval gives the right result with survival terms without frailty") {
    double const eps = std::sqrt(std::numeric_limits<double>::epsilon());

    subset_params params;
    params.add_marker({ 2, 2, 1 });
    params.add_marker({ 1, 4, 1 });
    params.add_surv({ 1, 2, {0, 0}, false });
    params.add_surv({ 5, 2, {1, 1}, true });
    params.add_surv({ 1, 4, {1, 1}, true });
    params.add_surv({ 5, 2, {1, 1}, true });
    params.add_surv({ 3, 4, {0, 0}, false });

    // create and fill parameter vector
    vajoint_uint const n_params_w_va = params.n_params_w_va();
    std::unique_ptr<double[]> par(new double[n_params_w_va]);
    std::fill(par.get(), par.get() + n_params_w_va, 0.);

    auto fill_par = [](double const *value, vajoint_uint const n_ele,
                       double *out){
      std::copy(value, value + n_ele, out);
    };
    fill_par(Xi, n_shared_surv * n_shared_surv,
             par.get() + params.vcov_surv());
    fill_par(Psi, n_shared * n_shared,
             par.get() + params.vcov_vary());
    fill_par(Omega, n_vars * n_vars,
             par.get() + params.va_vcov());
    fill_par(zeta, n_vars,
             par.get() + params.va_mean());

    // compute the kl term
    kl_term term(params);
    std::unique_ptr<double[]> mem(new double[term.n_wmem()]);
    term.setup(par.get(), mem.get());
    expect_true(pass_rel_err(term.eval(par.get(), mem.get()), true_kl_term));

    // check the gradient
    std::unique_ptr<double[]> gr(new double[n_params_w_va]),
                          gr_res(new double[params.n_params_w_va<true>()]);
    std::fill(gr.get(), gr.get() + n_params_w_va, 0.);
    std::fill(gr_res.get(), gr_res.get() + params.n_params_w_va<true>(), 0.);

    double val = term.grad(gr.get(), par.get(), mem.get());
    expect_true(pass_rel_err(val, true_kl_term));

    auto test_grad = [&](double const *Xi_deriv, double const *Psi_deriv,
                         double const *Omega_deriv, double const *zeta_deriv){
      {
        double *g_out = gr_res.get() + params.vcov_surv<true>();
        double const *g_in = gr.get() + params.vcov_surv();
        log_chol::dpd_mat::get(Xi_chol, n_shared_surv, g_out, g_in);

        for(vajoint_uint i = 0; i < dim_tri(n_shared_surv); ++i)
          expect_true(pass_rel_err(g_out[i], Xi_deriv[i]));
      }

      {
        double *g_out = gr_res.get() + params.vcov_vary<true>();
        double const *g_in = gr.get() + params.vcov_vary();
        log_chol::dpd_mat::get(Psi_chol, n_shared, g_out, g_in);

        for(vajoint_uint i = 0; i < dim_tri(n_shared); ++i)
          expect_true(pass_rel_err(g_out[i], Psi_deriv[i]));
      }

      {
        double *g_out = gr_res.get() + params.va_vcov<true>();
        double const *g_in = gr.get() + params.va_vcov();
        log_chol::dpd_mat::get(Omega_chol, n_vars, g_out, g_in);

        for(vajoint_uint i = 0; i < dim_tri(n_vars); ++i)
          expect_true(pass_rel_err(g_out[i], Omega_deriv[i]));
      }

      {
        double const *g_out = gr.get() + params.va_mean();
        for(vajoint_uint i = 0; i < n_vars; ++i)
          expect_true(pass_rel_err(g_out[i], zeta_deriv[i]));
      }
    };

    test_grad(Xi_deriv, Psi_deriv, Omega_deriv, zeta_deriv);

    term.setup(par.get(), mem.get(), lb_terms::markers);
    expect_true(
      pass_rel_err(term.eval(par.get(), mem.get()), true_kl_term_no_surv));

    std::fill(gr.get(), gr.get() + n_params_w_va, 0.);
    std::fill(gr_res.get(), gr_res.get() + params.n_params_w_va<true>(), 0.);

    val = term.grad(gr.get(), par.get(), mem.get());
    expect_true(pass_rel_err(val, true_kl_term_no_surv));

    test_grad(Xi_deriv_no_s, Psi_deriv_no_s, Omega_deriv_no_s, zeta_deriv_no_s);

    term.setup(par.get(), mem.get(), lb_terms::surv);
    expect_true(
      pass_rel_err(term.eval(par.get(), mem.get()), true_kl_term_no_marker));


    // clean up
    wmem::clear_all();
  }
}
