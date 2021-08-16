#include "marker-term.h"
#include "testthat-wrapper.h"
#include "wmem.h"
#include <iterator>

using std::begin;
using std::end;

context("marker_term is correct") {
  test_that("marker_term gives the correct result in the univariate case"){
    /* R code to reproduce the result
     raw_poly <- function(x, degree, intercept){
     if(intercept)
     drop(outer(x, 0:degree, `^`))
     else
     drop(outer(x, 1:degree, `^`))
     }

     # parameters
     set.seed(1)
     X <- matrix(c(1, -1, -2, 1, 2, .5), 3)
     ti <- c(1, 3)
     Sig <- .5^2
     gamma <- c(-1, -.25, -.4)
     beta <- .4
     Psi <- matrix(c(3, 1, .5, 1, 2, .25, .5, .25, 1), 3)
     dput(zeta <- drop(round(mvtnorm::rmvnorm(1, sigma = Psi), 2)))

     g1 <- function(x) raw_poly(x, 1, FALSE)
     m1 <- function(x) raw_poly(x, 2, TRUE)

     y <- sapply(1:2, function(i)
     rnorm(1, X[, i] %*% gamma + g1(ti[i]) %*% beta + m1(ti[i]) %*% zeta,
     sqrt(Sig)))
     dput(y <- round(y, 2))

     f <- function(x){
     # get the parameters
     get_next <- function(n){
     out <- head(x, n)
     x <<- tail(x, -n)
     out
     }
     g <- get_next(length(gamma))
     b <- get_next(length(beta))
     S <- get_next(1)
     z <- get_next(length(zeta))
     P <- matrix(x, NROW(Psi))

     # compute the output
     out <- 0
     for(i in 1:2){
     G <- g1(ti[i])
     M <- m1(ti[i])
     out <- out -
     dnorm(y[i], X[, i] %*% g + G %*% b + M %*% z, sqrt(S), log = TRUE) +
     M %*% P %*% M / S / 2
     }

     out
     }

     dput(f(c(gamma, beta, Sig, zeta, Psi)))
     dput(round(numDeriv::grad(f, c(gamma, beta, Sig, zeta, Psi)), 7))
     */

    // the bases
    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;
    // raw poly of degree one without an intercept
    bases_fix.emplace_back(new joint_bases::orth_poly(1, false));
    // raw poly of degree two with an intercept
    bases_rng.emplace_back(new joint_bases::orth_poly(2, true));

    constexpr vajoint_uint n_obs{2},
                         n_fixef{3};
    constexpr double obs[n_obs] { -0.84, -9.83 },
                obs_time[n_obs] {1, 3};
    constexpr int ids[n_obs] = {1, 1};

    // the design matrix
    double X[n_obs * n_fixef] {1, -1, -2, 1, 2, .5};
    constexpr double Sig[]{.25},
                     Psi[]{3, 1, .5, 1, 2, .25, .5, .25, 1},
                       g[]{-1, -.25, -.4},
                    beta[]{.4},
                    zeta[]{-1.15, -0.02, -0.92},
                  true_val{281.782782705289},
              // has to change if the order of the parameters change
             true_derivs[]{-3.84, 1.92, 6.08, -5.12, -1121.3248, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3.84, -5.1199999, -8.96, 4, 8, 20, 8, 20, 56, 20, 56, 164};

    // create the input data
    std::vector<marker::setup_marker_dat_helper> input_dat;
    input_dat.emplace_back(X, 3, n_obs, ids, obs_time, obs);

    subset_params par_idx;
    par_idx.add_marker({n_fixef, 1, 3});

    marker::marker_dat comp_obj = marker::get_comp_dat
      (input_dat, par_idx, bases_fix, bases_rng);

    // test the basic members
    expect_true(comp_obj.n_obs == n_obs);
    expect_true(comp_obj.n_markers == 1);

    // the function returns the right value
    std::vector<double> par(par_idx.get_n_parms_w_va());

    std::fill(par.begin(), par.end(), 0);
    std::copy(begin(g), end(g), par.begin() + par_idx.get_fixef_idx_marker(0));
    std::copy(begin(beta), end(beta),
              par.begin() + par_idx.get_varying_idx_marker(0));
    std::copy(begin(Sig), end(Sig), par.begin() + par_idx.get_idx_error_term());
    std::copy(begin(zeta), end(zeta), par.begin() + par_idx.get_idx_va_mean());
    std::copy(begin(Psi), end(Psi), par.begin() + par_idx.get_idx_va_vcov());

    {
      double *wk_mem = wmem::get_double_mem(comp_obj.get_n_wmem());
      comp_obj.setup(par.data(), wk_mem);

      double res = comp_obj.eval(par.data(), wk_mem, 0)
        + comp_obj.eval(par.data(), wk_mem, 1);

      expect_true(res == Approx(true_val).epsilon(1e-7));
    }

    // test the gradient
    cfaad::Number::tape->rewind();
    std::vector<cfaad::Number> ad_par(par.size());
    cfaad::convertCollection(par.begin(), par.end(), ad_par.begin());

    cfaad::Number * wk_mem = wmem::get_Number_mem(comp_obj.get_n_wmem());
    cfaad::Number res = comp_obj.eval(ad_par.data(), wk_mem, 0) +
      comp_obj.eval(ad_par.data(), wk_mem, 1);

    res.propagateToStart();
    expect_true(res.value() == Approx(true_val).epsilon(1e-6));

    expect_true(ad_par.size() == static_cast<size_t>(
      std::distance(begin(true_derivs), end(true_derivs))));

    for(size_t i = 0; i < ad_par.size(); ++i)
      expect_true(ad_par[i].adjoint() == Approx(true_derivs[i]).epsilon(1e-6));

    // clean up
    wmem::clear_all();
  }

  test_that("marker_term gives the correct result with three markers"){
    /* R code to reproduce the result
     raw_poly <- function(x, degree, intercept){
     if(intercept)
     drop(outer(x, 0:degree, `^`))
     else
     drop(outer(x, 1:degree, `^`))
     }

# parameters
     set.seed(1)
     Xs <- list(matrix(c(1:2, -1, -2), 2),
     matrix(-(1:2), 2),
     matrix(c(.5, .67), 1))
     gs <- list(c(-.25, .33), c(-1, 2), c(1.5))
     betas <- list(c(.2), c(.5, 3), c(2.5))
     dput(Sig <- round(drop(rWishart(1, 6, diag(3))), 2))
     stopifnot(all(eigen(Sig)$values > 0))
     dput(Psi <- round(drop(rWishart(1, 21, diag(7))), 2))
     dput(zeta <- drop(round(mvtnorm::rmvnorm(1, sigma = Psi), 2)))
     rng_idx <- list(1:3, 4:5, 6:7)

     g_funcs <- list(function(x) raw_poly(x, 1, FALSE),
     function(x) raw_poly(x, 2, FALSE),
     function(x) raw_poly(x, 1, FALSE))
     m_funcs <- list(function(x) raw_poly(x, 2, TRUE),
     function(x) raw_poly(x, 1, TRUE),
     function(x) raw_poly(x, 1, TRUE))

     obs <- list(
     list(is_obs = 1:3, ti = .5, idx = c(1L, 1L, 1L)),
     list(is_obs = c(1L, 3L), ti = 1.2, idx = c(2L, 2L)))

# generate a plausible
     obs <- lapply(obs, function(x){
# compute the mean
     ti <- x$ti
     mea <- mapply(function(which_idx, idx){
     Xs[[which_idx]][, idx] %*% gs[[which_idx]] +
     g_funcs[[which_idx]](ti) %*% betas[[which_idx]] +
     m_funcs[[which_idx]](ti) %*% zeta[rng_idx[[which_idx]]]
     }, x$is_obs, x$idx)

     x$y <- round(drop(
     mvtnorm::rmvnorm(1, mean = mea,
     sigma =  matrix(Sig, 3)[x$is_obs, x$is_obs])), 2)
     x
     })
     dput(obs)

# evaluate the lower bound
     f <- function(x){
     get_next <- function(n){
     out <- head(x, n)
     x <<- tail(x, -n)
     out
     }
     gs <- sapply(1:3, function(i) get_next(length(gs[[i]])))
     betas <- sapply(1:3, function(i) get_next(length(betas[[i]])))
     Sig <- matrix(get_next(length(Sig)), NROW(Sig))
     zeta <- get_next(length(zeta))
     Psi <- matrix(x, NROW(Psi))

     out <- 0
     for(x in obs){
     ti <- x$ti
     idx <- x$idx
     is_obs <- x$is_obs
     mea <- mapply(function(which_idx, idx){
     Xs[[which_idx]][, idx] %*% gs[[which_idx]] +
     g_funcs[[which_idx]](ti) %*% betas[[which_idx]] +
     m_funcs[[which_idx]](ti) %*% zeta[rng_idx[[which_idx]]]
     }, is_obs, idx)

     Sig <- (Sig + t(Sig)) / 2
     out <- out - mvtnorm::dmvnorm(x$y, mea, Sig[is_obs, is_obs], log = TRUE)

     M <- matrix(0, length(zeta), length(is_obs))
     for(i in seq_along(is_obs))
     M[rng_idx[[is_obs[i]]], i] <- m_funcs[[is_obs[i]]](ti)

     rng_indices <- do.call(c, rng_idx[is_obs])
     M <- M[rng_indices, ]
     out <- out + sum(diag(solve(
     Sig[is_obs, is_obs],
     crossprod(M, Psi[rng_indices, rng_indices]) %*% M))) / 2
     }

     out
     }

     dput(f(c(do.call(c, gs), do.call(c, betas), Sig, zeta, Psi)))
     gr <- numDeriv::grad(f, c(do.call(c, gs), do.call(c, betas), Sig, zeta, Psi))
     dput(head(gr, 18))
     dput(tail(gr, -18))
     */

    // the bases
    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;
    // raw poly of degree x without an intercept
    bases_fix.emplace_back(new joint_bases::orth_poly(1, false));
    bases_fix.emplace_back(new joint_bases::orth_poly(2, false));
    bases_fix.emplace_back(new joint_bases::orth_poly(1, false));
    // raw poly of degree x with an intercept
    bases_rng.emplace_back(new joint_bases::orth_poly(2, true));
    bases_rng.emplace_back(new joint_bases::orth_poly(1, true));
    bases_rng.emplace_back(new joint_bases::orth_poly(1, true));

    constexpr vajoint_uint n_obs[3]{2, 1, 2},
                         n_fixef[3]{2, 2, 1};
    constexpr double obs_1[n_obs[0]] {3.51, 8.34},
                     obs_2[n_obs[1]] {-3.41},
                     obs_3[n_obs[2]] {2.11, 8.69},

                obs_time_1[n_obs[0]] {.5, 1.2},
                obs_time_2[n_obs[1]] {.5},
                obs_time_3[n_obs[2]] {.5, 1.2};
    constexpr int ids_1[n_obs[0]] {1, 1},
                  ids_2[n_obs[1]] {1},
                  ids_3[n_obs[2]] {1, 1};

    // the design matrix
    double X1[n_obs[0] * n_fixef[0]] {1, 2, -1, -2},
           X2[n_obs[1] * n_fixef[1]] {-1, -2},
           X3[n_obs[2] * n_fixef[2]] {.5, .67};
    constexpr double Sig[]{3.22, 2.28, -2.76, 2.28, 10.26, -4.69, -2.76, -4.69, 7.34},
                     Psi[]{18.18, -1.3, 1.66, -1.75, -5.28, -0.24, -1, -1.3, 22.81, -3.08, 1.33, -0.69, 2.42, -2.52, 1.66, -3.08, 28.75, -5.05, 1.66, 5.43, -2.06, -1.75, 1.33, -5.05, 7.57, 0.46, -2.58, -1.31, -5.28, -0.69, 1.66, 0.46, 20.3, -5.26, 3.29, -0.24, 2.42, 5.43, -2.58, -5.26, 23.29, 3.9, -1, -2.52, -2.06, -1.31, 3.29, 3.9, 17.06},
                    zeta[]{3.42, -2.39, 7.86, -1.78, 6.86, 1.91, -1.06},
                      g1[]{-.25, .33},
                      g2[]{-1, 2},
                      g3[]{1.5},
            b1[n_fixef[0]]{.2},
            b2[n_fixef[1]]{.5, 3},
            b3[n_fixef[2]]{2.5},
             // has to change if the order of the parameters change
                  true_val{51.6018122187957},
             true_derivs[]{-0.193634481951502, -0.387268966740059, -0.486880407147759, -0.973760814811482, 0.206928716203386, 1.33555215993641, 0.243440204153009, 0.121720101804142, 0.082623776620504, -14.9974646929942, 0.0558312462478272, -6.53595373246103, 0.0558312462478272, -0.249914581141997, -0.295612614003981, -6.53595373246103, -0.295612614003981, -3.95101766793635,
                           // the random effect covariance matrix
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           1.49150599048705, 1.33555216123769, 1.37553507964905, 0.48688040733508, 0.243440203762803, 0.493600223878725, 0.0826237765032879, 0.462179286253018, 0.391479723534188, 0.388207958653328, -0.0165711893301057, -0.0082855945695312, 0.163201086019816, 0.141910708592852, 0.391479723534198, 0.388207958691539, 0.425065695562744, -0.00828559460639699, -0.00414279678863834, 0.141910708216482, 0.143327554016351, 0.388207958653328, 0.425065695562061, 0.489686907141934, -0.00414279732987994, -0.00207139866704958, 0.143327553996527, 0.158510416789727, -0.0165711893686029, -0.00828559460639699, -0.00414279733029684, 0.0700179520454573, 0.0350089755230614, 0.0385078627925807, 0.0192539311401717, -0.00828559458229066, -0.00414279738069372, -0.00207139866704958, 0.0350089755230614, 0.017504488005654, 0.0192539314160383, 0.00962696563434957, 0.163201082887329, 0.141910708357457, 0.143327554071765, 0.0385078626334194, 0.0192539314160383, 0.222212107192184, 0.181467914669852, 0.141910708527588, 0.143327554069825, 0.158510416791771, 0.0192539311915992, 0.00962696563434957, 0.181467914669852, 0.175168190474574 };

    // create the data objects
    std::vector<marker::setup_marker_dat_helper> input_dat;
    input_dat.emplace_back
      (X1, n_fixef[0], n_obs[0], ids_1, obs_time_1, obs_1);
    input_dat.emplace_back
      (X2, n_fixef[1], n_obs[1], ids_2, obs_time_2, obs_2);
    input_dat.emplace_back
      (X3, n_fixef[2], n_obs[2], ids_3, obs_time_3, obs_3);

    subset_params par_idx;
    for(size_t i = 0; i < 3; ++i)
      par_idx.add_marker({n_fixef[i], bases_fix[i]->get_n_basis(),
                         bases_rng[i]->get_n_basis()});

    marker::marker_dat comp_obj = marker::get_comp_dat
      (input_dat, par_idx, bases_fix, bases_rng);

    // test the basic members
    expect_true(comp_obj.n_obs == 2);
    expect_true(comp_obj.n_markers == 3);

    // the function returns the right value
    std::vector<double> par(par_idx.get_n_parms_w_va());

    std::fill(par.begin(), par.end(), 0);
    std::copy(
      begin(g1), end(g1), par.begin() + par_idx.get_fixef_idx_marker(0));
    std::copy(
      begin(g2), end(g2), par.begin() + par_idx.get_fixef_idx_marker(1));
    std::copy(
      begin(g3), end(g3), par.begin() + par_idx.get_fixef_idx_marker(2));

    std::copy(
      begin(b1), end(b1), par.begin() + par_idx.get_varying_idx_marker(0));
    std::copy(
      begin(b2), end(b2), par.begin() + par_idx.get_varying_idx_marker(1));
    std::copy(
      begin(b3), end(b3), par.begin() + par_idx.get_varying_idx_marker(2));

    std::copy(begin(Sig), end(Sig), par.begin() + par_idx.get_idx_error_term());
    std::copy(begin(zeta), end(zeta), par.begin() + par_idx.get_idx_va_mean());
    std::copy(begin(Psi), end(Psi), par.begin() + par_idx.get_idx_va_vcov());

    {
      double *wk_mem = wmem::get_double_mem(comp_obj.get_n_wmem());
      comp_obj.setup(par.data(), wk_mem);

      double res = comp_obj.eval(par.data(), wk_mem, 0) +
        comp_obj.eval(par.data(), wk_mem, 1);
      expect_true(res == Approx(true_val).epsilon(1e-7));
    }

    // test the gradient
    cfaad::Number::tape->rewind();
    std::vector<cfaad::Number> ad_par(par.size());
    cfaad::convertCollection(par.begin(), par.end(), ad_par.begin());

    cfaad::Number * wk_mem = wmem::get_Number_mem(comp_obj.get_n_wmem());
    cfaad::Number res = comp_obj.eval(ad_par.data(), wk_mem, 0) +
      comp_obj.eval(ad_par.data(), wk_mem, 1);

    res.propagateToStart();
    expect_true(res.value() == Approx(true_val).epsilon(1e-6));

    expect_true(ad_par.size() == static_cast<size_t>(
      std::distance(begin(true_derivs), end(true_derivs))));

    for(size_t i = 0; i < ad_par.size(); ++i)
      expect_true(ad_par[i].adjoint() == Approx(true_derivs[i]).epsilon(1e-6));

    // clean up
    wmem::clear_all();
  }
}
