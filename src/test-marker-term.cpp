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
             true_derivs[]{-3.84, 1.92, 6.08, -5.12,
                           // the survival parameters
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the error term covariance matrix
                           -1121.3248,
                           // the random effect covariance matrix
                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the frailty effects covariance matrix
                           0, 0, 0, 0,
                           // the VA mean
                           -3.84, -5.1199999, -8.96, 0, 0,
                           // the VA covariance matrix
                           4, 8, 20, 0, 0,
                           8, 20, 56, 0, 0,
                           20, 56, 164, 0, 0,
                           0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0};

    // create the input data
    std::vector<marker::setup_marker_dat_helper> input_dat;
    input_dat.emplace_back(X, 3, n_obs, ids, obs_time, obs,
                           nullptr, 0, nullptr, 0);

    subset_params par_idx;
    par_idx.add_marker({n_fixef, 1, 3});
    par_idx.add_surv({2, 1, {1}, true});
    par_idx.add_surv({3, 2, {1}, true});

    auto dat_n_idx = marker::get_comp_dat
      (input_dat, par_idx, bases_fix, bases_rng);
    marker::marker_dat &comp_obj = dat_n_idx.dat;

    // test the basic members
    expect_true(comp_obj.n_obs() == n_obs);
    expect_true(comp_obj.n_markers() == 1);
    expect_true(dat_n_idx.id.size() == n_obs);
    for(int x : dat_n_idx.id)
      expect_true(x == 1);

    // the function returns the right value
    std::vector<double> par(par_idx.n_params_w_va());

    std::fill(par.begin(), par.end(), 0);
    std::copy(begin(g), end(g), par.begin() + par_idx.fixef_marker(0));
    std::copy(begin(beta), end(beta),
              par.begin() + par_idx.fixef_vary_marker(0));
    std::copy(begin(Sig), end(Sig), par.begin() + par_idx.vcov_marker());
    std::copy(begin(zeta), end(zeta), par.begin() + par_idx.va_mean());

    const vajoint_uint n_shared_p_surv
      {par_idx.n_shared() + par_idx.n_shared_surv()};
    for(vajoint_uint j = 0; j < par_idx.n_shared(); ++j)
        std::copy(
          Psi + j * par_idx.n_shared(),
          Psi + (j + 1) * par_idx.n_shared(),
          par.begin() + par_idx.va_vcov() + j * n_shared_p_surv);

    {
      double *wk_mem = wmem::get_double_mem(comp_obj.n_wmem());
      comp_obj.setup(par.data(), wk_mem);

      double res = comp_obj(par.data(), wk_mem, 0)
        + comp_obj(par.data(), wk_mem, 1);

      expect_true(pass_rel_err(res, true_val, 1e-7));
    }

    // test the gradient
    cfaad::Number::tape->rewind();
    std::vector<cfaad::Number> ad_par(par.size());
    cfaad::convertCollection(par.begin(), par.end(), ad_par.begin());

    cfaad::Number * wk_mem = wmem::get_Number_mem(comp_obj.n_wmem());
    cfaad::Number res = comp_obj(ad_par.data(), wk_mem, 0) +
      comp_obj(ad_par.data(), wk_mem, 1);

    res.propagateToStart();
    expect_true(pass_rel_err(res.value(), true_val, 1e-6));

    expect_true(ad_par.size() == static_cast<size_t>(
      std::distance(begin(true_derivs), end(true_derivs))));

    for(size_t i = 0; i < ad_par.size(); ++i)
      expect_true(pass_rel_err(ad_par[i].adjoint(), true_derivs[i], 1e-6));

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
     matrix(c(-(1:2), .23, .1), 2),
     matrix(c(.5, .67, -.5, 4.2, .7), 1))
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
     list(is_obs = c(1L, 3L), ti = 1.2, idx = c(2L, 2L)),
     list(is_obs = 3L, ti = 2, idx = 3L),
     list(is_obs = 2:3, ti = 2.5, idx = c(2L, 4L)),
     list(is_obs = 3L, ti = 3, idx = 5L))

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
     mvtnorm::rmvnorm(
     1, mean = mea,
     sigma =  matrix(Sig, 3)[x$is_obs, x$is_obs, drop = FALSE])), 2)
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
     out <- out - mvtnorm::dmvnorm(x$y, mea, Sig[is_obs, is_obs, drop = FALSE],
     log = TRUE)

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

    constexpr vajoint_uint n_obs[3]{2, 2, 5},
                         n_fixef[3]{2, 2, 1};
    constexpr double obs_1[n_obs[0]] {3.51, 8.34},
                     obs_2[n_obs[1]] {-3.41, 34.42},
                     obs_3[n_obs[2]] {2.11, 8.69, 6.29, 12.68, 6.26},

                obs_time_1[n_obs[0]] {.5, 1.2},
                obs_time_2[n_obs[1]] {.5, 2.5},
                obs_time_3[n_obs[2]] {.5, 1.2, 2, 2.5, 3};
    constexpr int ids_1[n_obs[0]] {1, 1},
                  ids_2[n_obs[1]] {1, 1},
                  ids_3[n_obs[2]] {1, 1, 1, 1, 1};

    // the design matrix
    double X1[n_obs[0] * n_fixef[0]] {1, 2, -1, -2},
           X2[n_obs[1] * n_fixef[1]] {-1, -2, .23, .1},
           X3[n_obs[2] * n_fixef[2]] {.5, .67, -.5, 4.2, .7};
    constexpr double Sig[]{3.22, 2.28, -2.76, 2.28, 10.26, -4.69, -2.76, -4.69, 7.34},
                     Psi[]{18.18, -1.3, 1.66, -1.75, -5.28, -0.24, -1, -1.3, 22.81, -3.08, 1.33, -0.69, 2.42, -2.52, 1.66, -3.08, 28.75, -5.05, 1.66, 5.43, -2.06, -1.75, 1.33, -5.05, 7.57, 0.46, -2.58, -1.31, -5.28, -0.69, 1.66, 0.46, 20.3, -5.26, 3.29, -0.24, 2.42, 5.43, -2.58, -5.26, 23.29, 3.9, -1, -2.52, -2.06, -1.31, 3.29, 3.9, 17.06},
                    zeta[]{3.42, -2.39, 7.86, -1.78, 6.86, 1.91, -1.06},
                      g1[n_fixef[0]]{-.25, .33},
                      g2[n_fixef[1]]{-1, 2},
                      g3[n_fixef[2]]{1.5},
                      b1[]{.2},
                      b2[]{.5, 3, .11},
                      b3[]{2.5, .03, -.05, .15},
             // has to change if the order of the parameters change
                  true_val{104.600753185115},
             true_derivs[]{-0.193634482490612, -0.387268970881493, -0.475350703945104, -0.968747900047675, 0.0941826527723997, 1.3355521630108, 0.368763064148612, 0.435027255805212, -0.32980616617052,
                           // the survival parameters
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the error term covariance matrix
                           -14.997464692592, 0.0558312458857624, -6.53595373214058, 0.0558312458857624, -2.07540609242216, -2.36906621721635, -6.53595373214058, -2.36906621721635, -9.95390817846328,
                           // the random effect covariance matrix
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the frailty term covariance matrix
                           0, 0, 0, 0,
                           // the VA mean
                           1.49150599050613, 1.33555216121127, 1.37553507975353, 0.537009551710603, 0.368763065310235, 0.239527429124941, -0.329806166426263, 0, 0,
                           // the VA covariance matrix
                           0.462179286279196, 0.391479722800414, 0.388207958157325, -0.0165711898018102, -0.00828559438766372, 0.163201084019404, 0.141910708180098, 0, 0,
                           0.391479722800414, 0.38820795867049, 0.425065695716578, -0.00828559398890858, -0.00414279778357616, 0.14191070799113, 0.143327554369823,0, 0,
                           0.388207958157325, 0.425065695716578, 0.489686907168083, -0.00414279749333691, -0.00207139883325205, 0.143327553908884, 0.158510416423749, 0, 0,
                           -0.0165711899557987, -0.00828559398890858, -0.00414279749333691, 0.138857600708758, 0.207108096726541, 0.0824939595780647, 0.129219173816305, 0, 0,
                           -0.00828559438766372, -0.00414279777137145, -0.00207139883325205, 0.207108096726541, 0.447752291924522, 0.129219173997635, 0.284540072044498, 0, 0,
                           0.163201084054493, 0.141910708781202, 0.143327553958512, 0.0824939595241982, 0.12921917399764, 0.454677338135144, 0.762630992085864, 0, 0,
                           0.141910709822643, 0.143327554369823, 0.158510416423749, 0.129219174546085, 0.284540072044498, 0.762630992085871, 1.66213582927789, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    // create the data objects
    std::vector<marker::setup_marker_dat_helper> input_dat;
    input_dat.emplace_back
      (X1, n_fixef[0], n_obs[0], ids_1, obs_time_1, obs_1,
       nullptr, 0, nullptr, 0);
    input_dat.emplace_back
      (X2, n_fixef[1], n_obs[1], ids_2, obs_time_2, obs_2,
       nullptr, 0, nullptr, 0);
    input_dat.emplace_back
      (X3, n_fixef[2], n_obs[2], ids_3, obs_time_3, obs_3,
       nullptr, 0, nullptr, 0);

    subset_params par_idx;
    for(size_t i = 0; i < 3; ++i)
      par_idx.add_marker({n_fixef[i], bases_fix[i]->n_basis(),
                         bases_rng[i]->n_basis()});
    par_idx.add_surv({2, 1, {1, 1, 1}, true});
    par_idx.add_surv({3, 2, {1, 1, 1}, true});

    auto dat_n_idx = marker::get_comp_dat
      (input_dat, par_idx, bases_fix, bases_rng);
    marker::marker_dat &comp_obj = dat_n_idx.dat;

    // test the basic members
    expect_true(comp_obj.n_obs() == 5);
    expect_true(comp_obj.n_markers() == 3);
    expect_true(dat_n_idx.id.size() == 5);
    for(int x : dat_n_idx.id)
      expect_true(x == 1);

    // the function returns the right value
    std::vector<double> par(par_idx.n_params_w_va());

    std::fill(par.begin(), par.end(), 0);
    std::copy(
      begin(g1), end(g1), par.begin() + par_idx.fixef_marker(0));
    std::copy(
      begin(g2), end(g2), par.begin() + par_idx.fixef_marker(1));
    std::copy(
      begin(g3), end(g3), par.begin() + par_idx.fixef_marker(2));

    std::copy(
      begin(b1), end(b1), par.begin() + par_idx.fixef_vary_marker(0));
    std::copy(
      begin(b2), end(b2), par.begin() + par_idx.fixef_vary_marker(1));
    std::copy(
      begin(b3), end(b3), par.begin() + par_idx.fixef_vary_marker(2));

    std::copy(begin(Sig), end(Sig), par.begin() + par_idx.vcov_marker());
    std::copy(begin(zeta), end(zeta), par.begin() + par_idx.va_mean());
    const vajoint_uint n_shared_p_surv
    {par_idx.n_shared() + par_idx.n_shared_surv()};
    for(vajoint_uint j = 0; j < par_idx.n_shared(); ++j)
      std::copy(
        Psi + j * par_idx.n_shared(),
        Psi + (j + 1) * par_idx.n_shared(),
        par.begin() + par_idx.va_vcov() + j * n_shared_p_surv);

    {
      double *wk_mem = wmem::get_double_mem(comp_obj.n_wmem());
      comp_obj.setup(par.data(), wk_mem);

      double res{};
      for(vajoint_uint i = 0; i < comp_obj.n_obs(); ++i)
        res += comp_obj(par.data(), wk_mem, i);
      expect_true(pass_rel_err(res, true_val, 1e-7));
    }

    // test the gradient
    cfaad::Number::tape->rewind();
    std::vector<cfaad::Number> ad_par(par.size());
    cfaad::convertCollection(par.begin(), par.end(), ad_par.begin());

    cfaad::Number * wk_mem = wmem::get_Number_mem(comp_obj.n_wmem());

    cfaad::Number res{0};
    for(vajoint_uint i = 0; i < comp_obj.n_obs(); ++i)
      res += comp_obj(ad_par.data(), wk_mem, i);

    res.propagateToStart();
    expect_true(pass_rel_err(res.value(), true_val, 1e-6));

    expect_true(ad_par.size() == static_cast<size_t>(
      std::distance(begin(true_derivs), end(true_derivs))));

    for(size_t i = 0; i < ad_par.size(); ++i)
      expect_true(pass_rel_err(ad_par[i].adjoint(), true_derivs[i], 1e-6));

    // clean up
    wmem::clear_all();
  }

  test_that("marker_term gives the correct result with three markers and time-varying effects"){
    /* R code to reproduce the result
     raw_poly <- function(x, degree, intercept){
     if(intercept)
     outer(x, 0:degree, `^`)
     else
     outer(x, 1:degree, `^`)
     }

# parameters
     set.seed(1)
     Xs <- list(matrix(c(1:2, -1, -2), 2),
     matrix(c(-(1:2), .23, .1), 2),
     matrix(c(.5, .67, -.5, 4.2, .7), 1))
     gs <- list(c(-.25, .33), c(-1, 2), c(1.5))
     betas <- list(c(.2), c(.5, 3, .11), c(2.5, .03, -.05, .15))
     dput(Sig <- round(drop(rWishart(1, 6, diag(3))), 2))
     stopifnot(all(eigen(Sig)$values > 0))
     dput(Psi <- round(drop(rWishart(1, 22, diag(1/22, 11))), 2))
     dput(zeta <- drop(round(mvtnorm::rmvnorm(1, sigma = Psi), 2)))
     rng_idx <- list(1:3, 4:7, 8:11)

     g_funcs <- list(
     function(x, data)
     raw_poly(x, 1, FALSE),
     function(x, data)
     cbind(raw_poly(x, 2, FALSE),
     raw_poly(x, 1, FALSE) * data$w1),
     function(x, data)
     cbind(raw_poly(x, 1, FALSE),
     raw_poly(x, 2, FALSE) * data$w1,
     raw_poly(x, 1, FALSE) * data$w2))
     m_funcs <- list(
     function(x, data) raw_poly(x, 2, TRUE),
     function(x, data)
     cbind(raw_poly(x, 1, TRUE),
     raw_poly(x, 1, TRUE) * data$w1),
     function(x, data)
     cbind(raw_poly(x, 1, TRUE),
     raw_poly(x, 1, FALSE) * data$w1,
     raw_poly(x, 1, FALSE) * data$w1))

     obs <- list(
     list(is_obs = 1:3, ti = .5, idx = c(1L, 1L, 1L)),
     list(is_obs = c(1L, 3L), ti = 1.2, idx = c(2L, 2L)),
     list(is_obs = 3L, ti = 2, idx = 3L),
     list(is_obs = 2:3, ti = 2.5, idx = c(2L, 4L)),
     list(is_obs = 3L, ti = 3, idx = 5L))

     dats <- list(
     matrix(0., 2L, 0L) |> as.data.frame(),
     data.frame(w1 = c(0.071, -0.022)),
     data.frame(w1 = c(-0.024, 0.079, 0.029, 0.048, 0.021),
     w2 = c(0.081, -0.041, -0.062, 0.077, 0.001)))
     sapply(dats, \(x) as.matrix(x) |> t() |> dput()) |>
     invisible()

# generate a plausible
     obs <- lapply(obs, function(x){
# compute the mean
     ti <- x$ti
     mea <- mapply(function(which_idx, idx){
     Xs[[which_idx]][, idx] %*% gs[[which_idx]] +
     g_funcs[[which_idx]](ti, dats[[which_idx]][idx, , drop = FALSE]) %*%
     betas[[which_idx]] +
     m_funcs[[which_idx]](ti, dats[[which_idx]][idx, , drop = FALSE]) %*%
     zeta[rng_idx[[which_idx]]]
     }, x$is_obs, x$idx)

     x$y <- round(drop(
     mvtnorm::rmvnorm(
     1, mean = mea,
     sigma =  matrix(Sig, 3)[x$is_obs, x$is_obs, drop = FALSE])), 2)
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
     g_funcs[[which_idx]](ti, dats[[which_idx]][idx, , drop = FALSE]) %*%
     betas[[which_idx]] +
     m_funcs[[which_idx]](ti, dats[[which_idx]][idx, , drop = FALSE]) %*%
     zeta[rng_idx[[which_idx]]]
     }, is_obs, idx)

     Sig <- (Sig + t(Sig)) / 2
     out <- out - mvtnorm::dmvnorm(x$y, mea, Sig[is_obs, is_obs, drop = FALSE],
     log = TRUE)

     M <- matrix(0, length(zeta), length(is_obs))
     for(i in seq_along(is_obs))
     M[rng_idx[[is_obs[i]]], i] <-
     m_funcs[[is_obs[i]]](ti, dats[[is_obs[i]]][x$idx[i], , drop = FALSE])

     rng_indices <- do.call(c, rng_idx[is_obs])
     M <- M[rng_indices, ]
     out <- out + sum(diag(solve(
     Sig[is_obs, is_obs],
     crossprod(M, Psi[rng_indices, rng_indices]) %*% M))) / 2
     }

     out
     }

     dput(f(c(do.call(c, gs), do.call(c, betas), Sig, zeta, Psi)))
     gr <- numDeriv::grad(f, c(do.call(c, gs), do.call(c, betas), Sig, zeta, Psi),
     method.args = list(d = 1e-2, r = 6))
     dput(head(gr, -11 * 12))
     va_mean_n_vcov <- tail(gr, 11 * 12)
     dput(c(head(va_mean_n_vcov, 11), 0, 0))

     dum <- matrix(0., 13, 13)
     dum[1:11, 1:11] <- tail(va_mean_n_vcov, -11)
     dput(dum)

     */

    // the bases
    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;
    // raw poly of degree x without an intercept
    bases_fix.emplace_back(new joint_bases::orth_poly(1, false));
    {
      joint_bases::bases_vector expansions;
      expansions.emplace_back(new joint_bases::orth_poly(2, false));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly> (1, false));
      bases_fix.emplace_back(new joint_bases::stacked_basis(expansions));
    }
    {
      joint_bases::bases_vector expansions;
      expansions.emplace_back(new joint_bases::orth_poly(1, false));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly> (2, false));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly> (1, false));
      bases_fix.emplace_back(new joint_bases::stacked_basis(expansions));
    }
    // raw poly of degree x with an intercept
    bases_rng.emplace_back(new joint_bases::orth_poly(2, true));
    {
      joint_bases::bases_vector expansions;
      expansions.emplace_back(new joint_bases::orth_poly(1, true));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly> (1, true));
      bases_rng.emplace_back(new joint_bases::stacked_basis(expansions));
    }
    {
      joint_bases::bases_vector expansions;
      expansions.emplace_back(new joint_bases::orth_poly(1, true));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly> (1, false));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly> (1, false));
      bases_rng.emplace_back(new joint_bases::stacked_basis(expansions));
    }

    constexpr vajoint_uint n_obs[3]{2, 2, 5},
                         n_fixef[3]{2, 2, 1};
    constexpr double obs_1[n_obs[0]] {2.13, .29},
                     obs_2[n_obs[1]] {-3.11, 21.42},
                     obs_3[n_obs[2]] {2.99, 2.10, 6.31, 11.70, 10.56},

                obs_time_1[n_obs[0]] {.5, 1.2},
                obs_time_2[n_obs[1]] {.5, 2.5},
                obs_time_3[n_obs[2]] {.5, 1.2, 2, 2.5, 3};
    constexpr int ids_1[n_obs[0]] {1, 1},
                  ids_2[n_obs[1]] {1, 1},
                  ids_3[n_obs[2]] {1, 1, 1, 1, 1};

    // the design matrix
    double X1[n_obs[0] * n_fixef[0]] {1, 2, -1, -2},
           X2[n_obs[1] * n_fixef[1]] {-1, -2, .23, .1},
           X3[n_obs[2] * n_fixef[2]] {.5, .67, -.5, 4.2, .7};
    constexpr double Sig[]{3.22, 2.28, -2.76, 2.28, 10.26, -4.69, -2.76, -4.69, 7.34},
                     Psi[]{0.87, -0.06, 0.08, -0.08, -0.25, -0.01, -0.05, -0.09, -0.14, 0.05, 0.04, -0.06, 1.09, -0.14, 0.06, -0.03, 0.11, -0.12, 0.28, 0.09, 0.13, -0.4, 0.08, -0.14, 1.36, -0.23, 0.08, 0.25, -0.1, -0.11, 0.17, -0.06, 0.42, -0.08, 0.06, -0.23, 0.37, 0.02, -0.12, -0.06, 0.25, -0.03, -0.26, -0.06, -0.25, -0.03, 0.08, 0.02, 0.97, -0.25, 0.15, 0.13, 0.23, -0.29, 0.47, -0.01, 0.11, 0.25, -0.12, -0.25, 1.11, 0.18, -0.17, 0.08, 0.22, 0, -0.05, -0.12, -0.1, -0.06, 0.15, 0.18, 0.82, -0.23, -0.07, 0.03, -0.02, -0.09, 0.28, -0.11, 0.25, 0.13, -0.17, -0.23, 1.17, 0.12, -0.37, 0.08, -0.14, 0.09, 0.17, -0.03, 0.23, 0.08, -0.07, 0.12, 0.47, -0.07, 0.02, 0.05, 0.13, -0.06, -0.26, -0.29, 0.22, 0.03, -0.37, -0.07, 0.69, -0.37, 0.04, -0.4, 0.42, -0.06, 0.47, 0, -0.02, 0.08, 0.02, -0.37, 0.97},
                    zeta[]{0.46, -0.42, -0.2, 0.11, -0.68, -0.64, -0.28, 0.99, -1.1, 0.19, 0.2},
                      g1[n_fixef[0]]{-.25, .33},
                      g2[n_fixef[1]]{-1, 2},
                      g3[n_fixef[2]]{1.5},
                      b1[]{.2},
                      b2[]{.5, 3, .11},
                      b3[]{2.5, .03, -.05, .15},
             // has to change if the order of the parameters change
                  true_val{23.406006390383},
             true_derivs[]{-0.686930343431596, -1.37386068686004, -0.233465081814553, -0.290631135879994, -2.00569598187874, -0.491254696988613, -1.16388400801072, -3.03053961186787, 0.0312240241218763, -3.64151585302462, -0.0930940102372379, -0.268069660277636, -0.0497451358582203,
                           // the survival parameters
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the error term covariance matrix
                           -0.557565408275804, 0.0366041369901646, -0.255577218772604, 0.0366041369901646, -0.107608637686856, -0.135589573426681, -0.255577218772604, -0.135589573426681, -0.38796347753459,
                           // the random effect covariance matrix
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the frailty term covariance matrix
                           0, 0, 0, 0,
                           // the VA mean
                           -0.860800373148064, -0.49125469696673, -0.318652760967096, -0.368889929662617, -1.16388400800676, 0.0193527304916874, 0.0312240242039636, -1.56282354226133, -3.6415158530316, -0.0930940106065062, -0.0930940106291562, 0, 0,
                           // the VA covariance matrix
                           0.462179286289639, 0.391479723756652, 0.388207958617115, -0.0165711893079348, -0.00828559464539331, -0.00117655405009568, -0.000588277202943746, 0.16320108340842, 0.141910708243008, 0.00724319522849319, 0.00724319540775074, 0, 0, 0.391479723779554, 0.38820795870364, 0.425065695505314, -0.00828559457641688, -0.00414279720759609, -0.00058827718094002, -0.000294138581301791, 0.141910708301795, 0.143327554001576, 0.00933900143781742, 0.00933900142285715, 0, 0, 0.388207958617115, 0.425065695505314, 0.489686907132506, -0.00414279729576096, -0.00207139860990391, -0.000294138612663024, -0.000147069259472782, 0.143327553999435, 0.158510416860996, 0.0115303853198041, 0.0115303852743959, 0, 0, -0.0165711893079348, -0.00828559457641688, -0.00414279729576096, 0.138857600662365, 0.207108097640702, 0.00345680228851865, -0.00130054341148032, 0.0824939597805496, 0.129219174070462, 0.00481623727755596, 0.00481623736220063, 0, 0, -0.00828559464539331, -0.00414279720759609, -0.00207139860990391, 0.207108097640702, 0.447752291989805, -0.00130054339083194, -0.00822263302840264, 0.129219173939494, 0.28454007203036, 0.0129647819367505, 0.0129647819303895, 0, 0, -0.00117655405009568, -0.00058827718094002, -0.000294138612663024, 0.00345680228851865, -0.00130054339083194, 0.000386278883236921, 0.00025977622423747, 0.00176636409318064, -0.00105220621549495, -0.000148932004905833, -0.000148932168603116, 0, 0, -0.000588277202943746, -0.000294138581301791, -0.000147069259472782, -0.00130054341148032, -0.00822263302840264, 0.00025977622423747, 0.000296480063587383, -0.00105220621508964, -0.0053645738003836, -0.00030671245555132, -0.000306712593118933, 0, 0, 0.16320108340842, 0.141910708301795, 0.143327553999435, 0.0824939597805496, 0.129219173939494, 0.00176636409318064, -0.00105220621508964, 0.454677338112836, 0.762630991875527, 0.0278582250400917, 0.0278582251111442, 0, 0, 0.141910708243008, 0.143327554002083, 0.158510416860996, 0.129219174070462, 0.28454007203036, -0.00105220621549495, -0.0053645738003836, 0.762630991875527, 1.66213582930704, 0.0603488381019782, 0.0603488378982845, 0, 0, 0.00724319522849319, 0.00933900143781742, 0.0115303853198041, 0.00481623727755596, 0.0129647819367505, -0.000148932004904213, -0.00030671245555132, 0.0278582250623971, 0.0603488381013264, 0.00280604352141154, 0.00280604352067029, 0, 0, 0.00724319540775074, 0.00933900142285715, 0.0115303852743959, 0.00481623736220063, 0.0129647819303895, -0.000148932168603116, -0.000306712593118933, 0.0278582251111442, 0.0603488378982845, 0.00280604352067029, 0.00280604353653753, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    // create the data objects
    std::vector<marker::setup_marker_dat_helper> input_dat;
    input_dat.emplace_back
      (X1, n_fixef[0], n_obs[0], ids_1, obs_time_1, obs_1,
       nullptr, 0, nullptr, 0);

    double fixef_varying2[]{0.071, -0.022},
           rng_varying2[]{0.071, -0.022};
    input_dat.emplace_back
      (X2, n_fixef[1], n_obs[1], ids_2, obs_time_2, obs_2,
       fixef_varying2, 1, fixef_varying2, 1);

    double fixef_varying3[]{-0.024, 0.081, 0.079, -0.041, 0.029, -0.062, 0.048, 0.077, 0.021, 0.001},
             rng_varying3[]{-0.024, -0.024, 0.079, 0.079, 0.029, 0.029, 0.048, 0.048, 0.021,0.021};
    input_dat.emplace_back
      (X3, n_fixef[2], n_obs[2], ids_3, obs_time_3, obs_3,
       fixef_varying3, 2, rng_varying3, 2);

    subset_params par_idx;
    for(size_t i = 0; i < 3; ++i)
      par_idx.add_marker({n_fixef[i], bases_fix[i]->n_basis(),
                         bases_rng[i]->n_basis()});
    par_idx.add_surv({2, 1, {1, 1, 1}, true});
    par_idx.add_surv({3, 2, {1, 1, 1}, true});

    auto dat_n_idx = marker::get_comp_dat
      (input_dat, par_idx, bases_fix, bases_rng);
    marker::marker_dat &comp_obj = dat_n_idx.dat;

    // test the basic members
    expect_true(comp_obj.n_obs() == 5);
    expect_true(comp_obj.n_markers() == 3);
    expect_true(dat_n_idx.id.size() == 5);
    for(int x : dat_n_idx.id)
      expect_true(x == 1);

    // the function returns the right value
    std::vector<double> par(par_idx.n_params_w_va());

    std::fill(par.begin(), par.end(), 0);
    std::copy(
      begin(g1), end(g1), par.begin() + par_idx.fixef_marker(0));
    std::copy(
      begin(g2), end(g2), par.begin() + par_idx.fixef_marker(1));
    std::copy(
      begin(g3), end(g3), par.begin() + par_idx.fixef_marker(2));

    std::copy(
      begin(b1), end(b1), par.begin() + par_idx.fixef_vary_marker(0));
    std::copy(
      begin(b2), end(b2), par.begin() + par_idx.fixef_vary_marker(1));
    std::copy(
      begin(b3), end(b3), par.begin() + par_idx.fixef_vary_marker(2));

    std::copy(begin(Sig), end(Sig), par.begin() + par_idx.vcov_marker());
    std::copy(begin(zeta), end(zeta), par.begin() + par_idx.va_mean());
    const vajoint_uint n_shared_p_surv
    {par_idx.n_shared() + par_idx.n_shared_surv()};
    for(vajoint_uint j = 0; j < par_idx.n_shared(); ++j)
      std::copy(
        Psi + j * par_idx.n_shared(),
        Psi + (j + 1) * par_idx.n_shared(),
        par.begin() + par_idx.va_vcov() + j * n_shared_p_surv);

    {
      double *wk_mem = wmem::get_double_mem(comp_obj.n_wmem());
      comp_obj.setup(par.data(), wk_mem);

      double res{};
      for(vajoint_uint i = 0; i < comp_obj.n_obs(); ++i)
        res += comp_obj(par.data(), wk_mem, i);
      expect_true(pass_rel_err(res, true_val, 1e-7));
    }

    // test the gradient
    cfaad::Number::tape->rewind();
    std::vector<cfaad::Number> ad_par(par.size());
    cfaad::convertCollection(par.begin(), par.end(), ad_par.begin());


    cfaad::Number * wk_mem = wmem::get_Number_mem(comp_obj.n_wmem());

    cfaad::Number res{0};
    for(vajoint_uint i = 0; i < comp_obj.n_obs(); ++i)
      res += comp_obj(ad_par.data(), wk_mem, i);

    res.propagateToStart();
    expect_true(pass_rel_err(res.value(), true_val, 1e-6));

    size_t const expected_length =
      std::distance(begin(true_derivs), end(true_derivs));
    expect_true(ad_par.size() == expected_length);

    for(size_t i = 0; i < ad_par.size(); ++i){
      expect_true(pass_rel_err(ad_par[i].adjoint(), true_derivs[i], 1e-5));
    }

    // clean up
    wmem::clear_all();
  }

  test_that("marker_term gives the correct result with three markers and redudant survival terms"){
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

    constexpr vajoint_uint n_obs[3]{2, 2, 5},
                         n_fixef[3]{2, 2, 1};
    constexpr double obs_1[n_obs[0]] {3.51, 8.34},
                     obs_2[n_obs[1]] {-3.41, 34.42},
                     obs_3[n_obs[2]] {2.11, 8.69, 6.29, 12.68, 6.26},
                obs_time_1[n_obs[0]] {.5, 1.2},
                obs_time_2[n_obs[1]] {.5, 2.5},
                obs_time_3[n_obs[2]] {.5, 1.2, 2, 2.5, 3};
    constexpr int ids_1[n_obs[0]] {1, 1},
                  ids_2[n_obs[1]] {1, 1},
                  ids_3[n_obs[2]] {1, 1, 1, 1, 1};

    // the design matrix
    double X1[n_obs[0] * n_fixef[0]] {1, 2, -1, -2},
           X2[n_obs[1] * n_fixef[1]] {-1, -2, .23, .1},
           X3[n_obs[2] * n_fixef[2]] {.5, .67, -.5, 4.2, .7};
    constexpr double Sig[]{3.22, 2.28, -2.76, 2.28, 10.26, -4.69, -2.76, -4.69, 7.34},
                     Psi[]{18.18, -1.3, 1.66, -1.75, -5.28, -0.24, -1, -1.3, 22.81, -3.08, 1.33, -0.69, 2.42, -2.52, 1.66, -3.08, 28.75, -5.05, 1.66, 5.43, -2.06, -1.75, 1.33, -5.05, 7.57, 0.46, -2.58, -1.31, -5.28, -0.69, 1.66, 0.46, 20.3, -5.26, 3.29, -0.24, 2.42, 5.43, -2.58, -5.26, 23.29, 3.9, -1, -2.52, -2.06, -1.31, 3.29, 3.9, 17.06},
                    zeta[]{3.42, -2.39, 7.86, -1.78, 6.86, 1.91, -1.06},
                      g1[n_fixef[0]]{-.25, .33},
                      g2[n_fixef[1]]{-1, 2},
                      g3[n_fixef[2]]{1.5},
                      b1[]{.2},
                      b2[]{.5, 3},
                      b3[]{2.5},
             // has to change if the order of the parameters change
                  true_val{104.600753185115},
             true_derivs[]{-0.193634482490612, -0.387268970881493, -0.475350703945104, -0.968747900047675, 0.0941826527723997, 1.3355521630108, 0.368763064148612, 0.435027255805212, -0.32980616617052,
                           // the survival parameters
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the error term covariance matrix
                           -14.997464692592, 0.0558312458857624, -6.53595373214058, 0.0558312458857624, -2.07540609242216, -2.36906621721635, -6.53595373214058, -2.36906621721635, -9.95390817846328,
                           // the random effect covariance matrix
                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                           // the VA mean
                           1.49150599050613, 1.33555216121127, 1.37553507975353, 0.537009551710603, 0.368763065310235, 0.239527429124941, -0.329806166426263,
                           // the VA covariance matrix
                           0.462179286279196, 0.391479722800414, 0.388207958157325, -0.0165711898018102, -0.00828559438766372, 0.163201084019404, 0.141910708180098,
                           0.391479722800414, 0.38820795867049, 0.425065695716578, -0.00828559398890858, -0.00414279778357616, 0.14191070799113, 0.143327554369823,
                           0.388207958157325, 0.425065695716578, 0.489686907168083, -0.00414279749333691, -0.00207139883325205, 0.143327553908884, 0.158510416423749,
                           -0.0165711899557987, -0.00828559398890858, -0.00414279749333691, 0.138857600708758, 0.207108096726541, 0.0824939595780647, 0.129219173816305,
                           -0.00828559438766372, -0.00414279777137145, -0.00207139883325205, 0.207108096726541, 0.447752291924522, 0.129219173997635, 0.284540072044498,
                           0.163201084054493, 0.141910708781202, 0.143327553958512, 0.0824939595241982, 0.12921917399764, 0.454677338135144, 0.762630992085864,
                           0.141910709822643, 0.143327554369823, 0.158510416423749, 0.129219174546085, 0.284540072044498, 0.762630992085871, 1.66213582927789 };

    // create the data objects
    std::vector<marker::setup_marker_dat_helper> input_dat;
    input_dat.emplace_back
      (X1, n_fixef[0], n_obs[0], ids_1, obs_time_1, obs_1,
       nullptr, 0, nullptr, 0);
    input_dat.emplace_back
      (X2, n_fixef[1], n_obs[1], ids_2, obs_time_2, obs_2,
       nullptr, 0, nullptr, 0);
    input_dat.emplace_back
      (X3, n_fixef[2], n_obs[2], ids_3, obs_time_3, obs_3,
       nullptr, 0, nullptr, 0);

    subset_params par_idx;
    for(size_t i = 0; i < 3; ++i)
      par_idx.add_marker({n_fixef[i], bases_fix[i]->n_basis(),
                         bases_rng[i]->n_basis()});
    par_idx.add_surv({2, 1, {1, 1, 1}, false});
    par_idx.add_surv({3, 2, {1, 1, 1}, false});

    auto dat_n_idx = marker::get_comp_dat
      (input_dat, par_idx, bases_fix, bases_rng);
    marker::marker_dat &comp_obj = dat_n_idx.dat;

    // test the basic members
    expect_true(comp_obj.n_obs() == 5);
    expect_true(comp_obj.n_markers() == 3);
    expect_true(dat_n_idx.id.size() == 5);
    for(int x : dat_n_idx.id)
      expect_true(x == 1);

    // the function returns the right value
    std::vector<double> par(par_idx.n_params_w_va());

    std::fill(par.begin(), par.end(), 0);
    std::copy(
      begin(g1), end(g1), par.begin() + par_idx.fixef_marker(0));
    std::copy(
      begin(g2), end(g2), par.begin() + par_idx.fixef_marker(1));
    std::copy(
      begin(g3), end(g3), par.begin() + par_idx.fixef_marker(2));

    std::copy(
      begin(b1), end(b1), par.begin() + par_idx.fixef_vary_marker(0));
    std::copy(
      begin(b2), end(b2), par.begin() + par_idx.fixef_vary_marker(1));
    std::copy(
      begin(b3), end(b3), par.begin() + par_idx.fixef_vary_marker(2));

    std::copy(begin(Sig), end(Sig), par.begin() + par_idx.vcov_marker());
    std::copy(begin(zeta), end(zeta), par.begin() + par_idx.va_mean());
    const vajoint_uint n_shared_p_surv
    {par_idx.n_shared() + par_idx.n_shared_surv()};
    for(vajoint_uint j = 0; j < par_idx.n_shared(); ++j)
      std::copy(
        Psi + j * par_idx.n_shared(),
        Psi + (j + 1) * par_idx.n_shared(),
        par.begin() + par_idx.va_vcov() + j * n_shared_p_surv);

    {
      double *wk_mem = wmem::get_double_mem(comp_obj.n_wmem());
      comp_obj.setup(par.data(), wk_mem);

      double res{};
      for(vajoint_uint i = 0; i < comp_obj.n_obs(); ++i)
        res += comp_obj(par.data(), wk_mem, i);
      expect_true(pass_rel_err(res, true_val, 1e-7));
    }

    // test the gradient
    cfaad::Number::tape->rewind();
    std::vector<cfaad::Number> ad_par(par.size());
    cfaad::convertCollection(par.begin(), par.end(), ad_par.begin());


    cfaad::Number * wk_mem = wmem::get_Number_mem(comp_obj.n_wmem());

    cfaad::Number res{0};
    for(vajoint_uint i = 0; i < comp_obj.n_obs(); ++i)
      res += comp_obj(ad_par.data(), wk_mem, i);

    res.propagateToStart();
    expect_true(pass_rel_err(res.value(), true_val, 1e-6));

    expect_true(ad_par.size() == static_cast<size_t>(
      std::distance(begin(true_derivs), end(true_derivs))));

    for(size_t i = 0; i < ad_par.size(); ++i)
      expect_true(pass_rel_err(ad_par[i].adjoint(), true_derivs[i], 1e-6));

    // clean up
    wmem::clear_all();
  }
}
