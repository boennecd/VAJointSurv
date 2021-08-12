#include "marker-term.h"
#include "testthat-wrapper.h"
#include "wmem.h"
#include <iterator>

using std::begin;
using std::end;

context("marker_term is correct") {
  test_that("the number basic members are correct"){
    /** R code to reproduce the result
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
}
