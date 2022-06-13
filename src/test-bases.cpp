#include "testthat-wrapper.h"
#include "bases.h"
#include <array>

template<class Basis, size_t N>
void run_test(double const xx_val, std::array<double, N> const &yy_val,
              std::array<double, N> const &dx_val, bool const intercept,
              std::array<double, N> const ix_val = {},
              bool const comp_eval = false, double const xx_lower = 0){
  arma::vec bk = { 0., 1. },
            ik = { 0.333333333333333, 0.666666666666667 };
  int const order(4);

  auto bas = Basis(bk, ik, intercept, order);
  expect_true(bas.n_weights() == 0);
  arma::vec y = bas(xx_val, wmem::get_double_mem(bas.n_wmem()), nullptr, 0);

  expect_true(y.size() == yy_val.size());
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(y[i], yy_val[i]));

  arma::vec dx = bas(xx_val, wmem::get_double_mem(bas.n_wmem()), nullptr, 1);
  expect_true(dx.size() == dx_val.size());
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(dx[i], dx_val[i]));

  // work when a pointer is passed
  y.zeros();
  bas(y.memptr(), wmem::get_double_mem(bas.n_wmem()), xx_val, nullptr, 0);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(y[i], yy_val[i]));

  dx.zeros();
  bas(dx.memptr(), wmem::get_double_mem(bas.n_wmem()), xx_val, nullptr, 1);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(dx[i], dx_val[i]));

  if(!comp_eval)
    return;

  // test the integral
  bas.set_lower_limit(xx_lower);
  arma::vec ix = bas(xx_val, wmem::get_double_mem(bas.n_wmem()), nullptr, -1);
  expect_true(ix.size() == ix_val.size());
  for(unsigned i = 0; i < ix.size(); ++i)
    expect_true(pass_rel_err(ix[i], ix_val[i], 1e-6));
}

template<class Basis, size_t N>
void run_test_use_log(double const xx_val, std::array<double, N> const &yy_val,
                      std::array<double, N> const &dx_val, bool const intercept){
  arma::vec bk = { 1, 3 },
            ik = { 1.5, 2.5 };
  int const order(4);

  auto bas = Basis(bk, ik, intercept, order, true);
  arma::vec y = bas(xx_val, wmem::get_double_mem(bas.n_wmem()), nullptr, 0);

  expect_true(y.size() == yy_val.size());
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(y[i], yy_val[i]));

  arma::vec dx = bas(xx_val, wmem::get_double_mem(bas.n_wmem()), nullptr, 1);
  expect_true(dx.size() == dx_val.size());
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(dx[i], dx_val[i]));

  // work when a pointer is passed
  y.zeros();
  bas(y.memptr(), wmem::get_double_mem(bas.n_wmem()), xx_val, nullptr, 0);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(y[i], yy_val[i]));

  dx.zeros();
  bas(dx.memptr(), wmem::get_double_mem(bas.n_wmem()), xx_val, nullptr, 1);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(dx[i], dx_val[i]));
}

context("test bs") {
  test_that("bs works (no intercept)") {
    /*
     library(splines)
     library(numDeriv)
     xs <- c(-1, .5, 2)
     dput(interior_knots <- (1:2)/3)
     dput(boundary_knots <- c(0.0, 1.0))
     f <- function(z)
     bs(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = FALSE)
     dput(t(f(xs)))
     dput(sapply(xs, function(zz) jacobian(f, zz)))
     f <- function(z, i)
     bs(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = FALSE)[, i]
     dput(sapply(1:5, function(i) integrate(f, 0, xs[1], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:5, function(i) integrate(f, 0, xs[2], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:5, function(i) integrate(f, 0, xs[3], i = i, rel.tol = 1e-12)$value))
     */
    bool const intercept(false);
    double xx_val(-1);
    std::array<double, 5> yy_val{-96.75, 38.25, -4.5, 0, 0},
                          dx_val{231.749999999124, -101.249999999758,
                                 13.4999999999143, 0, 0},
                          ix_val{29.8125, -10.6875, 1.125, 0, 0};
    run_test<joint_bases::bs>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = .5;
    yy_val = {0.03125, 0.46875, 0.46875, 0.03125, 0};
    dx_val = {-0.56249999999794, -1.68749999999074, 1.68749999998975,
              0.56249999999794, 0};
    ix_val = {0.165364583333333, 0.18359375, 0.06640625, 0.00130208333333333, 0};
    run_test<joint_bases::bs>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = 2;
    yy_val = {0, -4.5, 38.25, -96.75, 64};
    dx_val = {0, -13.4999999999394, 101.24999999954, -231.749999999544,
              143.999999999226};
    ix_val = {0.166666666666667, -0.875, 10.9375, -29.6458333333333, 21.3333333333333};
    run_test<joint_bases::bs>(xx_val, yy_val, dx_val, intercept, ix_val, true);
  }

  test_that("bs works (intercept)") {
    /*
     library(splines)
     library(numDeriv)
     xs <- c(-1, .5, 2)
     dput(interior_knots <- (1:2)/3)
     dput(boundary_knots <- c(0.0, 1.0))
     f <- function(z)
     bs(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = TRUE)
     dput(t(f(xs)))
     dput(sapply(xs, function(zz) jacobian(f, zz)))
     f <- function(z, i)
     bs(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = TRUE)[, i]
     dput(sapply(1:6, function(i) integrate(f, 0, xs[1], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:6, function(i) integrate(f, 0, xs[2], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:6, function(i) integrate(f, 0, xs[3], i = i, rel.tol = 1e-12)$value))
     */

    bool constexpr intercept(true);
    double xx_val(-1);
    std::array<double, 6> yy_val,
                          dx_val,
                          ix_val;
    yy_val = { 64, -96.75, 38.25, -4.5, 0, 0 };
    dx_val = {-143.999999999379, 231.749999999124, -101.249999999758,
              13.4999999999143, 0, 0 };
    ix_val = { -21.25, 29.8125, -10.6875, 1.125, 0, 0 };
    run_test<joint_bases::bs>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = .5;
    yy_val = { 0, 0.03125, 0.46875, 0.46875, 0.03125, 0};
    dx_val = { 0, -0.56249999999794, -1.68749999999074,
               1.68749999998975, 0.56249999999794, 0 };
    ix_val = { 0.0833333333333333, 0.165364583333333, 0.18359375, 0.06640625,
               0.00130208333333333, 0 };
    run_test<joint_bases::bs>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = 2;
    yy_val = { 0, 0, -4.5, 38.25, -96.75, 64 };
    dx_val = { 0, 0, -13.4999999999394,
               101.24999999954, -231.749999999544, 143.999999999226 };
    ix_val = { 0.0833333333333333, 0.166666666666667, -0.875, 10.9375, -29.6458333333333, 21.3333333333333 };
    run_test<joint_bases::bs>(xx_val, yy_val, dx_val, intercept, ix_val, true);
  }

  test_that("bs works (intercept) use_log = true") {
    /*
     library(splines)
     library(numDeriv)
     xs <- exp(c(.5, 1.9, 3.2))
     dput(interior_knots <- c(1.5, 2.5))
     dput(boundary_knots <- c(1, 3))
     f <- function(z)
     bs(log(z), knots = interior_knots, Boundary.knots = boundary_knots, intercept = TRUE)
     dput(t(f(xs)))
     dput(sapply(xs, function(zz) jacobian(f, zz)))
     */

    bool constexpr intercept(true);
    double xx_val(std::exp(.5));
    std::array<double, 6> yy_val,
                          dx_val,
                          ix_val;
    yy_val = { 8, -8.44444444444444, 1.52777777777778, -0.0833333333333333,
               0, 0 };
    dx_val = {-14.5567358332459, 18.6002735647145, -4.34680306132608,
              0.303265329860015, 0, 0 };
    run_test_use_log<joint_bases::bs>(xx_val, yy_val, dx_val, intercept);

    xx_val = std::exp(1.9);
    yy_val = { 0, 0.096, 0.503333333333333, 0.372222222222222, 0.0284444444444444,
               0 };
    dx_val = { 0, -0.0717929372264858, -0.0747843096103767,
               0.114669274737316, 0.0319079721007659, 0 };
    run_test_use_log<joint_bases::bs>(xx_val, yy_val, dx_val, intercept);

    xx_val = std::exp(3.2);
    yy_val = { 0, 0, -0.00533333333333335, 0.193777777777778, -1.93244444444445,
               2.744 };
    dx_val = { 0, 0, -0.00326097631834105,
               0.085872376383374, -0.561974918863904, 0.479363518798855 };
    run_test_use_log<joint_bases::bs>(xx_val, yy_val, dx_val, intercept);
  }

  test_that("bs works (intercept), ders == -1, one inner knot") {
    // this has failed failed before

    /* R code to reproduce the result
     do_int <- function(xs, f, lower_limit){
     n_terms <- length(f(xs[1]))
     g <- function(x, i) f(x)[, i]
     outer(xs, 1:n_terms, Vectorize(
     function(x, i) integrate(g, lower_limit, x, i = i, rel.tol = 1e-12)$value))
     }
     bk <- c(0.00110535163162413, 2.40161776050478)
     ik <- 0.593923754422965

     f <- function(x)
     splines::bs(x, knots = ik, Boundary.knots = bk, intercept = TRUE)
     x <- 1.1
     dput(jacobian(f, x))
     dput(do_int(x, f, bk[1]))
     dput(do_int(x, f, .5))
     */

   arma::vec const bk{0.00110535163162413, 2.40161776050478},
                   ik{0.593923754422965};
    constexpr double derivs[]{0, -0.487927569566309, -0.175929182016681, 0.533786372465602, 0.130070379118054},
                     int_1[]{0.148204600697835, 0.531240503932728, 0.321197390710507, 0.0954760916404746, 0.00277606138691671},
                     int_2[]{9.33849208982431e-05, 0.244822455236761, 0.261365633984753, 0.0909424644706684, 0.00277606138691876};

    // bk, ik, intercept, order
    joint_bases::bs bas(bk, ik, true, 4);

    arma::vec res = bas(1.1, wmem::get_double_mem(bas.n_wmem()), nullptr, 1);
    expect_true(res.size() == 5);
    for(vajoint_uint i = 0; i < res.size(); ++i)
      expect_true(pass_rel_err(res[i], derivs[i]));

    bas.set_lower_limit(bk[0]);
    res = bas(1.1, wmem::get_double_mem(bas.n_wmem()), nullptr, -1);
    expect_true(res.size() == 5);
    for(vajoint_uint i = 0; i < res.size(); ++i)
      expect_true(pass_rel_err(res[i], int_1[i], 1e-7));

    bas.set_lower_limit(.5);
    res = bas(1.1, wmem::get_double_mem(bas.n_wmem()), nullptr, -1);
    expect_true(res.size() == 5);
    for(vajoint_uint i = 0; i < res.size(); ++i)
      expect_true(pass_rel_err(res[i], int_2[i], 1e-7));
  }
}

context("test ns") {
  test_that("ns works (no intercept)") {
    /*
     library(splines)
     library(numDeriv)
     xs <- c(-1, .5, 2)
     dput(interior_knots <- (1:2)/3)
     dput(boundary_knots <- c(0.0, 1.0))
     f <- function(z)
     ns(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = FALSE)
     dput(t(f(xs)))
     dput(sapply(xs, function(zz) jacobian(f, zz)))
     f <- function(z, i)
     ns(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = FALSE)[, i]
     dput(sapply(1:3, function(i) integrate(f, 0, xs[1], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:3, function(i) integrate(f, 0, xs[2], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:3, function(i) integrate(f, 0, xs[3], i = i, rel.tol = 1e-12)$value))
     */

    bool constexpr intercept(false);
    double xx_val(-1);
    std::array<double, 3> yy_val,
                          dx_val,
                          ix_val;
    yy_val = { 0.760638829255665, -2.28191648776699, 1.52127765851133 };
    dx_val = { -0.760638829254791, 2.28191648776756, -1.52127765850958 };
    ix_val = {-0.380319414627832, 1.1409582438835, -0.760638829255665};
    run_test<joint_bases::ns>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = .5;
    yy_val = { 0.320473361597061, 0.476079915208815, -0.296553276805877 };
    dx_val = { 2.16289926827143, -0.863697804854082, 0.950798536564431 };
    ix_val = { 0.0014165944809377, 0.19627104989052, -0.129979311038125 };
    run_test<joint_bases::ns>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = 2;
    yy_val = { -3.35714285714286, 1.07142857142857, 3.28571428571428 };
    dx_val = { -3.21428571427524, 0.642857142843741, 2.57142857136958 };
    ix_val = { -1.57152009239773, 1.13122694385987, 1.94029314853787 };
    run_test<joint_bases::ns>(xx_val, yy_val, dx_val, intercept, ix_val, true);
  }

  test_that("ns works (intercept)") {
    /*
     library(splines)
     library(numDeriv)
     xs <- c(-1, .5, 2)
     dput(interior_knots <- (1:2)/3)
     dput(boundary_knots <- c(0.0, 1.0))
     f <- function(z)
     ns(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = TRUE)
     dput(t(f(xs)))
     dput(sapply(xs, function(zz) jacobian(f, zz)))
     f <- function(z, i)
     ns(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = TRUE)[, i]
     dput(sapply(1:4, function(i) integrate(f, 0, xs[1], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:4, function(i) integrate(f, 0, xs[2], i = i, rel.tol = 1e-12)$value))
     dput(sapply(1:4, function(i) integrate(f, 0, xs[3], i = i, rel.tol = 1e-12)$value))
     */

    bool constexpr intercept(true);
    double xx_val(-1);
    std::array<double, 4> yy_val,
                          dx_val,
                          ix_val;
    yy_val = {-3.92940171279715, -0.745179167993104, 2.23553750397931,
              -1.49035833598621};
    dx_val = {3.66214047089044, 0.530893453727897, -1.59268036116767,
              1.06178690745579};
    ix_val = { 2.09833147735479, 0.479732441139409, -1.43919732341823, 0.959464882278818 };
    run_test<joint_bases::ns>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = .5;
    yy_val = {0.451294593143432, 0.419616910760803, 0.17864926771759,
              -0.0982661784783934};
    dx_val = {-1.6874999999859, 1.8378344485655, 0.11149665426856,
              0.300668897151231};
    ix_val = {0.175868152221105, 0.0115521618050298, 0.165864347918313, -0.109708176389986};
    run_test<joint_bases::ns>(xx_val, yy_val, dx_val, intercept, ix_val, true);

    xx_val = 2;
    yy_val = {0, -3.35714285714286, 1.07142857142857,
              3.28571428571428};
    dx_val = {0, -3.21428571427524, 0.642857142843741, 2.57142857136958};
    ix_val = { 0.23936516774501, -1.54711799479435, 1.05802065104973, 1.98909734374463 };
    run_test<joint_bases::ns>(xx_val, yy_val, dx_val, intercept, ix_val, true);
  }

  test_that("ns works (intercept) use_log = true") {
    /*
     library(splines)
     library(numDeriv)
     xs <- exp(c(.5, 1.9, 3.2))
     dput(interior_knots <- c(1.5, 2.5))
     dput(boundary_knots <- c(1, 3))
     f <- function(z)
     ns(log(z), knots = interior_knots, Boundary.knots = boundary_knots, intercept = TRUE)
     dput(t(f(xs)))
     dput(sapply(xs, function(zz) jacobian(f, zz)))
     */

    bool constexpr intercept(true);
    double xx_val(std::exp(.5));
    std::array<double, 4> yy_val,
                          dx_val,
                          ix_val;
    yy_val = { -1.0750421367307, -0.254984030450728, 1.01993612180291,
               -0.764952091352183 };
    dx_val = { 1.06619113518837, 0.122686446044928, -0.490745784179713,
               0.368059338136731 };
    run_test_use_log<joint_bases::ns>(xx_val, yy_val, dx_val, intercept);

    xx_val = std::exp(1.9);
    yy_val = { 0.500443700571341, 0.341188087871783, 0.1525809818462,
               -0.0931024030513163 };
    dx_val = { -0.0799272289080084, 0.125213152356344, -0.0102675383769478,
               0.0316316328581655 };
    run_test_use_log<joint_bases::ns>(xx_val, yy_val, dx_val, intercept);

    xx_val = std::exp(3.2);
    yy_val = { 0, -0.438461538461539, 0.553846153846154,
               0.884615384615385 };
    dx_val = { 0, -0.0658466371974155, 0.0188133249128504,
               0.0470333122847401 };
    run_test_use_log<joint_bases::ns>(xx_val, yy_val, dx_val, intercept);
  }
}

context("test iSpline") {
  test_that("iSpline works (no intercept)") {
    /*
     library(splines2)
     library(numDeriv)
     xs <- c(-1, 0, .5, 1, 2)
     dput(interior_knots <- (1:2)/3)
     dput(boundary_knots <- c(0.0, 1.0))
     f <- function(z)
     iSpline(
     z, knots = interior_knots, Boundary.knots = boundary_knots,
     intercept = FALSE)
     for(x in xs){
     dput(x)
     dput(c(f(x)))
     dput(c(jacobian(f, x, side = -1, method.args=list(eps=1e-9))))
     dput(c(jacobian(f, x, side =  1, method.args=list(eps=1e-9))))
     dput(c(jacobian(f, x,         , method.args=list(eps=1e-9))))
     }
    */

    bool constexpr intercept(false);
    double xx_val(-1);
    std::array<double, 5> yy_val,
                          dx_val;
    /* TODO: this is what we should get outside (0, 1), right? */
    yy_val = { 0, 0, 0, 0, 0 };
    dx_val = { 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 0, 0, 0, 0, 0 };
    dx_val = { 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0.9921875, 0.734375, 0.265625, 0.0078125, 0 };
    dx_val = { 0.1875, 1.875, 1.875, 0.1875, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 12 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    /* TODO: this is what we should get outside (0, 1), right? */
    xx_val = 2;
    yy_val = { 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);
  }

  test_that("iSpline works (intercept)") {
    /*
    library(splines2)
    library(numDeriv)
    xs <- c(-1, 0, .5, 1, 2)
    dput(interior_knots <- (1:2)/3)
    dput(boundary_knots <- c(0.0, 1.0))
    f <- function(z)
    iSpline(
      z, knots = interior_knots, Boundary.knots = boundary_knots,
    intercept = TRUE)
    for(x in xs){
    dput(x)
    dput(c(f(x)))
    dput(c(jacobian(f, x, side = -1, method.args=list(eps=1e-9))))
    dput(c(jacobian(f, x, side =  1, method.args=list(eps=1e-9))))
    dput(c(jacobian(f, x,         , method.args=list(eps=1e-9))))
    }
    */

    bool constexpr intercept(true);
    double xx_val(-1);
    std::array<double, 6> yy_val,
                          dx_val;
    /* TODO: this is what we should get outside (0, 1), right? */
    yy_val = { 0, 0, 0, 0, 0, 0 };
    dx_val = { 0, 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 0, 0, 0, 0, 0, 0 };
    dx_val = { 12, 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 1, 0.9921875, 0.734375, 0.265625, 0.0078125, 0 };
    dx_val = { 0, 0.1875, 1.875, 1.875, 0.1875, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 1, 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 0, 12 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);

    /* TODO: this is what we should get outside (0, 1), right? */
    xx_val = 2;
    yy_val = { 1, 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline>(xx_val, yy_val, dx_val, intercept);
  }
}

context("test msSpline") {
  test_that("mSpline works (no intercept)") {
    /*
    library(splines2)
    library(numDeriv)
    xs <- c(-1, 0, .5, 1, 2)
    dput(interior_knots <- (1:2)/3)
    dput(boundary_knots <- c(0.0, 1.0))
    f <- function(z)
    mSpline(
      z, knots = interior_knots, Boundary.knots = boundary_knots,
    intercept = FALSE)
    for(x in xs){
    dput(x)
    dput(c(f(x)))
    dput(c(jacobian(f, x, side = -1, method.args=list(eps=1e-9))))
    dput(c(jacobian(f, x, side =  1, method.args=list(eps=1e-9))))
    dput(c(jacobian(f, x,         , method.args=list(eps=1e-9))))
    }
    */

    bool constexpr intercept(false);
    double xx_val(-1);
    std::array<double, 5> yy_val,
                          dx_val;
    yy_val = { -580.5, 153, -18, 0, 0 };
    dx_val = { 1390.5, -405, 54, 0, 0 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 0, 0, 0, 0, 0 };
    dx_val = { 54, 0, 0, 0, 0 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0.1875, 1.875, 1.875, 0.1875, 0 };
    dx_val = { -3.375, -6.75, 6.75, 3.375, 0 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 0, 0, 0, 0, 12 };
    dx_val = { 0, 0, 0, -54, 108 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = { 0, -18, 153, -580.5, 768 };
    dx_val = { 0, -54, 405, -1390.5, 1728 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);
  }

  test_that("mSpline works (intercept)") {
    /*
    library(splines2)
    library(numDeriv)
    xs <- c(-1, 0, .5, 1, 2)
    dput(interior_knots <- (1:2)/3)
    dput(boundary_knots <- c(0.0, 1.0))
    f <- function(z)
    mSpline(
      z, knots = interior_knots, Boundary.knots = boundary_knots,
    intercept = TRUE)
    for(x in xs){
    dput(x)
    dput(c(f(x)))
    dput(c(jacobian(f, x, side = -1, method.args=list(eps=1e-9))))
    dput(c(jacobian(f, x, side =  1, method.args=list(eps=1e-9))))
    dput(c(jacobian(f, x,         , method.args=list(eps=1e-9))))
    }
    */

    bool constexpr intercept(true);
    double xx_val(-1);
    std::array<double, 6> yy_val,
                          dx_val;
    yy_val = { 768, -580.5, 153, -18, 0, 0 };
    dx_val = { -1728, 1390.5, -405, 54, 0, 0 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 12, 0, 0, 0, 0, 0 };
    dx_val = { -108, 54, 0, 0, 0, 0 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0, 0.1875, 1.875, 1.875, 0.1875, 0 };
    dx_val = { 0, -3.375, -6.75, 6.75, 3.375, 0 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 0, 0, 0, 0, 0, 12 };
    dx_val = { 0, 0, 0, 0, -54, 108 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = { 0, 0, -18, 153, -580.5, 768 };
    dx_val = { 0, 0, -54, 405, -1390.5, 1728 };
    run_test<joint_bases::mSpline>(xx_val, yy_val, dx_val, intercept);
  }
}

context("test orth_poly") {
  test_that("orth_poly gives the same as poly in R") {
    /*
     set.seed(1)
     dput(x <- round(rnorm(6), 2))
     obj <- poly(x, degree = 3)
     dput(attr(obj, "coefs"))
     dput(cbind(1, predict(obj, x)))
    */
    arma::vec const x = { -0.63, 0.18, -0.84, 1.6, 0.33, -0.82 },
                alpha = { -0.03, 0.673718350183412, 0.388455148829439 },
                norm2 = { 0.166666666666667, 1, 0.745133333333333, 0.413676179913314,
                          0.0235062965408605 };
    arma::mat basis =
      { { 1, 1, 1, 1, 1, 1, -0.283764870210735, 0.0993177045737573,
          -0.383082574784492, 0.770894564072497, 0.170258922126441, -0.373623745777468,
          0.0235472844805746, -0.538774146147203, 0.305295082237885, 0.485387375757513,
          -0.551505552689155, 0.276049956360386, 0.78636696277879, 0.160707894213562,
          -0.375908167424548, 0.0573748191776611, -0.396941602521704, -0.23159990622376 } };
    basis.reshape(6L, 4L);

    arma::mat Xout;
    joint_bases::orth_poly const obj =
      joint_bases::orth_poly::poly_basis(x, 3L, Xout);

    expect_true(basis.n_cols == Xout.n_cols);
    expect_true(basis.n_rows == Xout.n_rows);
    for(unsigned j = 0; j < Xout.n_cols; ++j)
      for(unsigned i = 0; i < Xout.n_rows; ++i)
        expect_true(pass_rel_err(Xout.at(i, j), basis.at(i, j)));

    for(unsigned i = 0; i < Xout.n_rows; ++i){
      arma::vec const b = obj(x[i], wmem::get_double_mem(obj.n_wmem()),
                              nullptr);
      expect_true(b.n_elem == Xout.n_cols);
      for(unsigned j = 0; j < Xout.n_cols; ++j)
        expect_true(pass_rel_err(Xout.at(i, j), b[j]));
    }
  }

  test_that("orth_poly works with raw == true") {
    // with the intercept
    for(unsigned i = 1; i < 4; ++i){
      joint_bases::orth_poly const obj    {i, true},
                                   obj_log{i, true, true};

      constexpr double x{2};
      const double log_x{std::log(x)};
      arma::vec res = obj    (x, wmem::get_double_mem(obj    .n_wmem()),
                              nullptr),
            res_log = obj_log(x, wmem::get_double_mem(obj_log.n_wmem()),
                              nullptr);
      double true_val{1}, true_val_log{1};
      for(unsigned j = 0; j <= i; ++j){
        expect_true(pass_rel_err(res    [j], true_val    ));
        expect_true(pass_rel_err(res_log[j], true_val_log));
        true_val *= x;
        true_val_log *= log_x;
      }

      expect_true(obj    .n_basis() == i + 1);
      expect_true(obj_log.n_basis() == i + 1);
    }

    // without the intercept
    for(unsigned i = 1; i < 4; ++i){
      joint_bases::orth_poly const obj    {i, false},
                                   obj_log{i, false, true};

      constexpr double x{3};
      const double log_x{std::log(x)};
      arma::vec res = obj    (x, wmem::get_double_mem(obj    .n_wmem()),
                              nullptr),
            res_log = obj_log(x, wmem::get_double_mem(obj_log.n_wmem()),
                              nullptr);
      double true_val{x}, true_val_log{log_x};
      for(unsigned j = 0; j < i; ++j){
        expect_true(pass_rel_err(res    [j], true_val));
        expect_true(pass_rel_err(res_log[j], true_val_log));
        true_val *= x;
        true_val_log *= log_x;
      }

      expect_true(obj    .n_basis() == i);
      expect_true(obj_log.n_basis() == i);
    }
  }

  test_that("orth_poly works with raw == true and ders > 0") {
    constexpr double x{3.5},
                 log_x{1.25276296849537};
    constexpr double d1[] {0, 1, 2 * 3.5, 3 * 3.5 * 3.5, 4 * 3.5 * 3.5 * 3.5},
                     d2[] {0, 0,       2,   3 * 2 * 3.5, 4 * 3 * 3.5 * 3.5  },
                     d3[] {0, 0,       0,         3 * 2, 4 * 3 * 2 * 3.5    },
                     d4[] {0, 0,       0,             0, 4 * 3 * 2          },
                  d_log[] {0, 1/x, 2 * log_x / x, 3 * log_x * log_x / x, 4 * log_x * log_x * log_x / x};


    // with an intercept
    {
      joint_bases::orth_poly const obj{4, true},
                               obj_log{4, true, true};

      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 1);
      expect_true(res.size() == 5);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d1[i]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 2);
      expect_true(res.size() == 5);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d2[i]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 3);
      expect_true(res.size() == 5);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d3[i]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 4);
      expect_true(res.size() == 5);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d4[i]));

      res = obj_log(x, wmem::get_double_mem(obj_log.n_wmem()), nullptr, 1);
      expect_true(res.size() == 5);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d_log[i]));
    }

    // without an intercept
    {
      joint_bases::orth_poly const obj{4, false},
                               obj_log{4, false, true};

      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 1);
      expect_true(res.size() == 4);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d1[i + 1]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 2);
      expect_true(res.size() == 4);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d2[i + 1]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 3);
      expect_true(res.size() == 4);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d3[i + 1]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 4);
      expect_true(res.size() == 4);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d4[i + 1]));

      res = obj_log(x, wmem::get_double_mem(obj_log.n_wmem()), nullptr, 1);
      expect_true(res.size() == 4);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d_log[i + 1]));
    }

  }

  test_that("orth_poly works with raw == true and ders < 0") {
    constexpr double x{.9},
                 x_low{.1},
                    i1[] {.9         , .9 * .9 / 2           , .9 * .9 * .9 / 3           , .9 * .9 * .9 * .9 / 4},
                    i2[] {.9 * .9 / 2, .9 * .9 * .9 / (2 * 3), .9 * .9 * .9 * .9 / (3 * 4), .9 * .9 * .9 * .9 * .9 / (4 * 5) },
                i1_low[] {.1         , .1 * .1 / 2           , .1 * .1 * .1 / 3           , .1 * .1 * .1 * .1 / 4},
                i2_low[] {.1 * .1 / 2, .1 * .1 * .1 / (2 * 3), .1 * .1 * .1 * .1 / (3 * 4), .1 * .1 * .1 * .1 * .1 / (4 * 5)};

      { // with an intercept
        joint_bases::orth_poly obj{3, true};
        arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -1);
        expect_true(res.size() == 4);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i1[i]));

        res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -2);
        expect_true(res.size() == 4);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i2[i]));

        obj.set_lower_limit(x_low);
        res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -1);
        expect_true(res.size() == 4);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i1[i] - i1_low[i]));

        res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -2);
        expect_true(res.size() == 4);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i2[i] - i2_low[i]));

      }
      { // without an intercept
        joint_bases::orth_poly obj{3, false};
        arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -1);
        expect_true(res.size() == 3);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i1[i + 1]));

        res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -2);
        expect_true(res.size() == 3);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i2[i + 1]));

        obj.set_lower_limit(x_low);
        res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -1);
        expect_true(res.size() == 3);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i1[i + 1] - i1_low[i + 1]));

        res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -2);
        expect_true(res.size() == 3);
        for(vajoint_uint i = 0; i < res.size(); ++i)
          expect_true(pass_rel_err(res[i], i2[i + 1] - i2_low[i + 1]));
      }
  }

  test_that("orth_poly works with raw == false") {
    /*
     dput(obj <- poly((-2):3, degree = 4))
     dput(predict(obj, 1.5))
     */
    arma::vec alpha{0.5, 0.5, 0.5, 0.5},
              norm2{1, 6, 17.5, 37.3333333333333, 64.8, 82.2857142857144};
    constexpr double x{1.5},
               truth[]{0.239045721866879, -0.313688217214239, -0.503115294937452, -0.0797268810253839};
    constexpr size_t n_res{4};

    // with intercept
    {
      joint_bases::orth_poly const obj(alpha, norm2, true);
      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr);

      expect_true(res[0] == 1);
      for(size_t i = 0; i < n_res; ++i)
        expect_true(pass_rel_err(res[i + 1], truth[i]));

      expect_true(obj.n_basis() == n_res + 1);
    }
    // without an intercept
    {
      joint_bases::orth_poly const obj(alpha, norm2, false);
      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr);

      for(size_t i = 0; i < n_res; ++i)
        expect_true(pass_rel_err(res[i], truth[i]));

      expect_true(obj.n_basis() == n_res);
    }
  }

  test_that("orth_poly works with raw == false and ders > 0") {
    const arma::vec alpha{-0.3448548522, 0.658759503793261, 0.499538901313647},
                    norm2{1, 10, 6.76808675371663, 3.22319474219269, 1.18536899578627};
    constexpr double x{.75};

    constexpr double d1[]{0.384385539783045, 0.660657439345008, -0.694447299988832},
                     d2[]{-6.90283286798098e-13, 1.11400392933787, 2.63891808433639};

    { // without an intercept
      /* brute force computation
       alpha <- c(-0.3448548522, 0.658759503793261, 0.499538901313647)
       norm2 <- c(1, 10, 6.76808675371663, 3.22319474219269, 1.18536899578627)
       dput(numDeriv::jacobian(
       function(x) poly(x, coefs = list(alpha = alpha, norm2 = norm2),
       degree = 3),
       .75))
       dput(sapply(1:3, function(i)
       numDeriv::hessian(
       function(x) poly(x, coefs = list(alpha = alpha, norm2 = norm2),
       degree = 3)[i], .75)))
       */

      joint_bases::orth_poly const obj(alpha, norm2, false);
      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 1);
      expect_true(res.size() == 3);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d1[i]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 2);
      expect_true(res.size() == 3);
      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d2[i]));

    }
    { // with an intercept
      /* brute force computation
       alpha <- c(-0.3448548522, 0.658759503793261, 0.499538901313647)
       norm2 <- c(1, 10, 6.76808675371663, 3.22319474219269, 1.18536899578627)
       dput(numDeriv::jacobian(
       function(x) poly(x, coefs = list(alpha = alpha, norm2 = norm2),
       degree = 3),
       .75))
       dput(sapply(1:3, function(i)
       numDeriv::hessian(
       function(x) poly(x, coefs = list(alpha = alpha, norm2 = norm2),
       degree = 3)[i], .75)))
       */

      joint_bases::orth_poly const obj(alpha, norm2, true);

      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 1);
      expect_true(res.size() == 4);
      expect_true(res[0] == 0);
      for(vajoint_uint i = 1; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d1[i - 1]));

      res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, 2);
      expect_true(res.size() == 4);
      expect_true(res[0] == 0);
      for(vajoint_uint i = 1; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], d2[i - 1]));

    }
  }

  test_that("orth_poly works with raw == false and ders < 0") {
    /* R code to reproduce the result
     alpha <- c(-0.3448548522, 0.658759503793261, 0.499538901313647)
     norm2 <- c(1, 10, 6.76808675371663, 3.22319474219269, 1.18536899578627)
     f <- function(x, i)
     poly(x, coefs = list(alpha = alpha, norm2 = norm2),
     degree = 3)[, i]
     dput(sapply(1:3, function(i) integrate(f, 0, 2, i = i)$value))
    */
    const arma::vec alpha{-0.3448548522, 0.658759503793261, 0.499538901313647},
                    norm2{1, 10, 6.76808675371663, 3.22319474219269, 1.18536899578627};
    constexpr double x{2};

    constexpr double int_val[]{1.03388551658793, 0.128604618505322, -0.0379503122041634};

    { // no intercept
      joint_bases::orth_poly const obj(alpha, norm2, false);
      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -1);
      expect_true(res.size() == 3);

      for(vajoint_uint i = 0; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], int_val[i]));
    }
    { // with intercept
      joint_bases::orth_poly const obj(alpha, norm2, true);
      arma::vec res = obj(x, wmem::get_double_mem(obj.n_wmem()), nullptr, -1);
      expect_true(res.size() == 4);

      expect_true(res[0] == x);
      for(vajoint_uint i = 1; i < res.size(); ++i)
        expect_true(pass_rel_err(res[i], int_val[i - 1]));
    }
  }
}

context("testing weighted basis"){
  /*
   library(splines)
   library(numDeriv)
   xs <- 2
   dput(interior_knots <- c(1,3))
   dput(boundary_knots <- c(0, 5))
   f <- function(z)
   ns(z, knots = interior_knots, Boundary.knots = boundary_knots, intercept = FALSE)
   dput(t(f(xs)))
   dput(sapply(xs, function(zz) jacobian(f, zz)))
   */
  const arma::vec ik{1,3};
  const arma::vec bk{0,5} ;
  const double x{2};
  const arma::vec basis_at_x{0.214240418913762, 0.519778743258714, -0.325685828839143};
  const arma::vec jacobian{0.40546301825183, -0.128889054761228, 0.148426036506845};

  test_that("single weighted"){
    auto do_tests = [&](joint_bases::basisMixin &test_basis){
      expect_true(test_basis.n_basis()==3);
      expect_true(test_basis.n_weights()==1);

      std::vector<double> mem(test_basis.n_wmem());
      const double weight = 2;
      arma::vec out = test_basis(x,mem.data(), &weight, 0);
      expect_true(out.n_elem == 3);

      for(unsigned i=0; i<out.n_elem;i++) {
        expect_true(pass_rel_err(weight*basis_at_x[i],out[i]));
      }
    };

    joint_bases::weighted_basis<joint_bases::ns> weight_1(bk, ik, false);
    do_tests(weight_1);
    auto weight_1_clone = weight_1.clone();
    do_tests(*weight_1_clone);
  }

  test_that("weighted of weighted"){
    auto do_test = [&](const joint_bases::basisMixin &test_basis){
      expect_true(test_basis.n_basis()==3);
      expect_true(test_basis.n_weights()==2);

      std::vector<double> mem(test_basis.n_wmem());

      std::array<double,2> weights{2,4};

      arma::vec out = test_basis(x,mem.data(), weights.data(), 0);
      expect_true(out.n_elem == 3);

      for(unsigned i=0; i<out.n_elem; ++i) {
        expect_true(pass_rel_err(out[i], weights[0]*weights[1]*basis_at_x[i]));
      }
    };

    joint_bases::weighted_basis
      <joint_bases::weighted_basis<joint_bases::ns> > weight_1(bk, ik, false);

    do_test(weight_1);

    auto weight_1_clone = weight_1.clone();
    do_test(*weight_1_clone);
  }
}

context("stacked basis") {
  /*
   library(splines)
   library(numDeriv)
   xs <- 2
   dput(interior_knots <- c(1,3))
   dput(boundary_knots <- c(0, 5))
   f <- function(z)
   cbind(
   ns(z, knots = interior_knots, Boundary.knots = boundary_knots,
   intercept = TRUE),
   1, z, z^2, z^3)

   dput(t(f(xs)))
   g <- function(x, i) f(x)[, i]
   dput(sapply(1:8, function(i) integrate(g, 1, xs, i = i)$value))
 */
  const arma::vec ik{1,3};
  const arma::vec bk{0,5} ;
  const double x{2};
  const arma::vec basis_at_x{0.528510354544419, 0.333871840764373, 0.160884477706881,
                             -0.0864229851379206, 1, 2, 4, 8};
  const arma::vec integral_from_one_to_x
    {0.579222604083723, 0.148606548120739, 0.194805355637782, -0.124661903758521,
     1, 1.5, 2.33333333333333, 3.75};

  test_that("properly stacked"){
    auto do_test = [&](joint_bases::basisMixin &test_basis){
      expect_true(test_basis.n_weights()==2);
      expect_true(test_basis.n_basis()==8);

      arma::vec const weights{3,5};
      std::vector<double> mem(test_basis.n_wmem());

      arma::vec out = test_basis(x, mem.data(),weights.memptr());
      expect_true(out.n_elem == 8);

      for (arma::uword i = 0; i <out.n_elem;i++) {
        expect_true(pass_rel_err(out[i],
                                 weights[i > 3 ? 1 : 0] * basis_at_x[i]));
      }

      test_basis.set_lower_limit(1);
      test_basis(out.memptr(), mem.data(),x, weights.memptr(),-1);

      for (arma::uword i = 0; i <out.n_elem;i++) {
        expect_true
          (pass_rel_err(out[i],
                        weights[i > 3 ? 1 : 0] * integral_from_one_to_x[i]));
      }
    };

    std::vector<std::unique_ptr<joint_bases::basisMixin> > input_arg;
    input_arg.emplace_back
      (new joint_bases::weighted_basis<joint_bases::ns>(bk, ik, true));
    input_arg.emplace_back
      (new joint_bases::weighted_basis<joint_bases::orth_poly>(3, true));

    joint_bases::stacked_basis test_basis(input_arg);
    auto tests_basis_clone = test_basis.clone();
    do_test(test_basis);
    do_test(*tests_basis_clone);
  }
}
