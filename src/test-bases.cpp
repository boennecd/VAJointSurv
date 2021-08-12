#include "testthat-wrapper.h"
#include "bases.h"
#include <array>

template<class Basis, bool do_optimize, size_t N>
void run_test(double const xx_val, std::array<double, N> const &yy_val,
              std::array<double, N> const &dx_val, bool const intercept){
  arma::vec bk = { 0., 1. },
            ik = { 0.333333333333333, 0.666666666666667 };
  int const order(4);

  auto bas = Basis(bk, ik, intercept, order);
  arma::vec y = bas(xx_val, 0);

  expect_true(y .size() == yy_val.size());
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(y[i], yy_val[i]));

  arma::vec dx = bas(xx_val, 1);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(dx[i], dx_val[i]));

  // work when a pointer is passed
  y.zeros();
  bas(y.memptr(), xx_val, 0);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(y[i], yy_val[i]));

  dx.zeros();
  bas(dx.memptr(), xx_val, 1);
  for(unsigned i = 0; i < y.size(); ++i)
    expect_true(pass_rel_err(dx[i], dx_val[i]));
}

context("bases unit tests") {
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
     */
    bool const intercept(false);
    double xx_val(-1);
    std::array<double, 5> yy_val{-96.75, 38.25, -4.5, 0, 0},
                          dx_val{231.749999999124, -101.249999999758,
                                 13.4999999999143, 0, 0};
    run_test<joint_bases::bs, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::bs, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = {0.03125, 0.46875, 0.46875, 0.03125, 0};
    dx_val = {-0.56249999999794, -1.68749999999074, 1.68749999998975,
              0.56249999999794, 0};
    run_test<joint_bases::bs, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::bs, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = {0, -4.5, 38.25, -96.75, 64};
    dx_val = {0, -13.4999999999394, 101.24999999954, -231.749999999544,
              143.999999999226};
    run_test<joint_bases::bs, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::bs, true >(xx_val, yy_val, dx_val, intercept);
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
     */

    bool constexpr intercept(true);
    double xx_val(-1);
    std::array<double, 6> yy_val,
                          dx_val;
    yy_val = { 64, -96.75, 38.25, -4.5, 0, 0 };
    dx_val = {-143.999999999379, 231.749999999124, -101.249999999758,
              13.4999999999143, 0, 0 };
    run_test<joint_bases::bs, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::bs, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0, 0.03125, 0.46875, 0.46875, 0.03125, 0};
    dx_val = { 0, -0.56249999999794, -1.68749999999074,
               1.68749999998975, 0.56249999999794, 0 };
    run_test<joint_bases::bs, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::bs, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = { 0, 0, -4.5, 38.25, -96.75, 64 };
    dx_val = { 0, 0, -13.4999999999394,
               101.24999999954, -231.749999999544, 143.999999999226 };
    run_test<joint_bases::bs, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::bs, true >(xx_val, yy_val, dx_val, intercept);
  }

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
     */

    bool constexpr intercept(false);
    double xx_val(-1);
    std::array<double, 3> yy_val,
                          dx_val;
    yy_val = { 0.760638829255665, -2.28191648776699, 1.52127765851133 };
    dx_val = { -0.760638829254791, 2.28191648776756, -1.52127765850958 };
    run_test<joint_bases::ns, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::ns, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0.320473361597061, 0.476079915208815, -0.296553276805877 };
    dx_val = { 2.16289926827143, -0.863697804854082, 0.950798536564431 };
    run_test<joint_bases::ns, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::ns, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = { -3.35714285714286, 1.07142857142857, 3.28571428571428 };
    dx_val = { -3.21428571427524, 0.642857142843741, 2.57142857136958 };
    run_test<joint_bases::ns, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::ns, true >(xx_val, yy_val, dx_val, intercept);
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
     */

    bool constexpr intercept(true);
    double xx_val(-1);
    std::array<double, 4> yy_val,
                          dx_val;
    yy_val = {-3.92940171279715, -0.745179167993104, 2.23553750397931,
              -1.49035833598621};
    dx_val = {3.66214047089044, 0.530893453727897, -1.59268036116767,
              1.06178690745579};
    run_test<joint_bases::ns, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::ns, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = {0.451294593143432, 0.419616910760803, 0.17864926771759,
              -0.0982661784783934};
    dx_val = {-1.6874999999859, 1.8378344485655, 0.11149665426856,
              0.300668897151231};
    run_test<joint_bases::ns, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::ns, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = {0, -3.35714285714286, 1.07142857142857,
              3.28571428571428};
    dx_val = {0, -3.21428571427524, 0.642857142843741, 2.57142857136958};
    run_test<joint_bases::ns, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::ns, true >(xx_val, yy_val, dx_val, intercept);
  }

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
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 0, 0, 0, 0, 0 };
    dx_val = { 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0.9921875, 0.734375, 0.265625, 0.0078125, 0 };
    dx_val = { 0.1875, 1.875, 1.875, 0.1875, 0 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 12 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    /* TODO: this is what we should get outside (0, 1), right? */
    xx_val = 2;
    yy_val = { 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);
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
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 0, 0, 0, 0, 0, 0 };
    dx_val = { 12, 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 1, 0.9921875, 0.734375, 0.265625, 0.0078125, 0 };
    dx_val = { 0, 0.1875, 1.875, 1.875, 0.1875, 0 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 1, 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 0, 12 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);

    /* TODO: this is what we should get outside (0, 1), right? */
    xx_val = 2;
    yy_val = { 1, 1, 1, 1, 1, 1 };
    dx_val = { 0, 0, 0, 0, 0, 0 };
    run_test<joint_bases::iSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::iSpline, true >(xx_val, yy_val, dx_val, intercept);
  }

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
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 0, 0, 0, 0, 0 };
    dx_val = { 54, 0, 0, 0, 0 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0.1875, 1.875, 1.875, 0.1875, 0 };
    dx_val = { -3.375, -6.75, 6.75, 3.375, 0 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 0, 0, 0, 0, 12 };
    dx_val = { 0, 0, 0, -54, 108 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = { 0, -18, 153, -580.5, 768 };
    dx_val = { 0, -54, 405, -1390.5, 1728 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);
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
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 0;
    yy_val = { 12, 0, 0, 0, 0, 0 };
    dx_val = { -108, 54, 0, 0, 0, 0 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = .5;
    yy_val = { 0, 0.1875, 1.875, 1.875, 0.1875, 0 };
    dx_val = { 0, -3.375, -6.75, 6.75, 3.375, 0 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 1;
    yy_val = { 0, 0, 0, 0, 0, 12 };
    dx_val = { 0, 0, 0, 0, -54, 108 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);

    xx_val = 2;
    yy_val = { 0, 0, -18, 153, -580.5, 768 };
    dx_val = { 0, 0, -54, 405, -1390.5, 1728 };
    run_test<joint_bases::mSpline, false>(xx_val, yy_val, dx_val, intercept);
    run_test<joint_bases::mSpline, true >(xx_val, yy_val, dx_val, intercept);
  }

  test_that("orth_poly gives the same as poly in R") {
    /*
     set.seed(1)
     dput(x <- round(rnorm(6), 2))
     obj <- poly(x, degree = 3)
     attr(obj, "coefs")$norm2 <- attr(obj, "coefs")$norm2 / NROW(obj)
     dput(attr(obj, "coefs"))
     dput(cbind(1, predict(obj, x)))
    */
    arma::vec const x = { -0.63, 0.18, -0.84, 1.6, 0.33, -0.82 },
                alpha = { -0.03, 0.673718350183412, 0.388455148829439 },
                norm2 = { 0.166666666666667, 1, 0.745133333333333, 0.413676179913314,
                          0.0235062965408605 };
    arma::mat basis =
      { { 1, 1, 1, 1, 1, 1, -0.695079138943395, 0.243277698630188,
          -0.938356837573584, 1.88829832746289, 0.417047483366037, -0.915187532942137,
          0.0576788318055649, -1.31972174466434, 0.747817172463846, 1.18895139819447,
          -1.35090719440005, 0.676181536600509, 1.92619780939021, 0.393652338460405,
          -0.920783200334851, 0.14053903106972, -0.97230438386083, -0.567301594724647 } };
    basis.reshape(6L, 4L);

    arma::mat Xout;
    joint_bases::orth_poly const obj =
      joint_bases::orth_poly::get_poly_basis(x, 3L, Xout);

    expect_true(basis.n_cols == Xout.n_cols);
    expect_true(basis.n_rows == Xout.n_rows);
    for(unsigned j = 0; j < Xout.n_cols; ++j)
      for(unsigned i = 0; i < Xout.n_rows; ++i)
        expect_true(pass_rel_err(Xout.at(i, j), basis.at(i, j)));

    for(unsigned i = 0; i < Xout.n_rows; ++i){
      arma::vec const b = obj(x[i]);
      expect_true(b.n_elem == Xout.n_cols);
      for(unsigned j = 0; j < Xout.n_cols; ++j)
        expect_true(pass_rel_err(Xout.at(i, j), b[j]));
    }
  }

  test_that("orth_poly works with raw == true") {
    // with the intercept
    for(unsigned i = 1; i < 4; ++i){
      joint_bases::orth_poly const obj{i, true};

      constexpr double x{2};
      arma::vec res = obj(x);
      double true_val{1};
      for(unsigned j = 0; j <= i; ++j){
        expect_true(pass_rel_err(res[j], true_val));
        true_val *= x;
      }

      expect_true(obj.get_n_basis() == i + 1);
    }

    // without the intercept
    for(unsigned i = 1; i < 4; ++i){
      joint_bases::orth_poly const obj{i, false};

      constexpr double x{3};
      arma::vec res = obj(x);
      double true_val{x};
      for(unsigned j = 0; j < i; ++j){
        expect_true(pass_rel_err(res[j], true_val));
        true_val *= x;
      }

      expect_true(obj.get_n_basis() == i);
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
        arma::vec res = obj(x);

        expect_true(res[0] == 1);
        for(size_t i = 0; i < n_res; ++i)
          expect_true(pass_rel_err(res[i + 1], truth[i]));

        expect_true(obj.get_n_basis() == n_res + 1);
      }
      // without an intercept
      {
        joint_bases::orth_poly const obj(alpha, norm2, false);
        arma::vec res = obj(x);

        for(size_t i = 0; i < n_res; ++i)
          expect_true(pass_rel_err(res[i], truth[i]));

        expect_true(obj.get_n_basis() == n_res);
      }
    }
  }
}
