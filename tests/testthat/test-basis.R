test_that("The C++ version of poly gives the right result", {
  in_x <- 2:5
  obj_truth <- poly(in_x, degree = 3)
  out_x <- 1:6
  truth <- predict(obj_truth, out_x)

  # without an intercept
  obj_cpp <- poly_term(in_x, degree = 3)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(truth), ignore_attr = TRUE)

  # with an intercept
  obj_cpp <- poly_term(in_x, degree = 3, intercept = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, rbind(1, t(truth)), ignore_attr = TRUE)

  # without an intercept and raw is TRUE
  obj_cpp <- poly_term(in_x, degree = 3, raw = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(outer(out_x, 1:3, `^`)))

  # with an intercept and raw is TRUE
  obj_cpp <- poly_term(in_x, degree = 3, raw = TRUE, intercept = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(outer(out_x, 0:3, `^`)))

  # without an intercept and degree == 0
  obj_cpp <- poly_term(in_x, degree = 0)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, matrix(0, 0, length(out_x)))

  # with an intercept and degree == 0
  obj_cpp <- poly_term(in_x, degree = 0, intercept = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, matrix(1, 1, length(out_x)))
})

test_that("The C++ version of bs gives the right result", {
  # without an intercept
  in_x <- 2:5
  obj_truth <- bs(in_x, df = 4)
  out_x <- c(2, 2.5, 3, 3.5, 4, 4.5, 5)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- bs_term(in_x, df = 4)
  expect_s3_class(obj_cpp, "bs_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)

  # with an intercept
  obj_truth <- bs(in_x, df = 4, intercept = TRUE)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- bs_term(in_x, df = 4, intercept = TRUE)
  expect_s3_class(obj_cpp, "bs_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)
})

test_that("The C++ version of ns gives the right result", {
  # without an intercept
  in_x <- 2:5
  obj_truth <- ns(in_x, df = 4)
  out_x <- c(2, 2.5, 3, 3.5, 4, 4.5, 5)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- ns_term(in_x, df = 4)
  expect_s3_class(obj_cpp, "ns_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)

  # with an intercept
  obj_truth <- ns(in_x, df = 4, intercept = TRUE)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- ns_term(in_x, df = 4, intercept = TRUE)
  expect_s3_class(obj_cpp, "ns_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- obj_cpp$eval(out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)
})

test_that("The plot_surv works with one-dimensional basis", {
  g1_basis <- ns_term(knots = c(3.33, 6.67), Boundary.knots = c(0, 10))
  g2_basis <- ns_term(knots = c(3.33, 6.67), Boundary.knots = c(0, 10))
  m1_basis <- poly_term(degree = 1, raw = TRUE, intercept = TRUE)
  m2_basis <- poly_term(degree = 0, raw = TRUE, intercept = TRUE)

  vcov_vary <- structure(c(0.021875, 5e-04, -0.003125, 5e-04, 0.012, -0.0015, -0.003125, -0.0015, 0.02), .Dim = c(3L, 3L))

  # the survival parameters
  fixef_surv <- c(-3, .4)
  association <- c(-1, 2)
  fixef_vary_surv <- c(.5, .1, -.015)
  fvar <- matrix(1e-6^2, 1)

  b_basis <- poly_term(degree = 3, raw = TRUE)

  # this failed in one version
  res <- plot_surv(
    time_fixef = b_basis, time_rng = list(m1_basis, m2_basis),
    x_range = c(0, 10), fixef_vary = fixef_vary_surv,  vcov_vary = vcov_vary,
    frailty_var = fvar, ps = c(.1, .5, .9), log_hazard_shift = fixef_surv[1],
    associations = association)

  expect_snapshot_value(res, cran = TRUE, style = "serialize")
})

test_that("A weighted poly_term gives the right resutls", {
  w_term <- weighted_term(poly_term(degree=3,raw = TRUE),x)
  expect_s3_class(w_term,"weighted_term")

  s <- c(3,4)
  x <- 2:3
  dat <- data.frame(x = x)

  expect_equal(w_term$eval(s,newdata = dat),
               matrix(c(s*x,s^2*x,s^3*x),ncol = length(s),byrow = TRUE))

  expect_error(weighted_term(w_term,2),
               "weighted_term of weighted_term is not supported")
  expect_error(weighted_term(1,x))
})

test_that("A stacked_term gives the right resutls", {
  main_term <- stacked_term(
    poly_term(degree=1,raw=TRUE),
    weighted_term(poly_term(degree=2,raw=TRUE),x),
    weighted_term(poly_term(degree=3,raw = TRUE),y),
    stacked_term(
      poly_term(degree=3,raw=TRUE),
      weighted_term(poly_term(degree=2,raw=TRUE),x)),
    stacked_term(
      poly_term(degree=3,raw=TRUE),
      weighted_term(
        stacked_term(
          poly_term(degree=3,raw=TRUE),
          weighted_term(poly_term(degree=2,raw=TRUE),x)),
        z)))

  s <- c(3,4)
  x <- 2:3
  y <- 1:2
  z <- c(10,10)

  correct_main_term <- matrix(
    c(s,s*x,s^2*x,s*y,s^2*y,s^3*y,s,s^2,s^3,s*x,s^2*x,s,s^2,s^3,s*z,s^2*z,s^3*z,s*x*z,s^2*x*z),
    ncol=length(s), byrow = TRUE)

  dat <- data.frame(x = x, y = y, z = z)

  expect_equal(main_term$eval(s,newdata = dat),correct_main_term)

  expect_s3_class(main_term, "stacked_term")

  expect_error(stacked_term(),
               "stacked_term created with less than two arguments")
  expect_error(stacked_term(1),
               "stacked_term created with less than two arguments")
  expect_error(main_term$eval(s,newdata = dat[1, ]))
  expect_error(main_term$eval(s,
                              newdata = transform(dat, x = as.character(x))))
  expect_error(stacked_term(1,2,3))
})
