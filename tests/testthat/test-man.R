test_that("manual pages gives the same results as previously for joint_ms type functions", {
  # load in the data
  library(survival)
  data(pbc, package = "survival")

  # re-scale by year
  pbcseq <- transform(pbcseq, day_use = day / 365.25)
  pbc <- transform(pbc, time_use = time / 365.25)

  # create the marker terms
  m1 <- marker_term(
    log(bili) ~ 1, id = id, data = pbcseq,
    time_fixef = bs_term(day_use, df = 5L),
    time_rng = poly_term(day_use, degree = 1L, raw = TRUE, intercept = TRUE))
  m2 <- marker_term(
    albumin ~ 1, id = id, data = pbcseq,
    time_fixef = bs_term(day_use, df = 5L),
    time_rng = poly_term(day_use, degree = 1L, raw = TRUE, intercept = TRUE))

  # base knots on observed event times
  bs_term_knots <-
    with(pbc, quantile(time_use[status == 2], probs = seq(0, 1, by = .2)))

  boundary <- c(bs_term_knots[ c(1, length(bs_term_knots))])
  interior <- c(bs_term_knots[-c(1, length(bs_term_knots))])

  # create the survival term
  s_term <- surv_term(
    Surv(time_use, status == 2) ~ 1, id = id, data = pbc,
    time_fixef = bs_term(time_use, Boundary.knots = boundary, knots = interior))

  # create the C++ object to do the fitting
  model_ptr <- joint_ms_ptr(
    markers = list(m1, m2), survival_terms = s_term,
    max_threads = 2L, ders = list(0L, c(0L, -1L)))

  # find the starting values
  start_vals <- joint_ms_start_val(model_ptr)
  start_vals_to_test <- head(start_vals, 500)

  attributes(start_vals_to_test) <- attributes(start_vals)

  expect_snapshot_value(
    start_vals_to_test, style = "serialize",
    tolerance = 1e-3)

  expect_equal(attr(start_vals,"value"),
               joint_ms_lb(model_ptr,par = start_vals))

  expect_snapshot_value(joint_ms_format(model_ptr,start_vals),
                        style = "serialize",
                        tolerance = 1e-3)

  fit <- joint_ms_opt(object = model_ptr, par = start_vals, gr_tol = .1)
  hess <- joint_ms_hess(object = model_ptr,par = fit$par)

  expect_snapshot_value(fit[c("value", "info", "convergence")],
                        tolerance = 1e-5)


  expect_snapshot_value(joint_ms_format(model_ptr, fit$par),
                        style = "serialize", tolerance = 1e-3)

  expect_snapshot_value(hess$hessian,
                        style = "serialize", tolerance = 1e-3)
  skip_on_cran()
  se <- 0.131148235758747
  which_prof <- model_ptr$indices$survival[[1]]$associations[1]
  delta <- 2*se

  profile_CI <- joint_ms_profile(
    object = model_ptr, opt_out = fit, which_prof = which_prof,
    delta= delta, gr_tol = .1, verbose = FALSE)

  expect_snapshot_value(profile_CI[c("confs","xs","p_log_Lik")],
                        style = "json2",
                        tolerance = 1e-3)
})

test_that("test manual page example for bs_term", {
  vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
  spline_basis <- bs_term(vals,df = 3)

  library(splines)

    correct_basis <- bs(vals, df = 3)


  expect_equal(c(spline_basis$eval(0.5)),
                bs(0.5,Boundary.knots = attr(correct_basis, "Boundary.knots"),
                                           knots = attr(correct_basis,
                                                        "knots")),
               ignore_attr = TRUE)

  expect_equal(c(spline_basis$eval(0.5,der = 1)),
               c(-1.12, 0.853333333333333, 1.13777777777778))
})
