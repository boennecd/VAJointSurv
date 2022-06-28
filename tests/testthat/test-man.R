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

  VA_pars <- joint_ms_va_par(object = model_ptr,par = start_vals)

  expect_equal(length(VA_pars),length(unique(pbc$id)))

  expect_snapshot_value(VA_pars[1:10],
                        style = "serialize",
                        tolerance = 1e-3)

  vcov_vary <- diag(1:4)

  alter_pars <- joint_ms_set_vcov(
    object = model_ptr,
    vcov_vary = vcov_vary,
    vcov_surv = matrix(0,0,0))

  for(vcov in joint_ms_va_par(object = model_ptr, par = alter_pars))
    expect_equal(vcov$vcov,vcov_vary)

  expect_equal(
    joint_ms_format(object = model_ptr, par = alter_pars)$vcov$vcov_vary,
    vcov_vary)

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
               c(-1.1199999999968, 0.853333333335853, 1.13777777777691),
               tolerance = 1e-6)
})

test_that("test manual page example for ns_term", {
  vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
  spline_basis <- ns_term(vals,df = 3)

  library(splines)

  correct_basis <- ns(vals, df = 3)


  expect_equal(c(spline_basis$eval(0.5)),
               ns(0.5,Boundary.knots = attr(correct_basis, "Boundary.knots"),
                  knots = attr(correct_basis,
                               "knots")),
               ignore_attr = TRUE)

  expect_equal(c(spline_basis$eval(0.5,der = 1)),
               c(1.06750746213182, -0.889635441325693, 1.95627432066475),
               tolerance = 1e-6)
})

test_that("plots marker returns same type of output", {
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

  fixef_vary <- c(-0.1048, 0.2583, 1.0578, 2.4006, 2.9734)
  vcov_vary <- rbind(c(0.96580, 0.09543), c(0.09543,  0.03998))

  # plot marker's trajectory
  marker_plot <- plot_marker(
    time_fixef = m1$time_fixef,
    time_rng = m1$time_rng,
    fixef_vary = fixef_vary,
    vcov_vary = vcov_vary, x_range = c(0,5))

  expect_snapshot_value(marker_plot, style = "serialize")
})
