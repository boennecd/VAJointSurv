test_that("manual pages gives the same results as previously", {
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

  # TODO: do a fit
  # TODO: add the Hessian matrix
  # TODO: ...
})
