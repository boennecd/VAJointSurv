test_that(desc = "weighted terms work without delayed entries", code = {
  survival_data <- data.frame(
    id = 1:12,
    time = c(2, 6, 3, 4, 8, 3, 5, 7, 4, 5, 2, 6),
    event = rep(x = c(1, 0, 1, 0), times = 3L),
    weight = seq(from = 0.5, to = 1.5, length.out = 12L)
  )
  marker_data <- data.frame(
    id = rep(x = survival_data$id, each = 2L),
    time = rep(x = 0:1, times = 12L),
    y = c(3, 4, 2, 3, 5, 5, 3, 5, 2, 2, 4, 5,
          2, 4, 3, 3, 2, 3, 1, 3, 4, 4, 1, 2),
    weight = rep(x = survival_data$weight, each = 2L)
  )
  marker <- marker_term(
    formula = y ~ 1, id = id, data = marker_data,
    time_fixef = poly_term(x = time, degree = 1L, raw = TRUE),
    time_rng = weighted_term(
      x = poly_term(x = time, degree = 0L, raw = TRUE, intercept = TRUE),
      weight = weight)
  )
  terminal <- surv_term(
    formula = survival::Surv(time = time, event = event) ~ 1,
    id = id, data = survival_data,
    time_fixef = weighted_term(
      x = poly_term(x = time, degree = 1L, raw = TRUE), weight = weight)
  )

  expect_warning(
    object = model <- joint_ms_ptr(markers = marker, survival_terms = terminal),
    regexp = NA
  )
  bound <- joint_ms_lb(object = model, par = model$start_val)
  expect_true(object = is.finite(x = bound))
})
