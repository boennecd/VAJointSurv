test_that("marker_term fails with singular design matrices", {
  ti <- 1:10
  y <- numeric(length(ti))

  expect_error(
    marker_term(y ~ 1, id = ti, time_fixef = poly_term(ti, intercept = TRUE),
                time_rng = poly_term(ti, intercept = TRUE)),
    regexp = "Design matrix does not have full rank. Perhaps remove an intercept or a time-varying term from 'formula'")

  expect_error(
    marker_term(y ~ ti, id = ti, time_fixef = poly_term(ti, intercept = FALSE),
                time_rng = poly_term(ti, intercept = TRUE)),
    regexp = "Design matrix does not have full rank. Perhaps remove an intercept or a time-varying term from 'formula'")
})
