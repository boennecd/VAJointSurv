library(survival)
test_that("marker_term fails with singular design matrices", {
  ti <- 1:10
  y <- numeric(length(ti))
  y[1:3] <- 1
  dat <- data.frame(ti = ti, y = y)

  expect_error(
    surv_term(Surv(ti, y) ~ 1, id = ti, data = dat,
                time_fixef = poly_term(ti, intercept = TRUE)),
    regexp = "Design matrix does not have full rank. Perhaps remove an intercept or a time-varying term from 'formula'")
})
