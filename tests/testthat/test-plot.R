test_that(desc = "The grid is drawn before the curves without padding the axes", code = {
  grDevices::pdf(file = NULL)
  on.exit(expr = grDevices::dev.off())

  layers <- character()
  record_layer <- function(layer) layers <<- c(layers, layer)
  plot_env <- new.env(parent = environment(fun = plot_marker))
  plot_env$grid <- function() {
    record_layer(layer = "grid")
    graphics::grid()
  }
  plot_env$polygon <- function(...) {
    record_layer(layer = "band")
    graphics::polygon(...)
  }
  draw_marker <- plot_marker
  draw_surv <- plot_surv
  environment(draw_marker) <- plot_env
  environment(draw_surv) <- plot_env

  fixed <- poly_term(degree = 1, raw = TRUE)
  random <- poly_term(degree = 0, raw = TRUE, intercept = TRUE)
  time_range <- c(0, 5)
  vcov <- matrix(data = 1)

  draw_marker(
    time_fixef = fixed, time_rng = random, fixef_vary = .1,
    x_range = time_range, vcov_vary = vcov,
    panel.first = record_layer(layer = "before curve"),
    panel.last = record_layer(layer = "after curve"))
  expect_identical(object = layers,
                   expected = c("grid", "band", "before curve", "after curve"))
  expect_equal(object = graphics::par("usr")[1:2], expected = time_range)

  layers <- character()
  draw_surv(
    time_fixef = fixed, time_rng = list(random), fixef_vary = .1,
    x_range = time_range, vcov_vary = vcov, frailty_var = 0, associations = 1,
    panel.first = record_layer(layer = "before curve"),
    panel.last = record_layer(layer = "after curve"))
  expect_identical(object = layers, expected = c("grid", "before curve", "after curve"))
  expect_equal(object = graphics::par("usr")[1:3], expected = c(time_range, 0))
})
