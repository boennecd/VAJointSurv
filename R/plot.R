#' Plots a Markers Mean Curve with Pointwise Quantiles
#'
#' @inheritParams marker_term
#' @param fixef_vary fixed effect coefficients for \code{time_fixef}.
#' @param x_range 2D numeric vector with start and end points.
#' @param vcov_vary the covariance matrix for \code{time_rng}.
#' @param p coverage of the two quantiles.
#' @param xlab,ylab,... arguments passed to \code{\link{plot}}.
#'
#' @importFrom stats qnorm
#' @importFrom graphics polygon grid
#' @importFrom grDevices gray
#' @export
plot_marker <- function(time_fixef, time_rng, fixef_vary, x_range, vcov_vary,
                        p = .95, xlab = "Time", ylab = "Marker", ...){
  is_valid_expansion(time_fixef)
  is_valid_expansion(time_rng)
  stopifnot(length(x_range) == 2, is.numeric(x_range), all(is.finite(x_range)),
            diff(x_range) > 0,
            length(p) == 1, is.finite(p), p > 0, p < 1)

  xs <- seq(x_range[1], x_range[2], length.out = 200)
  mea <- drop(fixef_vary %*% time_fixef$eval(xs))
  M <- time_rng$eval(xs)
  sds <- sqrt(diag(crossprod(M, vcov_vary %*% M))) # very inefficient

  sds <- sds * qnorm((1 + p) / 2)
  lbs <- mea - sds
  ubs <- mea + sds
  plot(xs, mea, ylim = range(lbs, ubs), type = "l", bty = "l", xlab = xlab,
       ylab = ylab, xaxs = "i", ...)
  polygon(x = c(xs, rev(xs)), y = c(lbs, rev(ubs)), border = NA,
          col = gray(0, .1))
  grid()

  invisible(list(lbs = lbs, ubs = ubs, mea = mea))
}

#' Plots Quantilies of the Conditional Hazards
#' @export
plot_surv <- function(time_fixef, time_rng, x_range, fixef_vary, vcov_vary,
                      frailty_var, ps = c(.025, .5, .975), log_hazard_shift = 0,
                      associations, xlab = "Time", ylab = "Hazard",...){
  # checks
  is_valid_expansion(time_fixef)
  if(is.list(time_rng))
    for(i in seq_along(time_rng))
      is_valid_expansion(time_rng[[i]])
  else {
    is_valid_expansion(time_rng)
    time_rng <- list(time_rng)
  }
  stopifnot(
    is.numeric(ps), all(ps > 0), all(ps < 1), all(is.finite(ps)),
    is.numeric(x_range), all(x_range >= 0), all(is.finite(x_range)),
    length(x_range) == 2, length(frailty_var) == 1)
  x_range <- sort(x_range)

  # assign function to evaluate the hazard pointwise
  time_rngs <- function(x){
    bases <- lapply(time_rng, function(expansion) expansion$eval(x))
    do.call(c, mapply(`*`, bases, associations, SIMPLIFY = FALSE))
  }

  tis <- seq(x_range[1], x_range[2], length.out = 100)
  hazs <- t(sapply(tis, function(ti){
    log_haz <- log_hazard_shift + fixef_vary %*% time_fixef$eval(ti)
    ti_basis <- time_rngs(ti)
    log_haz_var <- drop(frailty_var) + ti_basis %*% vcov_vary %*% ti_basis

    exp(qnorm(p = ps, mean = log_haz, sd = sqrt(log_haz_var)))
  }))
  matplot(tis, hazs, lty = 1, type = "l", col = "black", bty = "l",
          xlab = xlab, xaxs = "i", yaxs = "i", ylab = ylab,
          ylim = range(hazs, 0))
  grid()

  invisible(list(time = tis, hazard = hazs))
}
