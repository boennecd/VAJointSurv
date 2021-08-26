#' Plots a Markers Mean Curve with Pointwise Quantiles
#' @importFrom stats qnorm
#' @importFrom graphics polygon grid
#' @export
plot_marker <- function(time_fixef, time_rng, fixef_vary, x_range, vcov_vary,
                        p = .95, xlab = "Time", ylab = "Marker", ...){
  is_valid_expansion(time_fixef)
  is_valid_expansion(time_rng)
  stopifnot(length(x_range) == 2, is.numeric(x_range), all(is.finite(x_range)),
            diff(x_range) > 0,
            length(p) == 1, is.finite(p), p > 0, p < 1)

  xs <- seq(x_range[1], x_range[2], length.out = 200)
  mea <- drop(fixef_vary %*% eval_expansion(time_fixef, xs))
  M <- eval_expansion(time_rng, xs)
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
