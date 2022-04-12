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
#' @importFrom graphics polygon grid points
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
  plot(xs, mea, ylim = range(lbs, ubs), type = "n", xlab = xlab,
       ylab = ylab, bty = "l", ...)
  grid()
  polygon(x = c(xs, rev(xs)), y = c(lbs, rev(ubs)), border = NA,
          col = gray(0, .1))
  points(xs, mea, ylim = range(lbs, ubs), type = "l")

  invisible(list(lbs = lbs, ubs = ubs, mea = mea))
}

#' Plots Quantiles of the Conditional Hazards
#' @inheritParams surv_term
#'
#' @param time_rng an expansion or a list of expansions for the time-varying
#' random effects of the markers. See \code{\link{marker_term}}.
#' @param x_range two dimensional numerical vector with the range the hazard
#' should be plotted in.
#' @param fixef_vary fixed effect coefficients for \code{time_fixef}.
#' @param vcov_vary covariance matrix for the expansion or expansions in
#' \code{time_rng}.
#' @param frailty_var variance of the frailty.
#' @param ps quantiles to plot.
#' @param log_hazard_shift possible shift on the log hazard.
#' @param associations association parameter for each \code{time_rng} or
#' possible multiple parameters for each \code{time_rng} if \code{ders} is
#' supplied.
#' @param xlab,ylab,... arguments passed to \code{\link{matplot}}.
#'
#'
#' @importFrom graphics matplot
#' @export
plot_surv <- function(time_fixef, time_rng, x_range, fixef_vary, vcov_vary,
                      frailty_var, ps = c(.025, .5, .975), log_hazard_shift = 0,
                      associations, xlab = "Time", ylab = "Hazard",
                      ders = NULL, ...){
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

  if(is.null(ders))
    ders <- as.list(rep(0, length(time_rng)))

  # split the association parameters by effect
  association_var <- unlist(mapply(rep, seq_along(ders), lengths(ders)))
  associations <- split(associations, association_var)

  # assign function to evaluate the hazard pointwise
  time_rngs <- function(x){
    bases <- mapply(function(expansion, ders, assoc) {
      basis_vecs <- sapply(ders, expansion$eval, x = x)
      res <- basis_vecs * rep(assoc, each = NROW(basis_vecs))
      if(is.matrix(res)) rowSums(res) else res
    }, expansion = time_rng, ders = ders, assoc = associations,
    SIMPLIFY = FALSE)
    do.call(c, bases)
  }

  tis <- seq(x_range[1], x_range[2], length.out = 100)
  hazs <- t(sapply(tis, function(ti){
    log_haz <- log_hazard_shift + fixef_vary %*% time_fixef$eval(ti)
    ti_basis <- time_rngs(ti)
    log_haz_var <- drop(frailty_var) + ti_basis %*% vcov_vary %*% ti_basis

    exp(qnorm(p = ps, mean = log_haz, sd = sqrt(log_haz_var)))
  }))
  matplot(tis, hazs, type = "n", xlab = xlab, bty = "l", ylab = ylab,
          ylim = range(hazs, 0), ...)
  grid()
  matplot(tis, hazs, lty = 1, type = "l", col = "black", add = TRUE, ...)

  invisible(list(time = tis, hazard = hazs))
}
