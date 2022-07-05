#' Plots a Markers Mean Curve with Pointwise Quantiles
#'
#' @inheritParams marker_term
#' @param fixef_vary fixed effect coefficients for \code{time_fixef}.
#' @param x_range 2D numeric vector with start and end points.
#' @param vcov_vary the covariance matrix for \code{time_rng}.
#' @param p coverage of the two quantiles.
#' @param xlab,ylab,... arguments passed to \code{\link{plot}}.
#' @param newdata \code{data.frame} with data for the weights if any.
#'
#' @return
#' A list containing data for plotting.
#'
#' @importFrom stats qnorm
#' @importFrom graphics polygon grid
#' @importFrom grDevices gray
#' @examples
#' # load in the data
#' library(survival)
#' data(pbc, package = "survival")
#'
#' # re-scale by year
#' pbcseq <- transform(pbcseq, day_use = day / 365.25)
#' pbc <- transform(pbc, time_use = time / 365.25)
#'
#' # create the marker terms
#' m1 <- marker_term(
#'   log(bili) ~ 1, id = id, data = pbcseq,
#'   time_fixef = bs_term(day_use, df = 5L),
#'   time_rng = poly_term(day_use, degree = 1L, raw = TRUE, intercept = TRUE))
#'
#' fixef_vary <- c(-0.1048, 0.2583, 1.0578, 2.4006, 2.9734)
#' vcov_vary <- rbind(c(0.96580, 0.09543), c(0.09543,  0.03998))
#'
#' # plot marker's trajectory
#' plot_marker(
#'   time_fixef = m1$time_fixef,
#'   time_rng = m1$time_rng,
#'   fixef_vary = fixef_vary,
#'   vcov_vary = vcov_vary, x_range = c(0,5))
#' @export
plot_marker <- function(time_fixef, time_rng, fixef_vary, x_range, vcov_vary,
                        p = .95, xlab = "Time", ylab = "Marker",
                        newdata = NULL, ...){

  is_valid_expansion(time_fixef)
  is_valid_expansion(time_rng)
  stopifnot(length(x_range) == 2, is.numeric(x_range), all(is.finite(x_range)),
            diff(x_range) > 0,
            length(p) == 1, is.finite(p), p > 0, p < 1)

  xs <- seq(x_range[1], x_range[2], length.out = 200)
  mea <- drop(
    fixef_vary %*% extend_newdata_n_eval(time_fixef, xs, newdata = newdata))
  M <- extend_newdata_n_eval(time_rng, xs, newdata = newdata)
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
#' @param newdata \code{data.frame} with data for the weights if any.
#' @param ders a \code{\link{list}} with \code{\link{integer}} vectors for how
#' the survival outcome is linked to the markers. 0 implies present values,
#' -1 is integral of, and 1 is the derivative. \code{NULL} implies the present
#' value of the random effect for all markers.
#'
#' @return
#' A list containing data for plotting.
#'
#' @importFrom graphics matplot
#' @examples
#' # load in the data
#' library(survival)
#' data(pbc, package = "survival")
#'
#' # re-scale by year
#' pbcseq <- transform(pbcseq, day_use = day / 365.25)
#' pbc <- transform(pbc, time_use = time / 365.25)
#'
#' # create the marker terms
#' m1 <- marker_term(
#'   log(bili) ~ 1, id = id, data = pbcseq,
#'   time_fixef = bs_term(day_use, df = 5L),
#'   time_rng = poly_term(day_use, degree = 1L, raw = TRUE, intercept = TRUE))
#' m2 <- marker_term(
#'   albumin ~ 1, id = id, data = pbcseq,
#'   time_fixef = bs_term(day_use, df = 5L),
#'   time_rng = poly_term(day_use, degree = 1L, raw = TRUE, intercept = TRUE))
#'
#' # base knots on observed event times
#' bs_term_knots <-
#'   with(pbc, quantile(time_use[status == 2], probs = seq(0, 1, by = .2)))
#'
#' boundary <- c(bs_term_knots[ c(1, length(bs_term_knots))])
#' interior <- c(bs_term_knots[-c(1, length(bs_term_knots))])
#'
#' # create the survival term
#' s_term <- surv_term(
#'   Surv(time_use, status == 2) ~ 1, id = id, data = pbc,
#'   time_fixef = bs_term(time_use, Boundary.knots = boundary, knots = interior))
#'
#' # expansion of time for the fixed effects in the survival term
#' time_fixef <- s_term$time_fixef
#' # expansion of time for the random effects in the marker terms
#' time_rng <- list(m1$time_rng, m2$time_rng)
#' # no frailty
#' frailty_var <- matrix(0L,1)
#' # var-covar matrix for time-varying random effects
#' vcov_vary <- c(0.9658, 0.0954, -0.1756, -0.0418, 0.0954, 0.04, -0.0276,
#'                -0.0128, -0.1756, -0.0276, 0.1189, 0.0077, -0.0418, -0.0128,
#'                0.0077, 0.0057) |> matrix(4L)
#' # coefficients for time-varying fixed effects
#' fixef_vary <- c(1.0495, -0.2004, 1.4167, 1.255, 2.5007, 4.8545, 4.7889)
#' # association parameters
#' associations <- c(0.8627, -3.2358, 0.1842)
#' # constant shift on the log-hazard scale
#' log_hazard_shift <- -4.498513
#' # specify how the survival outcome is linked with markers
#' ders = list(0L, c(0L, -1L))
#'
#' # plot the hazard with pointwise quantiles
#' plot_surv(
#'   time_fixef = time_fixef,
#'   time_rng = time_rng,
#'   x_range = c(0, 5), vcov_vary = vcov_vary, frailty_var = frailty_var,
#'   ps = c(.25, .5, .75), log_hazard_shift = log_hazard_shift,
#'   fixef_vary = fixef_vary, associations = associations, ders = ders)
#' @export
plot_surv <- function(time_fixef, time_rng, x_range, fixef_vary, vcov_vary,
                      frailty_var, ps = c(.025, .5, .975), log_hazard_shift = 0,
                      associations, xlab = "Time", ylab = "Hazard",
                      ders = NULL, newdata = NULL, ...){
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
      basis_vecs <- sapply(ders, function(der){
        extend_newdata_n_eval(
          term = expansion, x = x, newdata = newdata, der = der)
      })
      res <- basis_vecs * rep(assoc, each = NROW(basis_vecs))
      if(is.matrix(res)) rowSums(res) else res
    }, expansion = time_rng, ders = ders, assoc = associations,
    SIMPLIFY = FALSE)
    do.call(c, bases)
  }

  tis <- seq(x_range[1], x_range[2], length.out = 100)
  hazs <- t(sapply(tis, function(ti){
    log_haz <- log_hazard_shift + fixef_vary %*%
      extend_newdata_n_eval(time_fixef, ti, newdata=newdata)
    ti_basis <- time_rngs(ti)
    log_haz_var <- drop(frailty_var) + ti_basis %*% vcov_vary %*% ti_basis

    exp(qnorm(p = ps, mean = log_haz, sd = sqrt(log_haz_var)))
  }))
  matplot(tis, hazs, lty = 1, type = "l", col = "black", bty = "l",
          xlab = xlab, xaxs = "i", yaxs = "i", ylab = ylab,
          ylim = range(hazs, 0), ...)
  grid()

  invisible(list(time = tis, hazard = hazs))
}

extend_newdata_n_eval <- function(term, x, newdata, ...){
  is_valid_expansion(term)
  if(!is.null(newdata)){
    stopifnot(is.data.frame(newdata), nrow(newdata) == 1)
    newdata <- lapply(newdata, replicate, n = length(x))
    newdata <- as.data.frame(newdata)
  }
  term$eval(x = x, newdata = newdata, ...)
}
