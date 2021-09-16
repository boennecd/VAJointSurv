#' Creates Data for One Type of Marker
#'
#' @param formula a two-sided \code{\link{formula}} with the marker outcome
#' on the left-hand side and fixed effect covariates on the right-hand
#' side.
#' @param id the variable for the id of each individual.
#' @param data a \code{data.frame} or environment to look at up the
#' variables in.
#' @param time_fixef the time-varying fixed effects. See .e.g.
#' \code{\link{poly_term}}.
#' @param time_rng the time-varying random effects. See .e.g.
#' \code{\link{poly_term}}.
#'
#' @details
#' The \code{time_fixef} should likely not include an intercept as this is
#' often included in \code{formula}. Use
#' \code{poly_term(degree = 0, raw = TRUE, intercept = TRUE)} if you want only
#' a random intercept.
#'
#' @importFrom stats model.frame model.matrix model.response
#'
#' @export
marker_term <- function(formula, id, data, time_fixef, time_rng){
  # get the input data
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  id <- eval(substitute(id), data)
  X <- model.matrix(mt, mf)
  y <- model.response(mf, "numeric")

  time_fixef <- eval(substitute(time_fixef), data)
  time_rng <- eval(substitute(time_rng), data)

  # sanity checks
  stopifnot(NROW(X) == length(id),
            !is.matrix(y),
            is.matrix(X),
            all(is.finite(X)),
            NROW(X) == length(y),
            NROW(X) == length(time_fixef$time),
            all(is.finite(time_fixef$time)),
            NROW(X) == length(time_rng$time),
            all(is.finite(time_rng$time)),
            all(time_rng$time == time_fixef$time))
  is_valid_expansion(time_fixef)
  is_valid_expansion(time_rng)

  # prepare the data to return
  time_var <-  time_fixef$time
  time_fixef$time <- NULL
  time_rng$time <- NULL

  # we have to sort by id and type
  ord <- order(id, time_var)
  id <- id[ord]
  time_var <- time_var[ord]
  X <- t(X)[, ord, drop = FALSE]
  y <- y[ord]

  structure(list(time = time_var, X = X, y = y, id = id, mt = mt,
                 time_fixef = time_fixef, time_rng = time_rng),
            class = "marker_term")
}

# computes the starting values for the fixed effect coefficients
#' @importFrom stats lm.fit
#' @importFrom utils head tail
marker_term_start_value <- function(object){
  stopifnot(inherits(object, "marker_term"))

  X <- t(object$X)
  Z <- t(object$time_fixef$eval(object$time))
  n_X <- NCOL(X)
  XZ <- cbind(X, Z)
  rm(X, Z)
  gc()
  y <- object$y
  fit <- lm.fit(XZ, y)
  coefs <- fit$coefficients

  list(fixef = head(coefs, n_X), fixef_vary = tail(coefs, -n_X),
       var = mean(fit$residuals^2))
}
