#' Creates Data for One Type of Survival Outcome
#'
#' @inheritParams marker_term
#' @param formula a two-sided \code{\link{formula}} with the survival outcome
#' on the left-hand side and fixed effect covariates on the right-hand
#' side. The left-hand side needs to be a \code{\link{Surv}} object and can
#' be either right-censored and left-truncated.
#' @param time_fixef the time-varying fixed effects. See .e.g.
#' \code{\link{poly_term}}. This is for the baseline hazard. Note that many
#' basis expansions has boundary knots. It is important that these are set
#' to cover the full range of survival times including time zero for some
#' expansions.
#' @param data \code{\link{data.frame}} with at least the time variable.
#' @param ders a \code{\link{list}} with \code{\link{integer}} vectors for how
#' the survival outcome is linked to the markers. 0 implies present values,
#' -1 is integral of, and 1 is the derivative. \code{NULL} implies the present
#' value of the random effect for all markers.
#'
#' @details
#' The \code{time_fixef} should likely not include an intercept as this is
#' often included in \code{formula}.
#'
#' @importFrom stats model.frame model.matrix model.response
#'
#' @export
surv_term <- function(formula, id, data, time_fixef, ders = NULL){
  # get the input data
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  id <- eval(substitute(id), data)
  Z <- model.matrix(mt, mf)
  y <- model.response(mf)

  stopifnot(inherits(y, "Surv"),
            attr(y, "type") %in% c("counting", "right"))
  if(attr(y, "type") == "right")
    y <- cbind(0, y)
  stopifnot(all(y[, 3] %in% 0:1),
            any(y[, 3] == 0),
            any(y[, 3] == 1))

  time_fixef <- eval(substitute(time_fixef), data[y[, 3] == 1, ])

  # sanity checks
  stopifnot(NROW(Z) == length(id),
            is.matrix(Z),
            all(is.finite(Z)),
            all(is.finite(y)),
            NROW(Z) == NROW(y))
  is_valid_expansion(time_fixef)

  # check for a singular design matrix
  XZ <- cbind(Z, t(time_fixef$eval(y[, 2])))
  rk <- rankMatrix(XZ)
  if(rk < NCOL(XZ))
    stop("Design matrix does not have full rank. Perhaps remove an intercept or a time-varying term from 'formula'")

  # reorder the data and return
  ord <- order(id, y[, 2])
  y <- y[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  id <- id[ord]

  structure(list(y = y, Z = t(Z), time_fixef = time_fixef, id = id, mt = mt,
                 ders = ders),
            class = "surv_term")
}

# computes the starting values for the fixed effects
#' @importFrom stats optim
#' @importFrom utils head tail
surv_term_start_value <- function(object, quad_rule, va_var){
  stopifnot(inherits(object, "surv_term"))

  # create the object to perform the optimization
  Z <- object$Z
  surv <- t(object$y)
  comp <- ph_ll(time_fixef = object$time_fixef, Z = Z, surv = surv)

  # define the likelihood and gradient function. Then optimize
  ll <- function(par)
    ph_eval(ptr = comp$ptr, par = par, quad_rule = quad_rule, va_var = va_var)
  gr <- function(par)
    ph_grad(ptr = comp$ptr, par = par, quad_rule = quad_rule, va_var = va_var)

  # TODO: find a better starting value
  opt <- optim(numeric(comp$n_params), ll, gr, method = "BFGS")
  if(opt$convergence != 0)
    warning(sprintf(
      "optim returned convergence code %d when finding the starting values for the survival parameters",
      opt$convergence))

  list(fixef = head(opt$par, NROW(Z)), fixef_vary = tail(opt$par, -NROW(Z)))
}
