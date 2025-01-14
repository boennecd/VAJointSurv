#' Creates Data for One Type of Survival Outcome
#'
#' @inheritParams marker_term
#' @param formula a two-sided \code{\link{formula}} with the survival outcome
#' on the left-hand side and fixed effect covariates on the right-hand
#' side. The left-hand side needs to be a \code{\link[survival]{Surv}} object and can
#' be either right-censored and left-truncated.
#' @param time_fixef the time-varying fixed effects. See .e.g.
#' \code{\link{poly_term}}. This is for the baseline hazard. Note that many
#' basis expansions have boundary knots. It is important that these are set
#' to cover the full range of survival times including time zero for some
#' expansions.
#' @param data \code{\link{data.frame}} with at least the time variable.
#' @param with_frailty \code{TRUE} if there should be a frailty term.
#' @param delayed a vector with an entry which is \code{TRUE} if the
#' left-truncation time from the survival outcome is from a delayed entry.
#'
#' @details
#' The \code{time_fixef} should likely not include an intercept as this is
#' often included in \code{formula}.
#'
#' The \code{delayed} argument is to account for delayed entry with terminal
#' events when observations are sampled in a way such that they must not have
#' had the event prior to their left-truncation time. In this case, the proper
#' complete data likelihood is
#'
#' \deqn{\frac{a(u)h(t_{ij}\mid u)^{d_{ij}}S(t_{ij} \mid u)g(u)}{\int a(u) S(v_{ij} \mid u) du}}
#'
#' and not
#'
#' \deqn{a(u)h(t_{ij} \mid u)^{d_{ij}}\frac{S(t_{ij} \mid u)}{S(v_{ij} \mid u)}g(u)}
#'
#' where \eqn{h} is conditional hazard, \eqn{S} is the conditional survival
#' function, \eqn{g} is additional conditional likelihood factors from other
#' outcomes, \eqn{a} is the random effect distribution, \eqn{t_{ij}} is the
#' observed time, \eqn{d_{ij}} is an event indicator, and \eqn{v_{ij}} is the
#' left truncation time.
#'
#' The denominator in the proper complete likelihood becomes the expectation
#' over all delayed entries when a cluster has more than one delayed entry. See
#' van den Berg and Drepper (2016) and Crowther et al. (2016) for further
#' details.
#'
#' @references
#' Crowther MJ, Andersson TM, Lambert PC, Abrams KR & Humphreys K (2016).
#' \emph{Joint
#' modelling of longitudinal and survival data: incorporating delayed entry and
#' an assessment of model misspecification}. Stat Med,
#' 35(7):1193-1209. doi:10.1002/sim.6779
#'
#' van den Berg GJ & Drepper B (2016). \emph{Inference for Shared-Frailty
#' Survival Models with Left-Truncated Data}. Econometric Reviews, 35:6,
#' 1075-1098, doi: 10.1080/07474938.2014.975640
#'
#' @importFrom stats model.frame model.matrix model.response
#'
#' @return
#' An object of class \code{surv_term} with data required for survival outcome.
#'
#' @examples
#' # load in the data
#' library(survival)
#' data(pbc, package = "survival")
#'
#' # re-scale by year
#' pbcseq <- transform(pbcseq, day_use = day / 365.25)
#' pbc <- transform(pbc, time_use = time / 365.25)
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
#' @export
surv_term <- function(formula, id, data, time_fixef,
                      with_frailty = FALSE, delayed = NULL){
  # get the input data
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  id <- eval(substitute(id), data, parent.frame())
  Z <- model.matrix(mt, mf)
  y <- model.response(mf)

  stopifnot(inherits(y, "Surv"),
            attr(y, "type") %in% c("counting", "right"))
  if(attr(y, "type") == "right")
    y <- cbind(0, y)
  stopifnot(all(y[, 3] %in% 0:1),
            any(y[, 3] == 0),
            any(y[, 3] == 1))

  delayed <- eval(substitute(delayed), data, parent.frame())
  if(is.null(delayed))
    delayed <- logical(NROW(y))
  else
    delayed <- as.logical(delayed & y[, 1] > 0)

  time_fixef <- eval(substitute(time_fixef), data, parent.frame())

  # sanity checks
  stopifnot(NROW(Z) == length(id),
            is.matrix(Z),
            all(is.finite(Z)),
            all(is.finite(y)),
            NROW(Z) == NROW(y),
            is.logical(with_frailty), length(with_frailty) == 1,
            is.finite(with_frailty),
            length(delayed) == length(id),
            all(is.finite(delayed)))
  is_valid_expansion(time_fixef)

  # check for a singular design matrix
  XZ <- cbind(Z, t(time_fixef$eval(y[, 2],newdata = data)))
  rk <- rankMatrix(XZ)
  if(rk < NCOL(XZ))
    stop("Design matrix does not have full rank. Perhaps remove an intercept or a time-varying term from 'formula'")

  # reorder the data and return
  ord <- order(id, y[, 2])
  y <- y[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  id <- id[ord]
  delayed <- delayed[ord]
  data <- data[ord, ]

  fixef_design_varying <-
    bases_weights(time_fixef$weights_symbol,data,parent.frame(),nrow(y))

  structure(list(y = y, Z = t(Z), time_fixef = time_fixef, id = id, mt = mt, with_frailty = with_frailty, delayed = delayed,
                 fixef_design_varying = fixef_design_varying, data = data),
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
  comp <- ph_ll(
    time_fixef = object$time_fixef, Z = Z, surv = surv,
    fixef_design_varying = object$fixef_design_varying,
    rng_design_varying = object$rng_design_varying,
    object$with_frailty)

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
