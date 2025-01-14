# use Gauss-Legendre quadrature by default
#' @importFrom SimSurvNMarker get_gl_rule
default_quad_rule <- function(){
  rule <- get_gl_rule(25L)
  within(rule, {
    node <- node / 2 + .5
    weight <- weight / 2
  })
}

set_n_check_quad_rule <- function(quad_rule){
  if(is.null(quad_rule))
    quad_rule <- default_quad_rule()

  with(quad_rule,
       stopifnot(is.list(quad_rule),
                 is.numeric(node), all(is.finite(node)),
                 all(node >= 0), all(node <= 1),
                 is.numeric(weight), all(is.finite(weight)),
                 length(node) == length(weight)))

  quad_rule
}

set_n_check_gh_quad_rule <- function(quad_rule){
  if(is.null(quad_rule))
    quad_rule <- list(node = c(-1.65068012388578, -0.52464762327529, 0.52464762327529, 1.65068012388578),
                      weight = c(0.0813128354472451, 0.804914090005513, 0.804914090005512, 0.0813128354472453))

  with(quad_rule,
       stopifnot(is.list(quad_rule),
                 is.numeric(node), all(is.finite(node)),
                 is.numeric(weight), all(is.finite(weight)),
                 length(node) == length(weight)))

  quad_rule
}

check_n_threads <- function(object, n_threads){
  stopifnot(length(n_threads) == 1, n_threads > 0)
  if(object$max_threads < n_threads)
    warning("Maximum number of threads is less then n_threads")
}

#' Creates a joint_ms Object to Estimate a Joint Survival and Marker Model
#'
#' @param markers either an object from \code{\link{marker_term}} or a list
#' of such objects.
#' @param survival_terms either an object from \code{\link{surv_term}} or a list
#' of such objects.
#' @param max_threads maximum number of threads to use.
#' @param quad_rule list with nodes and weights for a quadrature rule for the
#' integral from zero to one.
#' @param cache_expansions \code{TRUE} if the expansions in the numerical
#' integration in the survival parts of the lower bound should be cached (not
#' recomputed). This requires more memory and may be an advantage
#' particularly with
#' expansions that take longer to compute (like \code{\link{ns_term}} and
#' \code{\link{bs_term}}). The computation time may be worse particularly if
#' you use more threads as the CPU cache is not well utilized.
#' @param gh_quad_rule list with two numeric vectors called node and weight
#' with Gauss–Hermite quadrature nodes and weights to handle delayed entry.
#' A low number of quadrature nodes and weights is used when \code{NULL} is
#' passed.
#' This seems to work well when delayed entry happens at time with large
#' marginal survival probabilities. The nodes and weights can be obtained e.g.
#' from \code{fastGHQuad::gaussHermiteData}.
#'
#' @param ders a \code{\link{list}} of \code{\link{list}}s with
#' \code{\link{integer}} vectors for how
#' the survival outcomes are linked to the markers. 0 implies present values,
#' -1 is integral of, and 1 is the derivative. \code{NULL} implies the present
#' value of the random effect for all markers. Note that the number of integer
#' vectors should be equal to the number of markers.
#'
#' @return
#' An object of \code{joint_ms} class with the needed C++ and R objects
#' to estimate the model.
#' @seealso \code{\link{joint_ms_opt}}, \code{\link{joint_ms_lb}},
#' \code{\link{joint_ms_hess}}, and \code{\link{joint_ms_start_val}}.
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L)
#' @export
joint_ms_ptr <- function(markers = list(), survival_terms = list(),
                         max_threads = 1L, quad_rule = NULL,
                         cache_expansions = TRUE, gh_quad_rule = NULL,
                         ders = NULL){
  stopifnot(
    length(max_threads) == 1, max_threads > 0,
    is.logical(cache_expansions), length(cache_expansions) == 1)

  # handle defaults
  if(inherits(markers, "marker_term"))
    markers <- list(markers)
  else
    stopifnot(all(sapply(markers, inherits, "marker_term")))
  if(inherits(survival_terms, "surv_term"))
    survival_terms <- list(survival_terms)
  else
    stopifnot(all(sapply(survival_terms, inherits, "surv_term")))

  quad_rule <- set_n_check_quad_rule(quad_rule)
  gh_quad_rule <- set_n_check_gh_quad_rule(gh_quad_rule)

  if(is.null(ders)) ders <- replicate(length(survival_terms),
                                      replicate(length(markers), 0L,
                                                simplify = FALSE),
                                      simplify = FALSE)

  if(length(survival_terms) == 1 && !is.list(ders[[1]])) ders <- list(ders)

  stopifnot(is.list(ders), length(ders)==length(survival_terms))

  for(der in ders) stopifnot(is.list(der), length(der) == length(markers),
                             all(sapply(der,is.integer)))

  # adding ders
  survival_terms <- Map(function(x, ders){
    x$ders <- ders
    x
  }, survival_terms, ders)

  # we need to to alter the start times and add new observations for the delayed
  # entries
  marker_symbols <- unlist(lapply(markers,function(marker) {
    marker$time_rng$weights_symbol
  }))
  bases_weights_enclose <- parent.frame()
  survival_terms <- lapply(survival_terms, function(x){
    x$rng_design_varying <- bases_weights(
      marker_symbols, x$data, bases_weights_enclose, length(x$delayed))
    x
  })

  survival_terms_org <- survival_terms
  survival_terms <- lapply(survival_terms, function(x){
    idx_delayed <- which(x$delayed)

    if(length(idx_delayed) > 0){
      y_delayed <- x$y[idx_delayed, "start"]
      x$y[idx_delayed, "start"] <- 0
    }
    else
      y_delayed <- numeric()

    id_delayed <- x$id[idx_delayed]
    Z_delayed <- x$Z[, idx_delayed, drop = FALSE]
    fixef_design_varying_delayed <- x$fixef_design_varying[
      , idx_delayed, drop = FALSE]
    rng_design_varying_delayed <- x$rng_design_varying[
      , idx_delayed, drop = FALSE]

    x$delayed_data <- list(
      y = y_delayed, id = id_delayed, Z = Z_delayed,
      fixef_design_varying = fixef_design_varying_delayed,
      rng_design_varying = rng_design_varying_delayed)
    x
  })

  # extract all the delayed entries
  delayed_terms <- lapply(survival_terms, `[[`, "delayed_data")

  # create the C++ object
  ptr <- .joint_ms_ptr(
    markers, survival_terms, max_threads = max_threads,
    delayed_terms = delayed_terms)
  param_names <- joint_ms_parameter_names(ptr)
  out <- list(param_names = param_names, ptr = ptr)
  indices <- joint_ms_parameter_indices(ptr)

  # create a vector with the unique ids
  ids <- c(lapply(markers, `[[`, "id"),
           lapply(survival_terms, `[[`, "id"))
  ids <- sort(unique(unlist(ids)))

  # compute starting values. Start with the markers
  start_val <- numeric(joint_ms_n_params(ptr))
  m_start_val <- lapply(markers, marker_term_start_value)
  sigma <- matrix(0, length(markers), length(markers))
  for(i in seq_along(markers)){
    start_val[indices$markers[[i]]$fixef] <- m_start_val[[i]]$fixef
    start_val[indices$markers[[i]]$fixef_vary] <- m_start_val[[i]]$fixef_vary
    sigma[i, i] <- m_start_val[[i]]$var
  }

  if(length(markers) > 0)
    start_val[indices$vcovs$vcov_marker] <- .log_chol(sigma)

  # compute the starting values for the survival outcomes. Set the covariance
  # matrix for the frailties to some low value
  def_frailty_var <- 1e-2
  s_start_val <- lapply(
    # important that we pass the original survival_terms here
    survival_terms_org, surv_term_start_value,
    quad_rule = quad_rule, va_var = def_frailty_var)
  for(i in seq_along(survival_terms)){
    start_val[indices$survival[[i]]$fixef] <- s_start_val[[i]]$fixef
    start_val[indices$survival[[i]]$fixef_vary] <- s_start_val[[i]]$fixef_vary
  }

  out <- structure(
    list(param_names = param_names, indices = indices, ptr = ptr,
         start_val = start_val, max_threads = max_threads,
         quad_rule = quad_rule, gh_quad_rule = gh_quad_rule,
         n_lb_terms = joint_ms_n_terms(ptr),
         cache_expansions = cache_expansions, ids = ids,
         markers = markers, survival_terms = survival_terms,
         delayed_terms = delayed_terms,
         survival_terms_org = survival_terms_org),
    class = "joint_ms")

  # set the remaining covariance parameters
  n_frailties <-
    if(length(survival_terms) > 0)
      sum(sapply(survival_terms, `[[`, "with_frailty")) else 0L
  start_val <- joint_ms_set_vcov(out,
                                 vcov_vary = diag(indices$va_dim - n_frailties),
                                 vcov_surv = diag(def_frailty_var, n_frailties))
  out$start_val <- start_val
  out
}


#' Sets the Covariance Parameters
#'
#' @description
#' Sets the covariance matrices to the passed values. The function also sets
#' covariance matrices for the variational distributions to the same values.
#'
#' @inheritParams joint_ms_format
#' @param vcov_vary the covariance matrix for the time-varying effects.
#' @param vcov_surv the covariance matrix for the frailties.
#' @param va_mean a matrix with the number of rows equal to the number of
#' random effects per observation and the number of columns is the number
#' of observations. The order for the observations needs to be the same as the
#' \code{id} element of \code{object}.
#'
#' @return
#' Numeric vector with model parameters.
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#' # compute var-covar matrices with the first set of starting values
#' joint_ms_format(object = model_ptr)$vcov
#' joint_ms_va_par(object = model_ptr)[[1]]
#'
#' # altering var-covar matrices
#' alter_pars <- joint_ms_set_vcov(
#'   object = model_ptr,
#'   vcov_vary = diag(1:4),
#'   vcov_surv = matrix(0,0,0))
#'
#' # altered var-covar matrices
#' joint_ms_format(object = model_ptr, par = alter_pars)$vcov
#' joint_ms_va_par(object = model_ptr, par = alter_pars)[[1]]
#' @export
joint_ms_set_vcov <- function(
  object, vcov_vary, vcov_surv, par = object$start_val, va_mean = NULL){
  stopifnot(inherits(object, "joint_ms"))
  stopifnot(is.null(va_mean) || (is.matrix(va_mean) && all(is.finite(va_mean))))

  # insert the model parameter
  indices <- object$indices
  if(length(indices$vcovs$vcov_vary) > 0)
    par[indices$vcovs$vcov_vary] <- .log_chol(vcov_vary)
  if(length(indices$vcovs$vcov_surv) > 0)
    par[indices$vcovs$vcov_surv] <- .log_chol(vcov_surv)

  # fill in default values for the VA parameters
  va_dim <- indices$va_dim
  va_vcov <- matrix(0, va_dim, va_dim)
  n_vary <- NCOL(vcov_vary)
  if(n_vary > 0){
    va_vcov[1:n_vary, 1:n_vary] <- vcov_vary
    if(NCOL(vcov_surv) > 0)
      va_vcov[-(1:n_vary), -(1:n_vary)] <- vcov_surv
  } else
    va_vcov <- vcov_surv

  if(is.null(va_mean))
    va_mean <- matrix(0, va_dim, length(object$ids))
  else
    stopifnot(all(dim(va_mean) == c(va_dim, length(object$ids))))

  va_default <- tapply(va_mean, rep(1:NCOL(va_mean), each = NROW(va_mean)),
                       c, .log_chol(va_vcov))
  va_default <- unlist(va_default)
  par[-(1:(indices$va_params_start - 1L))] <- va_default
  par
}

#' Quick Heuristic for the Starting Values
#'
#' @inheritParams joint_ms_opt
#'
#' @return
#' Numeric vector of starting values for the model parameters.
#'
#' @importFrom stats coef
#' @importFrom psqn psqn
#' @importFrom lme4 lmer VarCorr lmerControl
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#' # find the starting values
#' start_vals <- joint_ms_start_val(model_ptr)
#'
#' @export
joint_ms_start_val <- function(
  object, par = object$start_val,
  rel_eps = 1e-8,  max_it = 1000L, n_threads = object$max_threads,
  c1 = 1e-4, c2 = .9, use_bfgs = TRUE, trace = 0, cg_tol = 0.5,
  strong_wolfe = TRUE, max_cg = 0, pre_method = 3L,
  quad_rule = object$quad_rule, mask = integer(),
  cache_expansions = object$cache_expansions, gr_tol = -1,
  gh_quad_rule = object$gh_quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  quad_rule <- set_n_check_quad_rule(quad_rule)
  gh_quad_rule <- set_n_check_gh_quad_rule(gh_quad_rule)
  check_n_threads(object, n_threads)
  stopifnot(is.integer(mask), all(mask >= 0 & mask < length(par)))

  # set the marker parameters using lme4
  vcov_marker <- .log_chol_inv(par[object$indices$vcovs$vcov_marker])
  vcov_marker <- diag(NCOL(vcov_marker))
  vcov_vary <- .log_chol_inv(par[object$indices$vcovs$vcov_vary])
  vcov_surv <- if(length(object$survival_terms) == 0)
    matrix(nrow = 0, ncol = 0) else
      .log_chol_inv(par[object$indices$vcovs$vcov_surv])

  n_rng <- NCOL(vcov_vary) + NCOL(vcov_surv)
  blups <- matrix(0, n_rng, length(object$ids))

  offset <- 0L
  for(i in seq_along(object$markers)){
    mark <- object$markers[[i]]
    n_fixef <- NROW(mark$X)
    X_vary <- mark$time_fixef$eval(mark$time, newdata = mark$data)
    n_fixef_vary <- NROW(X_vary)
    X_rng <- mark$time_rng$eval(mark$time, newdata = mark$data)
    y <- mark$y
    id <- mark$id
    lmer_fit <- lmer(
      y ~ cbind(t(mark$X), t(X_vary)) - 1 + (t(X_rng) - 1 | id),
      control = lmerControl(
        optimizer = "nloptwrap", optCtrl = list(
          xtol_abs = rel_eps, ftol_abs = rel_eps),
        check.nobs.vs.nRE = "warning"))

    # extract the parameters
    par[object$indices$markers[[i]]$fixef] <-
      lmer_fit@beta[seq_len(n_fixef)]
    par[object$indices$markers[[i]]$fixef_vary] <-
      lmer_fit@beta[seq_len(n_fixef_vary) + n_fixef]

    vcov_est <- VarCorr(lmer_fit)
    resid_var <- attr(vcov_est, "sc")^2
    vcov_est <- vcov_est[["id"]]

    vcov_marker[i, i] <- resid_var
    vcov_indices <- offset + seq_len(NROW(X_rng))
    vcov_vary[vcov_indices, vcov_indices] <- vcov_est

    blups_i <- coef(lmer_fit)[["id"]][, 1:NCOL(vcov_est)]
    id_match <- match(rownames(blups_i), object$ids)
    stopifnot(all(is.finite(id_match)))
    blups[vcov_indices, id_match] <- t(blups_i)

    offset <- offset + NROW(X_rng)
  }

  # set the parameters
  if(length(object$markers) > 0){
    par[object$indices$vcovs$vcov_marker] <- .log_chol(vcov_marker)
    par[object$indices$vcovs$vcov_vary] <- .log_chol(vcov_vary)

    par <- joint_ms_set_vcov(object, vcov_vary = vcov_vary,
                             vcov_surv = vcov_surv, par = par,
                             va_mean = blups)
  }

  # optimizes the VA parameters
  do_opt_priv <- function(x){
    res <- opt_priv(
      val = x, ptr = object$ptr, rel_eps = rel_eps,
      max_it = max_it, n_threads = n_threads, c1 = c1, c2 = c2,
      quad_rule = quad_rule, cache_expansions = cache_expansions,
      gr_tol = gr_tol, gh_quad_rule = gh_quad_rule)

    # this may lead to semi-definite covariance matrices. We ad-hock adjust
    # these which may help
    va_dim <- object$indices$va_dim
    n_va_params <- object$indices$n_va_params
    va_params_start <- object$indices$va_params_start
    n_ids <- (length(x) - va_params_start + 1L) / n_va_params

    new_va_par <- tapply(
      res[-(1:(va_params_start - 1L))],
      rep(1:n_ids, each = n_va_params), function(va_pars){
        va_cov <- .log_chol_inv(va_pars[-(1:va_dim)])
        eg <- eigen(va_cov)
        threshold <- 1e6
        if(max(eg$values) < threshold * min(eg$values))
          return(va_pars)

        eg$values <- pmax(eg$values, max(eg$values) / threshold)
        va_cov <- tcrossprod(eg$vectors %*% diag(eg$values), eg$vectors)
        return(c(va_pars[1:va_dim], .log_chol(va_cov)))
      })

    c(res[1:(va_params_start - 1L)], unlist(new_va_par))
  }

  # optimizes the marker parameters
  opt_marker_par <- function(x){
    res <- joint_ms_opt_lb(val = x, ptr = object$ptr, rel_eps = rel_eps,
                           max_it = max_it, n_threads = n_threads, c1 = c1,
                           c2 = c2, use_bfgs = use_bfgs, trace = trace,
                           cg_tol = cg_tol, strong_wolfe = strong_wolfe,
                           max_cg = max_cg, pre_method = pre_method,
                           quad_rule = quad_rule, mask = mask,
                           cache_expansions = cache_expansions,
                           only_markers = TRUE, gr_tol = gr_tol,
                           gh_quad_rule = gh_quad_rule)

    if(NCOL(vcov_surv))
      res$par <- do_opt_priv(res$par)

    structure(
      res$par,
      value = joint_ms_lb(
        object = object, par = res$par, n_threads = n_threads,
        quad_rule = quad_rule, cache_expansions = cache_expansions),
      conv = res$conv)
  }

  res <- opt_marker_par(par)
  if(is.finite(attr(res, "value")) && attr(res, "conv"))
    return(res)

  # this might help
  opt_marker_par(do_opt_priv(par))
}

#' Evaluates the Lower Bound or the Gradient of the Lower Bound
#'
#' @return
#' \code{joint_ms_lb} returns a number scalar with the lower bound.
#'
#' @inheritParams joint_ms_ptr
#' @param object a joint_ms object from \code{\link{joint_ms_ptr}}.
#' @param par parameter vector for where the lower bound is evaluated at.
#' @param n_threads number of threads to use. This is not supported on Windows.
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#' # find the starting values
#' start_vals <- joint_ms_start_val(model_ptr)
#'
#' # same lower bound
#' all.equal(attr(start_vals,"value"),joint_ms_lb(model_ptr,par = start_vals))
#' @export
joint_ms_lb <- function(object, par, n_threads = object$max_threads,
                        quad_rule = object$quad_rule,
                        cache_expansions = object$cache_expansions,
                        gh_quad_rule = object$gh_quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  quad_rule <- set_n_check_quad_rule(quad_rule)
  check_n_threads(object, n_threads)
  gh_quad_rule <- set_n_check_gh_quad_rule(gh_quad_rule)

  joint_ms_eval_lb(val = par, ptr = object$ptr, n_threads = n_threads,
                   quad_rule = quad_rule, cache_expansions = cache_expansions,
                   gh_quad_rule = gh_quad_rule)
}

#' @rdname joint_ms_lb
#'
#' @return
#' \code{joint_ms_lb_gr} returns a numeric vector with the gradient.
#'
#' @export
joint_ms_lb_gr <- function(object, par, n_threads = object$max_threads,
                           quad_rule = object$quad_rule,
                           cache_expansions = object$cache_expansions,
                           gh_quad_rule = object$gh_quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  quad_rule <- set_n_check_quad_rule(quad_rule)
  gh_quad_rule <- set_n_check_gh_quad_rule(gh_quad_rule)
  check_n_threads(object, n_threads)

  joint_ms_eval_lb_gr(val = par, ptr = object$ptr, n_threads = n_threads,
                      quad_rule = quad_rule, cache_expansions = cache_expansions,
                      gh_quad_rule = gh_quad_rule)
}

#' Computes the Hessian
#'
#' @inheritParams joint_ms_lb
#' @param eps,scale,tol,order parameter to pass to psqn. See
#' \code{\link[psqn]{psqn_hess}}.
#'
#' @return
#' A list with the following two Hessian matrices:
#' \item{\code{hessian}}{Hessian matrix of the model parameters with the
#' variational parameters profiled out.}
#' \item{\code{hessian_all}}{Hessian matrix of the model and variational parameters.}
#'
#' @import methods
#' @importFrom Matrix solve
#' @importMethodsFrom Matrix solve

#' @examples
#' \donttest{# load in the data
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#' # find the starting values
#' start_vals <- joint_ms_start_val(model_ptr)
#'
#' # optimize lower bound
#' fit <- joint_ms_opt(object = model_ptr, par = start_vals, gr_tol = .01)
#'
#' # compute the Hessian
#' hess <- joint_ms_hess(object = model_ptr,par = fit$par)
#'
#' # standard errors of the parameters
#' library(Matrix)
#' sqrt(diag(solve(hess$hessian))) }
#' @export
joint_ms_hess <- function(
  object, par, quad_rule = object$quad_rule,
  cache_expansions = object$cache_expansions, eps = 1e-4, scale = 2,
  tol = .Machine$double.eps^(3/5), order = 4L,
  gh_quad_rule = object$gh_quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  quad_rule <- set_n_check_quad_rule(quad_rule)
  gh_quad_rule <- set_n_check_gh_quad_rule(gh_quad_rule)

  res <- .joint_ms_hess(val = par, ptr = object$ptr, quad_rule =  quad_rule,
                        cache_expansions = cache_expansions, eps = eps,
                        scale = scale, tol = tol, order = order,
                        gh_quad_rule = gh_quad_rule)

  # try to compute
  is_global <- 1:(object$indices$va_params_start - 1L)
  hess_model_par <- try({
    v1 <- res[is_global, -is_global] %*% solve(
      res[-is_global, -is_global], res[-is_global, is_global])
    res[is_global, is_global] - v1
  }, silent = TRUE)

  if(inherits(hess_model_par, "try-error")){
    warning("Failed to compute the hessian of the model parametes")
    hess_model_par <- matrix(NA, length(is_global), length(is_global))
  }

  list(hessian = hess_model_par, hessian_all = res)
}

#' Optimizes the Lower Bound
#'
#' @inheritParams joint_ms_lb
#' @param par starting value.
#' @param rel_eps,max_it,c1,c2,use_bfgs,trace,cg_tol,strong_wolfe,max_cg,pre_method,mask,gr_tol
#' arguments to pass to the C++ version of \code{\link[psqn]{psqn}}.
#'
#' @return
#' A list with the following elements:
#' \item{\code{par}}{numeric vector of estimated model parameters.}
#' \item{\code{value}}{numeric scalar with the value of optimized lower bound.}
#' \item{\code{counts}}{integer vector with the function counts and the number
#' of conjugate gradient iterations. See \code{\link[psqn]{psqn}}.}
#' \item{\code{convergence}}{logical for whether the optimization converged.}
#'
#' @examples
#' \donttest{# load in the data
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#' # find the starting values
#' start_vals <- joint_ms_start_val(model_ptr)
#'
#' # optimize lower bound
#' fit <- joint_ms_opt(object = model_ptr, par = start_vals, gr_tol = .01)
#'
#' # formatted maximum likelihood estimators
#' joint_ms_format(model_ptr, fit$par)}
#' @export
joint_ms_opt <- function(
  object, par = object$start_val, rel_eps = 1e-8, max_it = 1000L,
  n_threads = object$max_threads, c1 = 1e-4, c2 = .9, use_bfgs = TRUE,
  trace = 0L, cg_tol = .5, strong_wolfe = TRUE, max_cg = 0L,
  pre_method = 3L, quad_rule = object$quad_rule, mask = integer(),
  cache_expansions = object$cache_expansions, gr_tol = -1,
  gh_quad_rule = object$gh_quad_rule){
  stopifnot(inherits(object, "joint_ms"))
  quad_rule <- set_n_check_quad_rule(quad_rule)
  gh_quad_rule <- set_n_check_gh_quad_rule(gh_quad_rule)
  check_n_threads(object, n_threads)
  stopifnot(is.integer(mask), all(mask >= 0 & mask < length(par)))

  fit <- joint_ms_opt_lb(val = par, ptr = object$ptr, rel_eps = rel_eps,
                         max_it = max_it, n_threads = n_threads, c1 = c1, c2 = c2,
                         use_bfgs = use_bfgs, trace = trace, cg_tol = cg_tol,
                         strong_wolfe = strong_wolfe, max_cg = max_cg,
                         pre_method = pre_method, quad_rule = quad_rule,
                         mask = mask,
                         cache_expansions = cache_expansions,
                         only_markers = FALSE, gr_tol = gr_tol,
                         gh_quad_rule = gh_quad_rule)
  if(!fit$convergence)
    warning(sprintf("Fit did not converge but returned with code %d. Perhaps increase the maximum number of iterations",
                    fit$info))
  fit
}

#' Formats the Parameter Vector
#'
#' @description
#' Formats a parameter vector by putting the model parameters into a \code{list}
#' with elements for each type of parameter.
#'
#' @param object a joint_ms object from \code{\link{joint_ms_ptr}}.
#' @param par parameter vector to be formatted.
#'
#' @return
#' A list with the following elements:
#' \item{markers}{list with an element for each marker. The lists contains an
#' element called \code{fixef} for non-time-varying fixed effects and an
#' element called \code{fixef_vary} time-varying fixed effects.}
#' \item{survival}{list with an element for each survival outcome. The lists
#' contains an element called \code{fixef} for non-time-varying fixed effects,
#' an element called \code{fixef_vary} time-varying fixed effects, and an
#' element called \code{associations} for the association parameters.}
#' \item{vcov}{contains three covariance matrices called \code{vcov_marker},
#' \code{vcov_vary} and \code{vcov_surv} for the covariance matrix of the
#' markers error term, the time-varying random effects, and the frailties,
#' respectively.}
#'
#' @importFrom stats setNames
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#' # find the starting values
#' start_vals <- joint_ms_start_val(model_ptr)
#'
#' # format the starting values
#' joint_ms_format(model_ptr,start_vals)
#' @export
joint_ms_format <- function(object, par = object$start_val){
  # TODO: add tests for this function
  stopifnot(inherits(object, "joint_ms"))
  indices <- object$indices

  # handle the markers
  markers <- lapply(seq_along(indices$markers), function(i)
    list(fixef = par[indices$markers[[i]]$fixef],
         fixef_vary = par[indices$markers[[i]]$fixef_vary]))

  # handle the survival outcomes
  survival <- lapply(seq_along(indices$survival), function(i)
    list(fixef = par[indices$survival[[i]]$fixef],
         fixef_vary = par[indices$survival[[i]]$fixef_vary],
         associations = par[indices$survival[[i]]$associations]))

  # the covariance matrices
  vcov <-
    with(indices$vcovs, list(
      vcov_marker = .log_chol_inv(par[vcov_marker]),
      vcov_surv = .log_chol_inv(par[vcov_surv]),
      vcov_vary = .log_chol_inv(par[vcov_vary])))

  out <- list(markers = markers, survival = survival, vcov = vcov)

  # remove objects with length zero
  skip_if_length_zero <- function(x){
    if(length(x) == 0)
      NULL
    else if(is.list(x)){
      x <- lapply(x, skip_if_length_zero)
      lens <- lengths(x)
      x[lengths(x) > 0]
    } else
      x
  }
  skip_if_length_zero(out)
}

#' Approximate Likelihood Ratio based Confidence Intervals
#'
#' @inheritParams joint_ms_opt
#' @param opt_out maximum lower bound estimator from \code{\link{joint_ms_opt}}.
#' @param which_prof index of the parameter to profile.
#' @param delta numeric scalar greater than zero for the initial step size.
#' Steps are made of size \code{2^(iteration - 1) * delta}. A guess of the
#' standard deviation is a good value.
#' @param verbose logical for whether to print output during the construction
#' of the approximate profile likelihood curve.
#' @param max_step maximum number of steps to take in each direction when
#' constructing the approximate profile likelihood curve.
#' @param level confidence level.
#' @param hess the Hessian from \code{\link{joint_ms_hess}}. It is used to get
#' better starting values along the profile likelihood curve. Use \code{NULL}
#' if it is not passed.
#'
#' @return
#' A list with the following elements:
#' \item{confs}{profile likelihood based confidence interval.}
#' \item{xs}{the value of the parameter at which the profile likelihood is evaluated at.}
#' \item{p_log_Lik}{numeric scalar with the profile log-likelihood.}
#' \item{data}{list of lists of the output of each point where the profile
#' likelihood is evaluated with the optimal parameter values of the
#' other parameters given the
#' constrained value of the parameter that is being profiled and
#' the optimal value of the lower bound.}
#'
#' @importFrom stats approx qchisq splinefun qnorm spline
#' @importFrom utils head
#' @examples
#' \donttest{ # load in the data
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
#' # create the C++ object to do the fitting
#' model_ptr <- joint_ms_ptr(
#'   markers = list(m1, m2), survival_terms = s_term,
#'   max_threads = 2L, ders = list(0L, c(0L, -1L)))
#'
#'
#' # find the starting values
#' start_vals <- joint_ms_start_val(model_ptr)
#'
#' # optimize lower bound
#' fit <- joint_ms_opt(object = model_ptr, par = start_vals, gr_tol = .01)
#'
#' # compute the Hessian
#' hess <- joint_ms_hess(object = model_ptr,par = fit$par)
#'
#' # compute the standard errors
#' library(Matrix)
#' se <- sqrt(diag(solve(hess$hessian)))
#'
#' # find index for the first association parameter
#' which_prof <- model_ptr$indices$survival[[1]]$associations[1]
#'
#' # initial step size for finding the confidence interval limits
#' delta <- 2*se[which_prof]
#'
#' # compute profile likelihood based confidence interval
#' # for the first association parameter
#' profile_CI <- joint_ms_profile(
#'   object = model_ptr, opt_out = fit, which_prof = which_prof,
#'   delta= delta, gr_tol = .01)
#'
#' # comparison of CIs
#' profile_CI$confs
#' fit$par[which_prof]+c(-1,1)*qnorm(0.975)*se[which_prof] }
#' @export
joint_ms_profile <- function(
  object, opt_out, which_prof, delta, level = .95, max_step = 15L,
  rel_eps = 1e-8, max_it = 1000L, n_threads = object$max_threads, c1 = 1e-04,
  c2 = 0.9, use_bfgs = TRUE, trace = 0L, cg_tol = 0.5, strong_wolfe = TRUE,
  max_cg = 0L, pre_method = 3L, quad_rule = object$quad_rule, verbose = TRUE,
  mask = integer(), cache_expansions = object$cache_expansions,
  gr_tol = -1, hess = NULL){
  stopifnot(is.integer(which_prof), length(which_prof) == 1L,
            which_prof >= 1L, which_prof <= length(opt_out$par),
            is.numeric(delta), is.finite(delta), length(delta) == 1,
            delta > 0,
            length(level) == 1, is.numeric(level), is.finite(level),
            level > 0, level < 1,
            opt_out$convergence,
            is.null(hess) || (is.list(hess) && !is.null(hess$hessian_all)))

  # setup
  mask <- unique(c(mask, which_prof - 1L))
  n_global <- object$indices$va_params_start - 1L
  par <- opt_out$par
  max_lb <- -opt_out$value

  chi_val <- qchisq(level, 1)
  crit_value <- max_lb - chi_val / 2

  # formats the output from one estimation
  format_fit <- function(x, fit, dir, ub = NULL, lb = NULL){
    if(verbose){
      new_val <- -fit$value
      message(sprintf("LogLike: %.4f at %16f", new_val, x), appendLF = FALSE)
      if(!is.null(ub) && !is.null(lb)){
        if(new_val < crit_value)
          lb$value <- new_val
        else
          ub$value <- new_val
        message(sprintf(". Lb, target, ub: %.4f, %.4f, %.4f",
                        lb$value, crit_value, ub$value))

      } else
        message() # line feed
    }
    list(x = x, value = -fit$value, par = fit$par,
         z_val = sign(dir) * sqrt((max_lb + fit$value) * 2))
  }

  # fits the model
  fit_model <- function(x, dir, par, lb = NULL, ub = NULL){
    par[which_prof] <- x
    fit <- joint_ms_opt(
      object = object, par = par, rel_eps = rel_eps, max_it = max_it,
      n_threads = n_threads, c1 = c1, c2 = c2, use_bfgs = use_bfgs,
      trace = trace, cg_tol = cg_tol, strong_wolfe = strong_wolfe,
      max_cg = max_cg, pre_method = pre_method, quad_rule = quad_rule,
      mask = mask, cache_expansions = cache_expansions, gr_tol = gr_tol)
    if(!fit$convergence)
      stop(sprintf("Fit did not converge but returned with code %d. Perhaps increase the maximum number of iterations",
                   fit$info))

    format_fit(x, dir = dir, fit = fit, lb = lb, ub = ub)
  }

  # objects needed for to get the starting value when the Hessian is passed.
  # Relies on an approximate normal distribution where we work with the
  # precision matrix rather than the covariance matrix
  if(!is.null(hess)){
    if(verbose)
      message("Computing the vector to find the starting values along the profile likelihood curve")

    cond_mean_factor <- -solve(
      hess$hessian_all[-which_prof, -which_prof],
      hess$hessian_all[-which_prof, which_prof])
    cond_mean_factor <- drop(as.matrix(cond_mean_factor))
  }

  # gets the starting value
  get_start_val <- function(step_size, closest_par){
    if(is.null(hess))
      return(closest_par)

    # use approximate normal distribution
    out <- par
    out[-which_prof] <- out[-which_prof] + cond_mean_factor * step_size
    out
  }

  # find the points along the approximate profile likelihood curve
  get_points <- function(dir){
    dir <- sign(dir)
    if(verbose)
      message(sprintf(
        "\nFinding the %s limit of the approximate profile likelihood curve",
        if(dir < 0) "lower" else "upper"))

    step <- 0L
    out <- vector("list", max_step)
    prev <- max_lb
    did_fail <- FALSE

    closest_par <- numeric()
    closests_diff <- Inf

    update_closest_par <- function(new_fit){
      new_diff <- abs(new_fit$value - crit_value)
      if(new_diff < closests_diff){
        closest_par <<- new_fit$par
        closests_diff <<- new_diff
      }
    }
    tmp <- opt_out
    tmp$value <- -tmp$value
    update_closest_par(tmp)

    while(prev > crit_value && (step <- step + 1L) <= max_step){
      step_size <- dir * 2^(step - 1) * delta
      out[[step]] <- fit_model(
        par[which_prof] + step_size, dir = dir,
        par = get_start_val(step_size = step_size, closest_par = closest_par))
      if(out[[step]]$value > prev){
        warning("The lower bound did not decrease. Either the input is not an optimum or the precision needs to be increased")
        did_fail <- TRUE
        break
      }

      prev <- out[[step]]$value
      update_closest_par(out[[step]])
    }

    .report_failed <- function()
      if(verbose && step > max_step)
        warning(sprintf(
          "Failed to find the appropiate point in %d steps", max_step))

    .report_failed()

    if(did_fail || step > max_step)
      return(out[sapply(out, length) > 0])

    ub <- if(step == 1L)
      format_fit(par[which_prof], opt_out, 1) else out[[step - 1L]]
    lb <- out[[step]]

    while(ub$value - lb$value > chi_val / 6 && (step <- step + 1L) <= max_step){
      # compute the next value
      xs <- c(unlist(sapply(out, `[[`, "x"))    , par[which_prof])
      ys <- c(unlist(sapply(out, `[[`, "value")), max_lb)
      sp <- splinefun(ys, xs, method = "monoH.FC")
      y <- seq(lb$value, ub$value, length.out = 20)
      x <- sp(y)
      next_val <- approx(y, x, xout = crit_value)$y
      if(abs(next_val - ub$x) > abs(next_val - lb$x))
        next_val <- 8 / 9 * next_val + ub$x / 9
      else
        next_val <- 8 / 9 * next_val + lb$x / 9

      step_size <- next_val - par[which_prof]
      out[[step]] <- fit_model(
        next_val, dir = dir, lb = lb, ub = ub,
        get_start_val(step_size = step_size, closest_par = closest_par))

      if(out[[step]]$value > ub$value || out[[step]]$value < lb$value){
        warning("Log likelihood does not seem monotonic. Likely the precision needs to be increased")
        break
      }

      update_closest_par(out[[step]])

      if(out[[step]]$value < crit_value)
        lb <- out[[step]]
      else
        ub <- out[[step]]
    }

    .report_failed()
    out <- out[sapply(out, length) > 0]
    out[order(sapply(out, `[[`, "x"))]
  }

  res_up   <- get_points( 1)
  # get rid of the extra parameters to save space
  res_up <- lapply(res_up, function(x){
    x$par <- head(x$par, n_global)
    x
  })
  res_down <- get_points(-1)
  # get rid of the extra parameters to save space
  res_down <- lapply(res_down, function(x){
    x$par <- head(x$par, n_global)
    x
  })
  out <- c(
    res_down, list(format_fit(par[which_prof], opt_out, dir = 0)), res_up)

  # compute the confidence interval
  xs  <- sapply(out, `[[`, "x")
  zs  <- sapply(out, `[[`, "z_val")
  pls <- sapply(out, `[[`, "value")
  sp <- spline(xs, zs)
  alpha <- 1 - level
  pvs <- c(alpha / 2, 1 - alpha/2)
  confs <- setNames(approx(sp$y, sp$x, xout = qnorm(pvs))$y,
                    sprintf("%.2f pct.", 100 * pvs))

  return(list(confs = confs, xs = xs, p_log_Lik = pls, data = out))
}
