# use Gauss-Legendre quadrature by default
#' @importFrom SimSurvNMarker get_gl_rule
default_quad_rule <- function(){
  rule <- get_gl_rule(50)
  within(rule, {
    node <- node / 2 + .5
    weight <- weight / 2
  })
}

# check a quadrature rule the integral from zero to one
check_quad_rule <- function(quad_rule)
  with(quad_rule, {
    stopifnot(is.list(quad_rule),
              is.numeric(node), all(is.finite(node)),
              all(node >= 0), all(node <= 1),
              is.numeric(weight), all(is.finite(weight)),
              length(node) == length(weight))
  })


#' Creates a joint_ms Object to Estimate a Joint Survival and Marker Model
#'
#' @param markers either an object from \code{\link{marker_term}} or a list
#' of such objects.
#' @param survival_terms either an object from \code{\link{surv_term}} or a list
#' of such objects.
#' @param max_threads maximum number of threads to use.
#' @param quad_rule list with nodes and weights for a quadrature rule for the
#' integral from zero to one.
#'
#' @export
joint_ms_ptr <- function(markers = list(), survival_terms = list(),
                         max_threads = 1L, quad_rule = NULL){
  # handle defaults
  if(inherits(markers, "marker_term"))
    markers <- list(markers)
  else
    stopifnot(all(sapply(markers, inherits, "marker_term")))
  if(inherits(survival_terms, "surv_term"))
    survival_terms <- list(survival_terms)
  else
    stopifnot(all(sapply(survival_terms, inherits, "surv_term")))

  if(is.null(quad_rule))
    quad_rule <- default_quad_rule()
  check_quad_rule(quad_rule)

  # handle the default for ders
  survival_terms <- lapply(survival_terms, function(x){
    if(is.null(x$ders))
      x$ders <- replicate(length(markers), 0L, simplify = FALSE)
    x
  })

  ptr <- .joint_ms_ptr(markers, survival_terms, max_threads = max_threads)
  param_names <- joint_ms_parameter_names(ptr)
  out <- list(param_names = param_names, ptr = ptr)
  indices <- joint_ms_parameter_indices(ptr)

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
  if(length(survival_terms) > 0){
    Xi <- diag(def_frailty_var, length(survival_terms))
    start_val[indices$vcovs$vcov_surv] <- .log_chol(Xi)
  }

  s_start_val <- lapply(survival_terms, surv_term_start_value,
                        quad_rule = quad_rule, va_var = def_frailty_var)
  for(i in seq_along(survival_terms)){
    start_val[indices$survival[[i]]$fixef] <- s_start_val[[i]]$fixef
    start_val[indices$survival[[i]]$fixef_vary] <- s_start_val[[i]]$fixef_vary
  }

  # fill in default values for the VA parameters
  va_dim <- indices$va_dim
  va_vcov_default <- diag(va_dim)
  if(length(indices$survival) > 0)
    diag(va_vcov_default)[va_dim:(va_dim - length(indices$survival) + 1L)] <-
    def_frailty_var

  va_default <- c(numeric(va_dim), .log_chol(va_vcov_default))
  start_val[-(1:(indices$va_params_start - 1L))] <- va_default

  structure(
    list(param_names = param_names, indices = indices, ptr = ptr,
         start_val = start_val, max_threads = max_threads,
         quad_rule = quad_rule, n_lb_terms = joint_ms_n_terms(ptr)),
    class = "joint_ms")
}

#' Quick Heuristic for the Starting Values
#'
#' @inheritParams joint_ms_opt
#' @param rel_eps,c1,c2,max_it arguments to pass to the C++ version of
#' \code{\link{psqn}}.
#'
#' @importFrom psqn psqn
#'
#' @export
joint_ms_start_val <- function(object, n_threads = object$max_threads,
                               rel_eps = 1e-8, c1 = 1e-4, c2 = .9,
                               max_it = 100L, quad_rule = object$quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  if(is.null(quad_rule))
    quad_rule <- default_quad_rule()
  check_quad_rule(quad_rule)

  opt_priv(val = object$start_val, ptr = object$ptr, rel_eps = rel_eps,
           max_it = max_it, n_threads = n_threads, c1 = c1, c2 = c2,
           quad_rule = quad_rule)
}

#' Evaluates the Lower Bound or the Gradient of the Lower Bound
#'
#' @inheritParams joint_ms_ptr
#' @param object a joint_ms object from \code{\link{joint_ms_ptr}}.
#' @param par parameter vector the lower bound is evaluated at.
#' @param n_threads number of threads to use.
#'
#' @export
joint_ms_lb <- function(object, par, n_threads = object$max_threads,
                        quad_rule = object$quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  if(is.null(quad_rule))
    quad_rule <- default_quad_rule()
  check_quad_rule(quad_rule)

  joint_ms_eval_lb(val = par, ptr = object$ptr, n_threads = n_threads,
                   quad_rule = quad_rule)
}

#' @rdname joint_ms_lb
#' @export
joint_ms_lb_gr <- function(object, par, n_threads = object$max_threads,
                           quad_rule = object$quad_rule){
  stopifnot(inherits(object, "joint_ms"))

  if(is.null(quad_rule))
    quad_rule <- default_quad_rule()
  check_quad_rule(quad_rule)

  joint_ms_eval_lb_gr(val = par, ptr = object$ptr, n_threads = n_threads,
                      quad_rule = quad_rule)
}

#' Optimizes the Lower Bound
#'
#' @inheritParams joint_ms_lb
#' @param par starting value.
#' @param rel_eps,max_it,c1,c2,use_bfgs,trace,cg_tol,strong_wolfe,max_cg,pre_method,mask
#' arguments to pass to the C++ version of \code{\link{psqn}}.
#' @export
joint_ms_opt <- function(object, par = object$start_val, rel_eps = 1e-8,
                         max_it = 100L, n_threads = object$max_threads,
                         c1 = 1e-4, c2 = .9, use_bfgs = TRUE, trace = 0L,
                         cg_tol = .5, strong_wolfe = TRUE, max_cg = 0L,
                         pre_method = 1L, quad_rule = object$quad_rule,
                         mask = integer()){
  stopifnot(inherits(object, "joint_ms"))
  if(is.null(quad_rule))
    quad_rule <- default_quad_rule()
  check_quad_rule(quad_rule)
  stopifnot(is.integer(mask), all(mask >= 0 & mask < length(par)))

  joint_ms_opt_lb(val = par, ptr = object$ptr, rel_eps = rel_eps,
                  max_it = max_it, n_threads = n_threads, c1 = c1, c2 = c2,
                  use_bfgs = use_bfgs, trace = trace, cg_tol = cg_tol,
                  strong_wolfe = strong_wolfe, max_cg = max_cg,
                  pre_method = pre_method, quad_rule = quad_rule,
                  mask = mask)
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
#' @importFrom stats setNames
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
#'
#' @importFrom stats approx qchisq splinefun qnorm
#'
#' @export
joint_ms_profile <- function(
  object, opt_out, which_prof, delta, level = .95, max_step = 15L,
  rel_eps = 1e-8, max_it = 100L, n_threads = object$max_threads, c1 = 1e-04,
  c2 = 0.9, use_bfgs = TRUE, trace = 0L, cg_tol = 0.5, strong_wolfe = TRUE,
  max_cg = 0L, pre_method = 1L, quad_rule = object$quad_rule, verbose = TRUE){
  stopifnot(is.integer(which_prof), length(which_prof) == 1L,
            which_prof >= 1L, which_prof <= length(opt_out$par),
            is.numeric(delta), is.finite(delta), length(delta) == 1,
            delta > 0,
            length(level) == 1, is.numeric(level), is.finite(level),
            level > 0, level < 1,
            opt_out$convergence)

  # setup
  mask <- which_prof - 1L
  n_global <- comp_obj$indices$va_params_start - 1L
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
    list(x = x, value = -fit$value, par = head(fit$par, n_global),
         z_val = sign(dir) * sqrt((max_lb + fit$value) * 2))
  }

  # fits the model
  fit_model <- function(x, dir, lb = NULL, ub = NULL){
    par[which_prof] <- x
    fit <- joint_ms_opt(
      object = object, par = par, rel_eps = rel_eps, max_it = max_it,
      n_threads = n_threads, c1 = c1, c2 = c2, use_bfgs = use_bfgs,
      trace = trace, cg_tol = cg_tol, strong_wolfe = strong_wolfe,
      max_cg = max_cg, pre_method = pre_method, quad_rule = quad_rule,
      mask = mask)
    format_fit(x, dir = dir, fit = fit, lb = lb, ub = ub)
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

    while(prev > crit_value && (step <- step + 1L) <= max_step){
      out[[step]] <- fit_model(par[which_prof] + dir *  2^(step - 1) * delta,
                               dir = dir)
      if(out[[step]]$value > prev){
        warning("The lower bound did not decrease. Either the input is not an optimum or the precision needs to be increased")
        did_fail <- TRUE
        break
      }

      prev <- out[[step]]$value
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

      out[[step]] <- fit_model(next_val, dir = dir, lb = lb, ub = ub)

      if(out[[step]]$value > ub$value || out[[step]]$value < lb$value){
        warning("Log likelihood does not seem monotonic. Likely the precision needs to be increased")
        break
      }

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
  res_down <- get_points(-1)
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
