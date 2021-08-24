#' Creates a joint_ms to Estimate a Joint Survival and Marker Model
#' @export
joint_ms_ptr <- function(markers, max_threads = 1L){
  if(inherits(markers, "marker_term"))
    markers <- list(markers)
  else
    stopifnot(all(sapply(markers, inherits, "marker_term")))

  ptr <- .joint_ms_ptr(markers, max_threads = max_threads)
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

  start_val[indices$vcovs$vcov_marker] <- .log_chol(sigma)

  # TODO: set the starting values for the survival outcomes

  structure(
    list(param_names = param_names, indices = indices, ptr = ptr,
         start_val = start_val, max_threads = max_threads),
    class = "joint_ms")
}

#' Quick Heuristic for the Starting Values
#' @export
joint_ms_start_val <- function(object, n_threads = object$max_threads,
                               rel_eps = 1e-8, c1 = 1e-4, c2 = .9,
                               max_it = 100L){
  stopifnot(inherits(object, "joint_ms"))
  opt_priv(val = object$start_val, ptr = object$ptr, rel_eps = rel_eps,
           max_it = max_it, n_threads = n_threads, c1 = c1, c2 = c2)
}

#' Evaluates the Lower Bound
#' @export
joint_ms_lb <- function(object, par, n_threads = object$max_threads){
  stopifnot(inherits(object, "joint_ms"))
  joint_ms_eval_lb(val = par, ptr = object$ptr, n_threads = n_threads)
}

#' Evaluates the Gradient of the Lower Bound
#' @export
joint_ms_lb_gr <- function(object, par, n_threads = object$max_threads){
  stopifnot(inherits(object, "joint_ms"))
  joint_ms_eval_lb_gr(val = par, ptr = object$ptr, n_threads = n_threads)
}

#' Optimizes the Lower Bound
#' @export
joint_ms_opt <- function(object, par = object$start_val, rel_eps = 1e-8,
                         max_it = 100L, n_threads = object$max_threads,
                         c1 = 1e-4, c2 = .9, use_bfgs = TRUE, trace = 0L,
                         cg_tol = .5, strong_wolfe = TRUE, max_cg = 0L,
                         pre_method = 1L){
  stopifnot(inherits(object, "joint_ms"))

  joint_ms_opt_lb(val = par, ptr = object$ptr, rel_eps = rel_eps,
                  max_it = max_it, n_threads = n_threads, c1 = c1, c2 = c2,
                  use_bfgs = use_bfgs, trace = trace, cg_tol = cg_tol,
                  strong_wolfe = strong_wolfe, max_cg = max_cg,
                  pre_method = pre_method)
}

#' Formats the Parameter Vector
#' @importFrom stats setNames
#' @export
joint_ms_format <- function(object, par = object$start_val){
  # TODO: add tests for this function
  stopifnot(inherits(object, "joint_ms"))
  indices <- object$indices

  # handle the markers
  markers <- lapply(seq_along(indices$markers), function(i)
    list(fixef = par[indices$markers[[i]]$fixef],
         fixef_vary <- par[indices$markers[[i]]$fixef_vary]))

  # handle the survival outcomes
  survival <- lapply(seq_along(indices$survival), function(i)
    list(fixef = par[indices$survival[[i]]$fixef],
         fixef_vary <- par[indices$survival[[i]]$fixef_vary],
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
