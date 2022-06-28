# compute the log Cholesky decomposition
.log_chol <- function(x){
  L <- chol(x)
  diag(L) <- log(diag(L))
  L[upper.tri(L, TRUE)]
}

# computes the inverse of the log Cholesky decomposition
.log_chol_inv <- function(x){
  n <- round((sqrt(8 * length(x) + 1) - 1) / 2)
  out <- matrix(0, n, n)
  out[upper.tri(out, TRUE)] <- x
  diag(out) <- exp(diag(out))
  crossprod(out)
}

#' Extracts the Variational Parameters
#'
#' @description
#' Computes the estimated variational parameters for each individual.
#'
#' @inheritParams joint_ms_format
#'
#' @returns
#' A list with one list for each individual with the estimated mean and
#' covariance matrix.
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
#' # extract variational parameters for each individual
#' VA_pars <- joint_ms_va_par(object = model_ptr,par = start_vals)
#'
#' # number of sets of variational parameters is equal to the number of subjects
#' length(VA_pars)
#' length(unique(pbc$id))
#'
#' # mean and var-covar matrix for 1st individual
#' VA_pars[[1]]
#' @export
joint_ms_va_par <- function(object, par = object$start_val){
  stopifnot(inherits(object, "joint_ms"))
  va_params_start <- object$indices$va_params_start
  va_dim <- object$indices$va_dim

  dim_out <- (va_dim * (va_dim + 3L)) / 2L
  vcov_dim <- (va_dim * (va_dim + 1L)) / 2L
  n_ids <- (length(par) - va_params_start + 1L) / dim_out
  va_par <- lapply(
    1:n_ids - 1L, function(idx){
      start <- idx * dim_out + va_params_start
      mean <- par[start + 0:(va_dim - 1L)]
      vcov <- .log_chol_inv(par[start + va_dim + 0:(vcov_dim - 1L)])
      list(mean = mean, vcov = vcov)
    })

  setNames(va_par, object$ids)
}
