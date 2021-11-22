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
