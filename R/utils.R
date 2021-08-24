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
