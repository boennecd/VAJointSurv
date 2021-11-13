wrap_term <- function(term){
  ptr <- expansion_object(term)
  term$ptr <- ptr
  term$eval <- function(x, der = 0, lower_limit = 0)
    eval_expansion(ptr, x, der, lower_limit)
  term
}

#' Term for Orthogonal Polynomials
#'
#' @param x,degree,coefs,raw same as \code{\link{poly}}.
#' @param intercept \code{TRUE} if there should be an intercept.
#' @param use_log \code{TRUE} if the polynomials should be in the log of the
#' argument.
#'
#' @importFrom stats poly
#' @export
poly_term <- function(x = numeric(), degree = 1, coefs = NULL, raw = FALSE,
                      intercept = FALSE, use_log = FALSE){
  out <- if(is.null(coefs))
    if(degree > 0 && !raw)
      list(coefs = attr(poly(if(use_log) log(x) else x, degree = degree),
                        "coefs"))
    else
      list(coefs = list(alpha = numeric(degree), norm2 = rep(1, 2 + degree)))
  else list(coefs = coefs)

  out[c("time", "intercept", "raw", "use_log")] <-
    list(x, intercept, raw, use_log)
  wrap_term(structure(out, class = "poly_term"))
}

#' Term for a Basis Matrix for Natural Cubic Splines
#'
#' @param x,df,knots,intercept,Boundary.knots same as \code{\link{ns}}.
#' @param use_log \code{TRUE} if the polynomials should be in the log of the
#' argument.
#'
#' @importFrom splines ns
#' @export
ns_term <- function(x = numeric(), df = NULL, knots = NULL, intercept = FALSE,
                    Boundary.knots = range(if(use_log) log(x) else x),
                    use_log = FALSE){
  out <- if(is.null(knots)){
    tmp <- ns(if(use_log) log(x) else x, df = df, knots = knots,
              Boundary.knots = Boundary.knots, intercept = intercept)
    list(knots = attr(tmp, "knots"))
  }
  else list(knots = knots)

  out[c("Boundary.knots", "time", "degree", "intercept", "use_log")] <-
    list(Boundary.knots, x, 3L, intercept, use_log)
  wrap_term(structure(out, class = "ns_term"))
}

#' Term for a B-Spline Basis for Polynomial Splines
#'
#' @param x,df,knots,degree,intercept,Boundary.knots same as \code{\link{bs}}.
#' @param use_log \code{TRUE} if the polynomials should be in the log of the
#' argument.
#'
#' @importFrom splines bs
#' @export
bs_term <- function(x = numeric(), df = NULL, knots = NULL, degree = 3,
                    intercept = FALSE,
                    Boundary.knots = range(if(use_log) log(x) else x),
                    use_log = FALSE){
  stopifnot(degree == 3)

  out <- if(is.null(knots)){
    tmp <- bs(if(use_log) log(x) else x, df = df, knots = knots,
              Boundary.knots = Boundary.knots, intercept = intercept)
    list(knots = attr(tmp, "knots"))

  } else list(knots = knots)

  out[c("Boundary.knots", "time", "degree", "intercept", "use_log")] <-
    list(Boundary.knots, x, degree, intercept, use_log)
  wrap_term(structure(out, class = "bs_term"))
}

is_valid_expansion <- function(x)
  stopifnot(inherits(x, c("poly_term", "ns_term", "bs_term")))
