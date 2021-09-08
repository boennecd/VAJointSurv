#' Term for Orthogonal Polynomials
#'
#' @param x,degree,coefs,raw same as \code{\link{poly}}.
#' @param intercept \code{TRUE} if there should be an intercept.
#'
#' @importFrom stats poly
#' @export
poly_term <- function(x = numeric(), degree = 1, coefs = NULL, raw = FALSE,
                      intercept = FALSE){
  out <- if(is.null(coefs))
    if(degree > 0 && !raw)
      list(coefs = attr(poly(x, degree = degree), "coefs"))
    else
      list(coefs = list(alpha = numeric(degree), norm2 = rep(1, 2 + degree)))
  else list(coefs = coefs)

  out[c("time", "intercept", "raw")] <- list(x, intercept, raw)
  structure(out, class = "poly_term")
}

#' Term for a Basis Matrix for Natural Cubic Splines
#'
#' @param x,df,knots,intercept,Boundary.knots same as \code{\link{ns}}.
#'
#' @importFrom splines ns
#' @export
ns_term <- function(x = numeric(), df = NULL, knots = NULL, intercept = FALSE,
                    Boundary.knots = range(x)){
  out <- if(is.null(knots)){
    tmp <- ns(x, df = df, knots = knots, Boundary.knots = Boundary.knots,
              intercept = intercept)
    list(knots = attr(tmp, "knots"))
  }
  else list(knots = knots)

  out[c("Boundary.knots", "time", "degree", "intercept")] <-
    list(Boundary.knots, x, 3L, intercept)
  structure(out, class = "ns_term")
}

#' Term for a B-Spline Basis for Polynomial Splines
#'
#' @param x,df,knots,degree,intercept,Boundary.knots same as \code{\link{bs}}.
#'
#' @importFrom splines bs
#' @export
bs_term <- function(x = numeric(), df = NULL, knots = NULL, degree = 3,
                    intercept = FALSE, Boundary.knots = range(x)){
  out <- if(is.null(knots)){
    tmp <- bs(x, df = df, knots = knots, Boundary.knots = Boundary.knots,
              intercept = intercept)
    list(knots = attr(tmp, "knots"))

  } else list(knots = knots)

  out[c("Boundary.knots", "time", "degree", "intercept")] <-
    list(Boundary.knots, x, degree, intercept)
  structure(out, class = "bs_term")
}

is_valid_expansion <- function(x)
  stopifnot(inherits(x, c("poly_term", "ns_term", "bs_term")))
