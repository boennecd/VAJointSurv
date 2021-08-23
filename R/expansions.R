#' @importFrom stats poly
#' @export
poly_term <- function(x, degree = 1, coefs = NULL, raw = FALSE,
                      intercept = FALSE){
  out <- if(is.null(coefs))
    if(degree > 0)
      list(coefs = attr(poly(x, degree = degree), "coefs"))
    else
      list(coefs = list(alpha = numeric(), norm2 = rep(1, 2)))
  else list()

  out[c("time", "intercept", "raw")] <- list(x, intercept, raw)
  structure(out, class = "poly_term")
}

#' @importFrom splines ns
#' @export
ns_term <- function(x, df = NULL, knots = NULL, intercept = FALSE,
                    Boundary.knots = range(x)){
  out <- if(is.null(knots)){
    tmp <- ns(x, df = df, knots = knots, Boundary.knots = Boundary.knots,
              intercept = intercept)
    list(knots = attr(tmp, "knots"))
  }
  else list()

  out[c("Boundary.knots", "time", "degree", "intercept")] <-
    list(Boundary.knots, x, 3L, intercept)
  structure(out, class = "ns_term")
}

#' @importFrom splines bs
#' @export
bs_term <- function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                    Boundary.knots = range(x)){
  out <- if(is.null(knots)){
    tmp <- bs(x, df = df, knots = knots, Boundary.knots = Boundary.knots,
              intercept = intercept)
    list(knots = attr(tmp, "knots"))

  } else list()

  out[c("Boundary.knots", "time", "degree", "intercept")] <-
    list(Boundary.knots, x, degree, intercept)
  structure(out, class = "bs_term")
}

is_valid_expansion <- function(x)
  stopifnot(inherits(x, c("poly_term", "ns_term", "bs_term")))
