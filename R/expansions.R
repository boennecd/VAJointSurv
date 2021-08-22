#' @importFrom stats poly
#' @export
poly_term <- function(x, degree = 1, coefs = NULL, raw = FALSE,
                      intercept = FALSE){
  out <- if(is.null(coefs))
    list(coefs = attr(poly(x, degree = degree), "coefs"))
  else list()

  out$time <- x
  out$intercept <- intercept
  structure(out, "poly_term")
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

  out$Boundary.knots <- Boundary.knots
  out$time <- x
  out$degree <- 3L
  out$intercept <- intercept
  structure(out, "ns_term")
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

  out$Boundary.knots <- Boundary.knots
  out$time <- x
  out$degree <- degree
  out$intercept <- intercept
  structure(out, "bs_term")
}

is_valid_expansion <- function(x)
  stopifnot(inherits(x, c("poly_term", "ns_term", "bs_term")))
