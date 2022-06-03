wrap_term <- function(term){
  ptr <- expansion_object(term)
  term$ptr <- ptr
  term$eval <- function(x, der = 0, lower_limit = 0, newdata = NULL) {
    if(is.null(term$weights_symbol)) {
      weights <- matrix(0., nrow = 0L, ncol = length(x))
    } else {
      weight_call <- as.call(c(list(as.name("rbind")), term$weights_symbol))
      weights <- eval(weight_call,newdata,parent.frame())
    }
    eval_expansion(ptr, x, weights, der, lower_limit)
  }
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

  out[c("time", "intercept", "raw", "use_log", "weights_symbol")] <-
    list(x, intercept, raw, use_log, NULL)
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

  out[c("Boundary.knots", "time", "degree", "intercept", "use_log", "weights_symbol")] <-
    list(Boundary.knots, x, 3L, intercept, use_log, NULL)
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

  out[c("Boundary.knots", "time", "degree", "intercept", "use_log", "weights_symbol")] <-
    list(Boundary.knots, x, degree, intercept, use_log, NULL)
  wrap_term(structure(out, class = "bs_term"))
}

#' Term for a Basis Matrix for Weighted Term
#'
#' @description
#' Creates a weighted basis matrix where the entries are weighted with a
#' numeric vector to e.g. create a varying-coefficient.
#'
#' @param x a term type from the package.
#' @param weight a numeric vector with weights for x.
#'
#' @seealso
#' \code{\link{poly_term}}, \code{\link{bs_term}}, \code{\link{ns_term}}, and
#' \code{\link{stacked_term}}
#'
#' @examples
#' # TODO: add an example
#'
#' @export
weighted_term <- function(x, weight){
  term <- x
  weights_symbol <- substitute(weight)

  out <- list(term=term, weights_symbol = c(weights_symbol,term$weights_symbol))
  wrap_term(structure(out, class = "weighted_term"))
}

#' Term for a Basis Matrix for of Different Types of Terms
#'
#' @description
#' Creates a basis matrix consisting of different types of terms.
#' E.g. to create a varying-coefficient.
#'
#' @param ... term objects from the package.
#'
#' @seealso
#' \code{\link{poly_term}}, \code{\link{bs_term}}, \code{\link{ns_term}}, and
#' \code{\link{weighted_term}}
#'
#' @examples
#' # TODO: add an example
#'
#' @export
stacked_term <- function(...){
  if(...length() < 2)
    stop("stacked_term created with less than two arguments")

  terms <- list(...)

  out <- list(terms = terms, weights_symbol = unlist(lapply(terms,`[[`, 'weights_symbol')))
  wrap_term(structure(out, class ="stacked_term"))
}

is_valid_expansion <- function(x)
  stopifnot(inherits(
    x, c("poly_term", "ns_term", "bs_term", "weighted_term", "stacked_term")))
