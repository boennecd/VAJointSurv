#' @title Expansions in the VAJointSurv package
#' @name VAJointSurv-terms
#'
#' @description
#' The VAJointSurv package uses different functions to allow for expansions in
#' time possibly with covariate interactions. The main usage of the functions
#' is internally but they do provide an element called `eval()` which is a
#' function to evaluate the expansion. These functions take the following
#' arguments:
#'
#' \itemize{
#' \item \code{x} numeric vector with points at which to evaluate the expansion.
#' \item \code{der} integer indicating whether to evaluate the expansion, its integral,
#' or the derivative.
#' \item \code{lower_limit} possible lower limit if integration is performed.
#' \item \code{newdata} a \code{\link{data.frame}} with new data if this is required.
#' E.g. for \code{\link{weighted_term}}.
#' }
#'
#' The supported terms are \code{\link{ns_term}}, \code{\link{bs_term}},
#' \code{\link{poly_term}}, \code{\link{weighted_term}}, and a
#' \code{\link{stacked_term}}.
NULL

wrap_term <- function(term){
  ptr <- expansion_object(term)
  term$ptr <- ptr
  term$eval <- function(x, der = 0, lower_limit = 0, newdata = NULL) {
    weights <- bases_weights(term$weights_symbol,newdata,parent.frame(),length(x))
    eval_expansion(ptr, x, weights, der, lower_limit)
  }
  term
}

bases_weights <- function(weights_symbol, newdata, enclos, expected_length) {
  if(is.null(weights_symbol)) {
    weights <- matrix(0., nrow = 0L, ncol = expected_length)
  } else {
    weight_call <- as.call(c(list(as.name("rbind")),weights_symbol))
    weights <- eval(weight_call,newdata,enclos)
  }
  stopifnot(ncol(weights) == expected_length,is.numeric(weights))
  weights
}

#' Term for Orthogonal Polynomials
#'
#' @param x,degree,coefs,raw same as \code{\link{poly}}.
#' @param intercept \code{TRUE} if there should be an intercept.
#' @param use_log \code{TRUE} if the polynomials should be in the log of the
#' argument.
#'
#' @return
#' A list like \code{\link{poly}} with an additional element called \code{eval}
#' to evaluate the basis. See \code{\link{VAJointSurv-terms}}.
#'
#' @seealso
#' \code{\link{bs_term}}, \code{\link{ns_term}}, \code{\link{weighted_term}}, and
#' \code{\link{stacked_term}}.
#'
#' @importFrom stats poly
#' @examples
#' vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
#' spline_basis <- poly_term(vals,degree = 3, raw = TRUE)
#' # evaluate spline basis at 0.5
#' spline_basis$eval(0.5)
#' # evaluate first derivative of spline basis at 0.5
#' spline_basis$eval(0.5, der = 1)
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
#' @return
#' A list like \code{\link{ns}} with an additional element called \code{eval}
#' to evaluate the basis. See \code{\link{VAJointSurv-terms}}.
#'
#' @seealso
#' \code{\link{poly_term}}, \code{\link{bs_term}}, \code{\link{weighted_term}}, and
#' \code{\link{stacked_term}}.
#'
#' @importFrom splines ns
#' @examples
#' vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
#' spline_basis <- ns_term(vals,df = 3)
#' # evaluate spline basis at 0.5
#' spline_basis$eval(0.5)
#' # evaluate first derivative of spline basis at 0.5
#' spline_basis$eval(0.5, der = 1)
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
#' @return
#' A list like \code{\link{bs}} with an additional element called \code{eval}
#' to evaluate the basis. See \code{\link{VAJointSurv-terms}}.
#'
#' @seealso
#' \code{\link{poly_term}}, \code{\link{ns_term}}, \code{\link{weighted_term}},
#' and \code{\link{stacked_term}}.
#'
#' @importFrom splines bs
#' @examples
#' vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
#' spline_basis <- bs_term(vals,df = 3)
#' # evaluate spline basis at 0.5
#' spline_basis$eval(0.5)
#' # evaluate first derivative of spline basis at 0.5
#' spline_basis$eval(0.5, der = 1)
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
#' @param weight a symbol for the weight. Notice that the symbol is first
#' first used when the \code{eval} function on the returned object is called.
#'
#' @return
#' A list with an element called \code{eval}
#' to evaluate the basis. See \code{\link{VAJointSurv-terms}}.
#'
#' @seealso
#' \code{\link{poly_term}}, \code{\link{bs_term}}, \code{\link{ns_term}}, and
#' \code{\link{stacked_term}}.
#'
#' @examples
#' vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
#'
#' spline_basis <- ns_term(vals, df = 3)
#' ws <- c(4,5)
#' # create a weighted term
#' w_term <- weighted_term(spline_basis, weights)
#'
#' # evaluate weighted basis at 0.5 and 0.7 with weights 4 and 5
#' w_term$eval(c(0.5,0.7), newdata = data.frame(weights = ws))
#' # evaluate the first derivative of weighted basis at 0.5 and 0.7
#' # with weights 4 and 5
#' w_term$eval(c(0.5,0.7), newdata = data.frame(weights = ws), der = 1)
#' @export
weighted_term <- function(x, weight){
  is_valid_expansion(x)
  term <- x
  weights_symbol <- substitute(weight)

  out <- list(
    term = term, weights_symbol = c(weights_symbol,term$weights_symbol),
    time = term$time)
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
#' @return
#' A list with an element called \code{eval}
#' to evaluate the basis. See \code{\link{VAJointSurv-terms}}.
#'
#' @seealso
#' \code{\link{poly_term}}, \code{\link{bs_term}}, \code{\link{ns_term}}, and
#' \code{\link{weighted_term}}.
#'
#' @examples
#' vals <- c(0.41, 0.29, 0.44, 0.1, 0.18, 0.65, 0.29, 0.85, 0.36, 0.47)
#'
#' spline_basis1 <- ns_term(vals, df = 3)
#' spline_basis2 <- bs_term(vals, df = 3)
#'
#' # create stacked term from two spline bases
#' stacked_basis <- stacked_term(spline_basis1, spline_basis2)
#'
#' # evaluate stacked basis at 0.5
#' stacked_basis$eval(0.5)
#' # evaluate first derivative of stacked basis at 0.5
#' stacked_basis$eval(0.5, der = 1)
#' @export
stacked_term <- function(...){
  if(...length() < 2)
    stop("stacked_term created with less than two arguments")

  terms <- list(...)
  for (x in terms) is_valid_expansion(x)
  stopifnot(all(sapply(terms, function(x) identical(x$time, terms[[1]]$time))))

  out <- list(
    terms = terms,
    weights_symbol = unlist(lapply(terms,`[[`, 'weights_symbol')),
    time = terms[[1]]$time)
  wrap_term(structure(out, class ="stacked_term"))
}

is_valid_expansion <- function(x)
  stopifnot(inherits(
    x, c("poly_term", "ns_term", "bs_term", "weighted_term", "stacked_term")))
