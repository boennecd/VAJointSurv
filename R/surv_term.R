#' @export
surv_term <- function(formula, id, data, time_fixef){
  # get the input data
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  id <- eval(substitute(id), data)
  Z <- model.matrix(mt, mf)
  y <- model.response(mf)

  stopifnot(inherits(y, "Surv"),
            attr(y, "type") %in% c("counting", "right"))
  if(attr(y, "type") == "right")
    y <- cbind(0, y)
  stopifnot(all(y[, 3] %in% 0:1),
            any(y[, 3] == 0),
            any(y[, 3] == 1))

  time_fixef <- eval(substitute(time_fixef), data[y[, 3] == 1, ])

  # sanity checks
  stopifnot(NROW(Z) == length(id),
            is.matrix(Z),
            all(is.finite(Z)),
            all(is.finite(y)),
            NROW(Z) == NROW(y))
  is_valid_expansion(time_fixef)

  # reorder the data and return
  ord <- order(id, y[, 2])
  y <- y[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]
  id <- id[ord]

  structure(list(y = y, Z = t(Z), time_fixef = time_fixef, id = id, mt = mt),
            class = "surv_term")
}

# computes the starting values for the fixed effects
#' @importFrom stats optim
#' @importFrom utils head tail
surv_term_start_value <- function(object, quad_rule, va_var){
  stopifnot(inherits(object, "surv_term"))

  # create the object to perform the optimization
  Z <- object$Z
  surv <- t(object$y)
  comp <- ph_ll(time_fixef = object$time_fixef, Z = Z, surv = surv)

  # define the likelihood and gradient function. Then optimize
  ll <- function(par)
    ph_eval(ptr = comp$ptr, par = par, quad_rule = quad_rule, va_var = va_var)
  gr <- function(par)
    ph_grad(ptr = comp$ptr, par = par, quad_rule = quad_rule, va_var = va_var)

  # TODO: find a better starting value
  opt <- optim(numeric(comp$n_params), ll, gr, method = "BFGS")
  if(opt$convergence != 0)
    warning(sprintf(
      "optim returned convergence code %d when finding the starting values for the survival parameters",
      opt$convergence))

  list(fixef = head(opt$par, NROW(Z)), fixef_vary = tail(opt$par, -NROW(Z)))
}
