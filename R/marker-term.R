#' @export
marker_term <- function(formula, id, data, time_fixef, time_rng){
  # get the input data
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  id <- eval(substitute(id), data)
  X <- model.matrix(mt, mf)
  y <- model.response(mf, "numeric")

  time_fixef <- eval(substitute(time_fixef), data)
  time_rng <- eval(substitute(time_rng), data)

  # sanity checks
  stopifnot(NROW(X) == length(id),
            NROW(X) == length(time_var),
            !is.matrix(y),
            is.matrix(X),
            all(is.finite(X)),
            NROW(X) == length(y),
            NROW(X) == length(time_fixef$time),
            all(is.finite(time_fixef$time)),
            NROW(X) == length(time_rng$time),
            all(is.finite(time_rng$time)),
            all(time_rng$time == time_fixef$time))
  is_valid_expansion(time_fixef)
  is_valid_expansion(time_rng)

  # prepare the data to return
  time_var <-  time_fixef$time
  time_fixef$time <- NULL
  time_rng$time <- NULL

  # we have to sort by id and type
  ord <- order(id, time_var)
  id <- id[ord]
  time_var <- time_var[ord]
  X <- t(X)[, ord]
  y <- y[ord]

  structure(list(time = time_var, X = X, y = y, id = id, mt = mt,
                 time_fixef = time_fixef, time_rng = time_rng),
            "marker_term")
}
