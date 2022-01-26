# settings for the simulation
library(splines)
g_funcs <- list(
  function(x)
    ns(x, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  function(x)
    # a raw polynomial
    outer(x, 1:2, `^`))
m_funcs <- list(
  function(x)
    ns(x, knots = numeric(), Boundary.knots = c(0, 10), intercept = TRUE),
  function(x)
    # a raw polynomial
    outer(x, 0:1, `^`))

fixef_vary_marker <- list(c(1.4, 1.2, -2.1), c(.5, -.02)) # beta
fixef_marker <- list(c(-.5, 2), 1) # gamma

# Psi
vcov_vary <- structure(c(0.35, 0.08, -0.05, 0.01, 0.08, 1.92, -0.24, -0.04,
                         -0.05, -0.24, 0.32, 0.09, 0.01, -0.04, 0.09, 0.12),
                       .Dim = c(4L, 4L))
vcov_marker <- matrix(c(.6^2, .1, .1, .4^2), 2)

# the survival parameters
vcov_surv <- matrix(c(.2^2, .15^2, .15^2, .25^2), 2) # Xi

fixef_surv <- list(c(-1, .25), .2)
associations <- list(c(.6, -.4), c(-.7, .2))
fixef_vary_surv <- list(c(.5, .1, -.2, .11),
                        c(-1, -.25))

b_funcs <- list(
  function(x) bs(x, knots = 5, Boundary.knots = c(0, 10)),
  function(x) ns(x, knots = 5, Boundary.knots = c(0, 10)))

library(SimSurvNMarker)

## ----delayed_sim_fit_obs_process_markers_and_recurrent, cache = 1------------------------------------------------------------------------
library(mvtnorm)
# simulates from the model by sampling a given number of individuals
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # sample the delayed entry time
    delayed_entry <- pmax(runif(1, -1, 5), 0)

    # sample the censoring time
    cens <- -Inf
    while(cens < delayed_entry)
      cens <- min(rexp(1, rate = 1/10), 10)

    # sample the terminal event time and the random effects
    Z1 <- c(1, runif(1, -1, 1))

    y_terminal <- -Inf
    while(y_terminal < delayed_entry){
      U <- drop(rmvnorm(1, sigma = vcov_vary))
      frailties <- drop(rmvnorm(1, sigma = vcov_surv))
      log_haz_offset <- sum(Z1 * fixef_surv[[1]]) + frailties[1]

      # assign the conditional hazard function
      expansion <- function(x, b_func)
        cbind(b_func(x), m_funcs[[1]](x) %*% U[1:2],
              m_funcs[[2]](x) %*% U[3:4])
      surv_func <- function(ti, fixef_vary_surv, associations, b_func){
        formals(expansion)$b_func <- b_func
        eval_surv_base_fun(
          ti = ti, omega = c(fixef_vary_surv, associations), b_func = expansion,
          gl_dat = gl_dat, delta = log_haz_offset)
      }

      # sample the survival time
      rng <- runif(1)
      root_func <- function(x, rng)
        rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[1]],
                        associations = associations[[1]], b_func = b_funcs[[1]])

      if(root_func(cens, rng) < 0){
        # the observation is censored
        y_terminal <- cens
        event <- 0
      } else {
        # find the event time
        root <- uniroot(root_func, c(0, cens), tol = 1e-6, rng = rng)
        y_terminal <- root$root
        event <- 1

      }
    }

    terminal_outcome <- cbind(y = y_terminal, event = event, Z1 = Z1[2],
                              id = id, delayed_entry = delayed_entry)

    # clean up
    rm(list = setdiff(ls(), c(
      "y_terminal", "terminal_outcome", "expansion", "surv_func", "frailties",
      "U", "id", "delayed_entry")))

    # simulate the observation times
    Z2 <- 1
    log_haz_offset <- sum(Z2 * fixef_surv[[2]]) + frailties[2]

    root_func <- function(x, left_trunc_surv, rng)
      rng - surv_func(x, fixef_vary_surv = fixef_vary_surv[[2]],
                      associations = associations[[2]], b_func = b_funcs[[2]]) /
      left_trunc_surv

    max_sample <- 1000L
    left_trunc_surv <- 1
    Z2 <- matrix(rep(Z2, each = max_sample), max_sample)
    event <- y <- lf_trunc <- rep(NA_real_, max_sample)
    lf_trunc_i <- 0
    for(i in 1:max_sample){
      # sample a random uniform variable and invert the survival function
      rng_i <- runif(1)
      lf_trunc[i] <- lf_trunc_i

      if(root_func(y_terminal, left_trunc_surv, rng_i) < 0){
        # the observation is right-censored and we can exit
        y[i] <- y_terminal
        event[i] <- 0
        break
      }

      # we need to invert the survival function to find the observation time
      root <- uniroot(root_func, c(lf_trunc_i, y_terminal), tol = 1e-6,
                      left_trunc_surv = left_trunc_surv, rng = rng_i)
      lf_trunc_i <- y[i] <- root$root
      event[i] <- 1
      left_trunc_surv <- surv_func(
        y[i], fixef_vary_surv = fixef_vary_surv[[2]], associations = associations[[2]],
        b_func = b_funcs[[2]])
    }

    colnames(Z2) <- paste0("Z", 1:NCOL(Z2) - 1L)
    obs_process <- cbind(lf_trunc = lf_trunc[1:i], y = y[1:i],
                         event = event[1:i], Z2[1:i, -1, drop = FALSE],
                         id = id)

    # account for the delayed entry
    obs_process[, "lf_trunc"] <- pmax(delayed_entry, obs_process[, "lf_trunc"])
    obs_process <- obs_process[
      obs_process[, "y"] > delayed_entry, , drop = FALSE]

    # clean up
    rm(list = setdiff(ls(), c("terminal_outcome", "U", "id",
                              "obs_process", "delayed_entry")))

    # sample the number of outcomes and the fixed effect covariates
    obs_time <- c(delayed_entry, obs_process[obs_process[, "event"] == 1, "y"])

    n_obs <- length(obs_time)
    X1 <- cbind(1, rnorm(n_obs))
    X2 <- matrix(1, n_obs)
    colnames(X1) <- paste0("X1_", 1:NCOL(X1) - 1L)
    colnames(X2) <- paste0("X2_", 1:NCOL(X2) - 1L)
    X <- list(X1, X2)

    # sample the outcomes
    eta <- sapply(1:2, function(i)
      X[[i]] %*% fixef_marker[[i]] +
        drop(g_funcs[[i]](obs_time)) %*% fixef_vary_marker[[i]] +
        drop(m_funcs[[i]](obs_time)) %*% U[1:2 + (i == 2) * 2])

    ys <- eta + rmvnorm(n_obs, sigma = vcov_marker)
    colnames(ys) <- paste0("Y", 1:2)

    # mask some observations
    do_mask <- sample.int(3L, n_obs, replace = TRUE)
    ys[do_mask == 2, 1] <- NA
    ys[do_mask == 3, 2] <- NA

    X <- do.call(cbind, lapply(X, function(x) x[, -1, drop = FALSE]))
    marker_data <- cbind(ys, X, time = obs_time, id = id)

    return(list(marker_data = marker_data, obs_process = obs_process,
                terminal_outcome = terminal_outcome))
  })

  # combine the data and return
  marker_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "marker_data")))
  marker_data$id <- as.integer(marker_data$id)
  # the order does not matter
  marker_data <- marker_data[sample.int(NROW(marker_data)), ]

  obs_process <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "obs_process")))
  obs_process$id <- as.integer(obs_process$id)
  # the order does not matter
  obs_process <- obs_process[sample.int(NROW(obs_process)), ]

  terminal_outcome <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "terminal_outcome")))
  terminal_outcome$id <- as.integer(terminal_outcome$id)
  # the order does not matter
  terminal_outcome <- terminal_outcome[sample.int(NROW(terminal_outcome)), ]

  list(marker_data = marker_data, obs_process = obs_process,
       terminal_outcome = terminal_outcome)
}

# sample a moderate sized data set
set.seed(4)
dat <- sim_dat(100)
save.image("bench-delayed.RData")
