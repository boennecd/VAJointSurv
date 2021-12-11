options(digits = 4)

# settings for the simulation. You can skip this
library(VAJointSurv)
library(splines)
g_basis <- ns_term(knots = c(3, 5), Boundary.knots = c(1, 8))
g_funcs <- \(x) ns(x, knots = c(3, 5), Boundary.knots = c(1, 8))
m_basis <- ns_term(knots = numeric(), Boundary.knots = c(1, 8),
                   intercept = TRUE)
m_funcs <- \(x) ns(x, knots = numeric(), Boundary.knots = c(1, 8),
                   intercept = TRUE)

fixef_vary_marker <- c(1.4, -1.2, .4) # beta
fixef_marker <- c(-.5, 1) # gamma

vcov_vary <- structure(c(0.35, 0.004, 0.004, 0.048), .Dim = c(2L, 2L))
vcov_marker <- matrix(1, 1)

# the survival parameters
fixef_surv <- c(-3.5, .4)
association <- -2
fixef_vary_surv <- c(.5, .1, -.015)
fvar <- matrix(1e-6^2, 1) # not supported by the other models

b_basis <- poly_term(degree = 3, raw = TRUE)
b_func <- \(x) cbind(x, x^2, x^3)

# plot the population marker curve
par(mar = c(5, 5, 1, 1), cex = 1.2)
plot_marker(
  time_fixef = g_basis, time_rng = m_basis,
  fixef_vary = fixef_vary_marker, x_range = c(0, 10),
  vcov_vary = vcov_vary, ylab = "Marker 1")

# plot the conditional hazard
plot_surv(time_fixef = b_basis, time_rng = list(m_basis),
          x_range = c(0, 10), fixef_vary = fixef_vary_surv,
          vcov_vary = vcov_vary, frailty_var = fvar, ps = c(.1, .5, .9),
          log_hazard_shift = fixef_surv[1], associations = association)

rm(g_basis, m_basis, b_basis)

# assign function that simulates from the model by sampling a given number of
# individuals. You can skip this
library(mvtnorm)
library(SimSurvNMarker)
sim_dat <- function(n_ids){
  # simulate the outcomes
  gl_dat <- get_gl_rule(100L)
  dat <- lapply(1:n_ids, function(id){
    # draw the censoring time and the random effects
    cens <- min(10, rexp(1, 1/8))
    U <- drop(rmvnorm(1, sigma = vcov_vary))

    # simulate the survival outcome
    rng_surv <- rnorm(1, sd = sqrt(fvar))
    Z <- c(1, runif(1, -1, 1))
    log_haz_offset <- sum(Z * fixef_surv) + rng_surv

    expansion <- function(x)
      cbind(b_func(x), m_funcs(x) %*% U)

    # the conditional survival function
    surv_func <- function(ti)
      eval_surv_base_fun(
        ti = ti, omega = c(fixef_vary_surv, association), b_func = expansion,
        gl_dat = gl_dat, delta = log_haz_offset)

    # simulate the event
    rng_i <- runif(1)
    root_func <- function(x) rng_i - surv_func(x)
    if(root_func(cens) < 0){
      y <- cens
      event <- 0
    } else {
      root <- uniroot(root_func, c(0, cens), tol = 1e-6)
      y <- root$root
      event <- 1
    }

    # format the data
    Z <- matrix(Z, 1)
    colnames(Z) <- paste0("Z", 1:NCOL(Z) - 1L)

    surv_data <- cbind(y = y, event = event, Z[, -1, drop = FALSE], id = id)

    # handle the markers
    # sample the observations times
    obs_time <- cumsum(c(0, rexp(20, .5)))
    obs_time <- obs_time[obs_time < y]
    n_obs <- length(obs_time)

    # sample the fixed effects
    X <- cbind(rep(1, n_obs), X = rnorm(1))
    colnames(X) <- paste0("X", 1:NCOL(X) - 1L)

    # sample the outcomes
    eta <- X %*% fixef_marker +
      g_funcs(obs_time) %*% fixef_vary_marker +
      m_funcs(obs_time) %*% U

    y <- eta + rnorm(n_obs, sd = sqrt(vcov_marker))

      marker_data <- cbind(
        Y = drop(y), X[, -1, drop = FALSE], time = obs_time, id = id)

    list(marker_data = marker_data, surv_data = surv_data)
  })

  # combine the data and return
  marker_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "marker_data")))
  marker_data$id <- as.integer(marker_data$id)
  # the order does not matter
  marker_data <- marker_data[sample.int(NROW(marker_data)), ]

  surv_data <- as.data.frame(do.call(
    rbind, lapply(dat, `[[`, "surv_data")))
  surv_data$id <- as.integer(surv_data$id)
  # the order does not matter
  surv_data <- surv_data[sample.int(NROW(surv_data)), ]

  list(marker_data = marker_data, surv_data = surv_data)
}

# fit the model w/ the VA method
set.seed(1)
dat <- sim_dat(400L)

mean(dat$surv_data$event) # event rate
# quantiles of the uncensored outcomes
subset(dat$surv_data, event == 1, y, TRUE) |>
  quantile(probs = seq(0, 1, length.out = 11))

# distribution in the number of observed markers
quantile(table(dat$marker_data$id),
         probs = seq(0, 1, length.out = 11))
mean(table(dat$marker_data$id))

# add the covariate from the marker to the survival outcome
dat$surv_data <- merge(
  dat$surv_data, subset(dat$marker_data, !duplicated(id), c(id, X1)),
  by = "id")

# fit the model with the variational approximation
n_nodes <- 32L
VA_time <- system.time({
  marker_1 <- marker_term(
      Y ~ X1, id = id, dat$marker_data,
      time_fixef = ns_term(time, knots = c(3, 5), Boundary.knots = c(1, 8)),
      time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(1, 8),
                         intercept = TRUE))

  bks <-  range(dat$surv_data$y)
  iks <- head(quantile(dat$surv_data$y, length.out = 9)[-1], -1)
  surv_obj <- surv_term(
    Surv(y, event) ~ Z1 + X1, id = id, dat$surv_data,
    time_fixef = ns_term(y, Boundary.knots = bks, knots = iks),
    with_frailty = FALSE)

  # use the same number of nodes as JMbayes
  library(SimSurvNMarker)
  gl_rule <- within(get_gl_rule(n_nodes), {
    node <- node/2 + 0.5
    weight <- weight/2
  })

  comp_obj <- joint_ms_ptr(markers = marker_1,
                           survival_terms = surv_obj, max_threads = 4L,
                           quad_rule = gl_rule)

  # get the starting values
  start_val <- joint_ms_start_val(comp_obj)

  # find the maximum lower bound
  opt_out <- joint_ms_opt(comp_obj, par = start_val, max_it = 1000L,
                          pre_method = 3L, cg_tol = .2, c2 = .1,
                          rel_eps = 1e-12)
})
VA_time # estimation time
sqrt(sum(joint_ms_lb_gr(comp_obj, opt_out$par)^2)) # gradient norm

# check the results
opt_out$convergence # did it converge?
opt_out$counts
fmt_par <- joint_ms_format(comp_obj, opt_out$par)

fmt_par$markers
fixef_marker
fixef_vary_marker

fmt_par$survival
c(fixef_surv, 0)
association

fmt_par$vcov$vcov_vary
vcov_vary

fmt_par$vcov$vcov_marker
vcov_marker

# starting values
joint_ms_format(comp_obj, start_val)

# construct an approximate profile likelihood based confidence interval for the
# association parameter
system.time(joint_pl <- joint_ms_profile(
  comp_obj, opt_out,
  which_prof = comp_obj$indices$survival[[1]]$associations, delta = .5,
  rel_eps = 1e-12, cg_tol = .2, c2 = .1))
joint_pl$confs # the confidence interval

# plot the log profile likelihood
with(joint_pl, {
  plot(xs, p_log_Lik, pch = 16, bty = "l",
       xlab = expression(alpha), ylab = "Log profile likelihood")

  smooth_est <- smooth.spline(xs, p_log_Lik)
  lines(predict(smooth_est, seq(min(xs), max(xs), length.out = 100)))
  abline(v = confs, lty = 2)
})

# get the observed information matrix instead
system.time(joint_obs_mat <- joint_ms_hess(comp_obj, opt_out$par))
joint_obs_mat <- -joint_obs_mat$hessian
joint_vcov <- solve(-joint_obs_mat)

# fit the model with JMbayes
library(JMbayes)

JMbayes_time <- system.time({
  # data has to be sorted the same way
  dat <- within(dat, {
    marker_data <- with(marker_data, marker_data[order(id, time), ])
    surv_data <- with(surv_data, surv_data[order(id, y), ])
  })

  # setup initial objects
  lme_fit <- lme(
    Y ~ X1 + ns(time, knots = c(3, 5), Boundary.knots = c(1, 8)),
    random = ~ ns(time, knots = numeric(), Boundary.knots = c(1, 8),
                  intercept = TRUE) - 1 | id,
    data = dat$marker_data,
    control = lmeControl(
      maxIter = 1000L, msMaxIter = 1000L, msMaxEval = 10000L))

  surv_fit <- coxph(Surv(y, event) ~ Z1 + X1, dat$surv_data, x = TRUE)

  # run the MCMC method
  JMbayes_fit <- jointModelBayes(lme_fit, surv_fit, timeVar = "time",
                                 # we use the same number Gauss-Legendre nodes
                                 GQsurv = "GaussLegendre", GQsurv.k = n_nodes,
                                 # avoid the penalization to ensure we should get
                                 # roughly an identical model
                                 baseHaz = "regression-splines",
                                 verbose = FALSE)
})

JMbayes_time # estimation time
JMbayes_time["elapsed"] / VA_time["elapsed"] # relative estimation time

# compare estimates
summary(JMbayes_fit)

association
fixef_surv[2] # Z1

# the starting values
summary(lme_fit)
summary(surv_fit)

plot(JMbayes_fit, ask = FALSE) # any mixing problems?

# fit the model with JM
library(JM)

JM_time <- system.time({
  # data has to be sorted the same way
  dat <- within(dat, {
    marker_data <- with(marker_data, marker_data[order(id, time), ])
    surv_data <- with(surv_data, surv_data[order(id, y), ])
  })

  # setup initial objects
  lme_fit <- lme(
    Y ~ X1 + ns(time, knots = c(3, 5), Boundary.knots = c(1, 8)),
    random = ~ ns(time, knots = numeric(), Boundary.knots = c(1, 8),
                  intercept = TRUE) - 1 | id,
    data = dat$marker_data,
    control = lmeControl(
      maxIter = 1000L, msMaxIter = 1000L, msMaxEval = 10000L))

  surv_fit <- coxph(Surv(y, event) ~ Z1 + X1, dat$surv_data, x = TRUE)

  # the number of quadrature nodes is not that comparable. Perhaps fewer
  # Gauss-Hermite quadrature nodes is enough though?
  JM_fit <- jointModel(lme_fit, surv_fit, timeVar = "time",
                       method = "spline-PH-aGH",
                       typeGH = "adaptive", GHk = 15, GKk = 15)
})

JM_time # the estimation time
JM_time["elapsed"] / VA_time["elapsed"] # relative estimation time

# the estimates
summary(JM_fit)

# the two are not directly comparable because of different baseline hazards
print(JM_fit$logLik, digits = 6) # maximum log likelihood estimate
print(-opt_out$value, digits = 6) # the lower bound on the maximum log likelihood

# fit the model with joineRML
library(joineRML)
set.seed(1)
gc()
joineRML_time <- system.time(
  joineRML_fit <- mjoint(
    formLongFixed =
      Y ~ X1 + ns(time, knots = c(3, 5), Boundary.knots = c(1, 8)),
    formLongRandom = Y ~ ns(time, knots = numeric(), Boundary.knots = c(1, 8),
                            intercept = TRUE) - 1 | id,
    formSurv = Surv(y, event) ~ Z1 + X1,
    data = dat$marker_data,
    survData = dat$surv_data,
    timeVar = "time",
    pfs = FALSE))

joineRML_time # the estimation time
joineRML_time["elapsed"] / VA_time["elapsed"] # relative estimation time

# look at the estimates
summary(joineRML_fit)

# gather all the estimators and show them together. We use the posterior mean
# for the MCMC method
true_vals <- list(
  fixef_marker = fixef_marker,
  fixef_vary_marker = fixef_vary_marker,
  # we use the parameterization from JM and JMbayes
  fixef_surv = c(fixef_surv[2], 0 - association * fixef_marker[2]),
  assoc = association,
  err_std = sqrt(vcov_marker),
  vcov_vary_diag = diag(vcov_vary))

mcmc_est <- with(summary(JMbayes_fit), list(
  fixef_marker = `CoefTable-Long`[1:2, "Value"],
  fixef_vary_marker = `CoefTable-Long`[-(1:2), "Value"],
  fixef_surv = `CoefTable-Event`[c("Z1", "X1"), "Value"],
  assoc = `CoefTable-Event`["Assoct", "Value"],
  err_std = sigma,
  vcov_vary_diag = diag(D)))

jm_est <- with(summary(JM_fit), list(
  fixef_marker = `CoefTable-Long`[1:2, "Value"],
  fixef_vary_marker = `CoefTable-Long`[-(1:2), "Value"],
  fixef_surv = `CoefTable-Event`[c("Z1", "X1"), "Value"],
  assoc = `CoefTable-Event`["Assoct", "Value"],
  err_std = sigma,
  vcov_vary_diag = diag(D)))

va_est <- with(fmt_par, list(
  fixef_marker = markers[[1]]$fixef,
  fixef_vary_marker = markers[[1]]$fixef_vary,
  # we use the parameterization from JM and JMbayes
  fixef_surv = with(
    survival[[1]],
    c(fixef[2], fixef[3] - associations * markers[[1]]$fixef[2])),
  assoc =  survival[[1]]$associations,
  err_std = sqrt(vcov$vcov_marker),
  vcov_vary_diag = diag(vcov$vcov_vary)))

joineRML_est <- with(coef(joineRML_fit), list(
  fixef_marker = beta[1:2],
  fixef_vary_marker = beta[3:5],
  # we use the parameterization from JM and JMbayes
  fixef_surv = c(gamma[1], gamma[2] - gamma[3] * beta[2]),
  assoc = gamma[3],
  err_std = sqrt(sigma2),
  vcov_vary_diag = diag(D)))

# all parameters
rbind(JMbayes = unlist(mcmc_est),
      JM = unlist(jm_est),
      VAJointSurv = unlist(va_est),
      joineRML = unlist(joineRML_est),
      Truth = unlist(true_vals))

# the estimated association parameter with bounds
rbind(
  JMbayes = summary(JMbayes_fit)$`CoefTable-Event`[
    "Assoct", c("Value", "2.5%", "97.5%")],
  JM = {
    ests <- summary(JM_fit)$`CoefTable-Event`["Assoct", ]
    ests["Value"] + c(0, -1.96, 1.96) * ests["Std.Err"]
  },
  VAJointSurv = c(
    fmt_par$survival[[1]]$associations, joint_pl$confs),
  joineRML = c(joineRML_est$assoc, NA, NA),
  truth = c(association, NA, NA))

# compare the covariance matrix estimates
JM_vcov <- vcov(JM_fit)
JM_vcov_keep <- which(!grepl(
  "(^T\\.bs\\d|\\.sigma$|^B\\.D\\d)", rownames(JM_vcov)))
JM_vcov <- JM_vcov[JM_vcov_keep, JM_vcov_keep]

VA_vcov <- joint_vcov
indices <- with(
  comp_obj$indices, c(
    markers[[1]]$fixef, markers[[1]]$fixef_vary,
    survival[[1]]$fixef[-1], # the intercept
    survival[[1]]$associations))
VA_vcov <- VA_vcov[indices, indices]
nams <- unlist(comp_obj$param_names$param_names)
dimnames(VA_vcov) <- list(nams[indices], nams[indices])

# compare the covariance matrix estimates (the X1 should not match because of
# different parameterizations)
JM_vcov
VA_vcov

# account for the different parameterizations (i.e. transform one)
# TODO: very hard coded
R <- diag(NCOL(VA_vcov))
R[NCOL(R) - 1L, 2] <- -fmt_par$survival[[1]]$associations
R[NCOL(R) - 1L, NCOL(R)] <- -fmt_par$markers[[1]]$fixef[2]
VA_vcov <- tcrossprod(R %*% VA_vcov, R)

# comparison is now valid
JM_vcov
VA_vcov

norm(JM_vcov - VA_vcov, "F") / norm(JM_vcov, "F") # relative diff

# compare the standard errors (the X1 should not match because of different
# parameterizations)
rbind(JM = sqrt(diag(JM_vcov)),
      VAJointSurv = sqrt(diag(VA_vcov)))

# estimation times
rbind(JMbayes = JMbayes_time,
      JM = JM_time,
      VAJointSurv = VA_time,
      joineRML = joineRML_time)[, 1:3]

sessionInfo()
