library(VAJointSurv)
load("bench-delayed.RData")

# estimate the model with this package. Get the object we need for the
# optimization while NOT account and accounting for the delayed entry
marker_1 <- marker_term(
  Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
  time_fixef = ns_term(time, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 10),
                     intercept = TRUE))
marker_2 <- marker_term(
  Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
  time_fixef = poly_term(time, degree = 2, raw = TRUE),
  time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

# the right way
library(survival)
surv_terminal <- surv_term(
  Surv(delayed_entry, y, event) ~ Z1, id = id, dat$terminal_outcome,
  time_fixef = bs_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE,
  # some have delayed entry
  delayed = delayed_entry > 0)

surv_obs <- surv_term(
  Surv(lf_trunc, y, event) ~ 1, id = id, dat$obs_process,
  time_fixef = ns_term(y, knots = 5, Boundary.knots = c(0, 10)),
  with_frailty = TRUE)

# the right way
comp_obj <- joint_ms_ptr(
  markers = list(marker_1, marker_2),
  survival_terms = list(surv_terminal, surv_obs),
  max_threads = 4L)

bench::mark(joint_ms_lb(comp_obj, comp_obj$start_val))
bench::mark(joint_ms_lb_gr(comp_obj, comp_obj$start_val))

system.time(opt_out <- joint_ms_opt(
  comp_obj, par = comp_obj$start_val, max_it = 100L, pre_method = 3L,
  cg_tol = .2, c2 = .1, gr_tol = 1))

opt_out$value
head(opt_out$par, 37)
