library(VAJointSurv)
load("bench-all.RData")

marker_1 <- marker_term(
  Y1 ~ X1_1, id = id, subset(dat$marker_data, !is.na(Y1)),
  time_fixef = ns_term(time, knots = c(3.33, 6.67), Boundary.knots = c(0, 10)),
  time_rng = ns_term(time, knots = numeric(), Boundary.knots = c(0, 10),
                     intercept = TRUE))
marker_2 <- marker_term(
  Y2 ~ 1, id = id, subset(dat$marker_data, !is.na(Y2)),
  time_fixef = poly_term(time, degree = 2, raw = TRUE),
  time_rng = poly_term(time, degree = 1, raw = TRUE, intercept = TRUE))

library(survival)
surv_terminal <- surv_term(
  Surv(y, event) ~ Z1, id = id, dat$terminal_outcome,
  time_fixef = bs_term(y, knots = 5, Boundary.knots = c(0, 10)),
  ders = ders[[1]])
surv_obs <- surv_term(
  Surv(lf_trunc, y, event) ~ 1, id = id, dat$obs_process,
  time_fixef = ns_term(y, knots = 5, Boundary.knots = c(0, 10)),
  ders = ders[[2]])

comp_obj <- joint_ms_ptr(markers = list(marker_1, marker_2),
                         survival_terms = list(surv_terminal, surv_obs),
                         max_threads = 4L)

system.time(start_val <- joint_ms_start_val(comp_obj))

library(microbenchmark)
microbenchmark(
  `fn 1` = joint_ms_lb   (comp_obj, start_val, n_threads = 1L),
  `fn 2` = joint_ms_lb   (comp_obj, start_val, n_threads = 2L),
  `fn 4` = joint_ms_lb   (comp_obj, start_val, n_threads = 4L),
  `gr 1` = joint_ms_lb_gr(comp_obj, start_val, n_threads = 1L),
  `gr 2` = joint_ms_lb_gr(comp_obj, start_val, n_threads = 2L),
  `gr 4` = joint_ms_lb_gr(comp_obj, start_val, n_threads = 4L), times = 25)

# without caching
microbenchmark(
  `fn 1` = joint_ms_lb   (comp_obj, start_val, n_threads = 1L, cache_expansions = FALSE),
  `fn 2` = joint_ms_lb   (comp_obj, start_val, n_threads = 2L, cache_expansions = FALSE),
  `fn 4` = joint_ms_lb   (comp_obj, start_val, n_threads = 4L, cache_expansions = FALSE),
  `gr 1` = joint_ms_lb_gr(comp_obj, start_val, n_threads = 1L, cache_expansions = FALSE),
  `gr 2` = joint_ms_lb_gr(comp_obj, start_val, n_threads = 2L, cache_expansions = FALSE),
  `gr 4` = joint_ms_lb_gr(comp_obj, start_val, n_threads = 4L, cache_expansions = FALSE),
  times = 25)
