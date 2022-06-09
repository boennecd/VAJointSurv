#include <testthat.h>
#include "ghq-delayed-entry.h"
#include <vector>

/* The tests are run using this package
 *
 *   https://github.com/boennecd/ghq-cpp/tree/main/ghqCpp
 *
 * The commit is 61f0506c21 as of this writing.
 */

using survival::delayed_dat;

namespace {

/*
 dput(ghq_dat <- fastGHQuad::gaussHermiteData(10))
 gl_dat <- with(SimSurvNMarker::get_gl_rule(30),
 list(node = (node + 1) / 2, weight = weight / 2))
 dput(gl_dat)

 # util functions
 upper_to_full <- \(x){
 dim <- (sqrt(8 * length(x) + 1) - 1) / 2
 out <- matrix(0, dim, dim)
 out[upper.tri(out, TRUE)] <- x
 out[lower.tri(out)] <- t(out)[lower.tri(out)]
 out
 }

 get_n_remove <- \(x, n){
 out <- x[1:n]
 eval(substitute(out <- out[-(1:n)], list(out = substitute(x), n = n)),
 parent.frame())
 out
 }

 d_upper_to_full <- \(x){
 dim <- (sqrt(8 * length(x) + 1) - 1) / 2
 out <- matrix(0, dim, dim)
 out[upper.tri(out, TRUE)] <- x
 out[upper.tri(out)] <- out[upper.tri(out)] / 2
 out[lower.tri(out)] <- t(out)[lower.tri(out)]
 out
 }
 */
constexpr size_t n_ghq{10};
constexpr double ghq_nodes[]{-3.43615911883774, -2.53273167423279, -1.75668364929988, -1.03661082978951, -0.342901327223705, 0.342901327223705, 1.03661082978951, 1.75668364929988, 2.53273167423279, 3.43615911883774},
                 ghq_weights[]{7.6404328552326e-06, 0.00134364574678124, 0.0338743944554811, 0.240138611082314, 0.610862633735326, 0.610862633735326, 0.240138611082315, 0.033874394455481, 0.00134364574678124, 7.64043285523265e-06};

constexpr size_t n_gl{30};
constexpr double gl_nodes[]{0.998446742037325, 0.991834061639874, 0.980010932484154, 0.963100023714637, 0.941280267896026, 0.914782881191384, 0.883888716052413, 0.848925247396658, 0.810263091494621, 0.76831207407101, 0.723516884769045, 0.676352362765439, 0.627318463083945, 0.576934956804292, 0.525735921277659, 0.474264078722341, 0.423065043195708, 0.372681536916055, 0.323647637234561, 0.276483115230955, 0.23168792592899, 0.189736908505379, 0.151074752603342, 0.116111283947587, 0.0852171188086158, 0.0587197321039737, 0.0368999762853629, 0.0199890675158463, 0.00816593836012641, 0.00155325796267525},
              gl_wewights[]{0.0039840962480833, 0.00923323415554548, 0.0143923539416617, 0.0193995962848135, 0.024201336415297, 0.0287465781088095, 0.0329871149410902, 0.0368779873688526, 0.0403779476147101, 0.0434498936005415, 0.0460612611188931, 0.0481843685873221, 0.0497967102933976, 0.0508811948742027, 0.0514263264467794, 0.0514263264467794, 0.0508811948742027, 0.0497967102933976, 0.0481843685873221, 0.0460612611188931, 0.0434498936005415, 0.0403779476147101, 0.0368779873688526, 0.0329871149410902, 0.0287465781088095, 0.024201336415297, 0.0193995962848135, 0.0143923539416617, 0.00923323415554548, 0.0039840962480833};

}

context("delayed_dat functions") {
  test_that("works with two survival outcomes of the same type with fraitly") {
    /*
     library(ghqCpp)

     set.seed(2)
     k1 <- 3
     k2 <- 1
     S1 <- rWishart(1, k1, diag(1/k1, k1)) |> drop() |> round(3)
     S2 <- rWishart(1, k2, diag(1/k2, k2)) |> drop() |> round(3)
     Sigma <- matrix(0, k1 + k2, k1 + k2)
     Sigma[1:k1, 1:k1] <- S1
     Sigma[1:k2 + k1, 1:k2 + k1] <- S2

     dput(S1)
     dput(S2)

     gamma <- c(.15, .33)
     delayed <- c(.5, 1)
     alpha <- c(.67, .25)

     eta <- \(x) cbind(1, log(x)) %*% gamma |> drop()
     M <- \(x) cbind(alpha[1], alpha[1] * x, alpha[2] * x, 1)

     expand <- \(x, fun) fun(x * gl_dat$node)

     M_comb <- lapply(delayed, expand, fun = M) |> do.call(what = rbind)
     eta_comb <- lapply(delayed, expand, fun = eta) |> unlist()
     ws <- c(gl_dat$weight %o% delayed)

     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log() |> dput()

     fn <- \(par){
     gamma <- get_n_remove(par, length(gamma))
     alpha <- get_n_remove(par, length(alpha))
     S1 <- upper_to_full(get_n_remove(par, k1 * (k1 + 1) / 2))
     S2 <- upper_to_full(par)

     Sigma <- matrix(0, k1 + k2, k1 + k2)
     Sigma[1:k1, 1:k1] <- S1
     Sigma[1:k2 + k1, 1:k2 + k1] <- S2

     environment(eta) <- environment()
     environment(M) <- environment()

     M_comb <- lapply(delayed, expand, fun = M) |> do.call(what = rbind)
     eta_comb <- lapply(delayed, expand, fun = eta) |> unlist()
     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log()
     }

     grad <- numDeriv::grad(fn, c(gamma, alpha, S1[upper.tri(S1, TRUE)], S2))

     get_n_remove(grad, length(gamma)) |> dput()
     get_n_remove(grad, length(alpha)) |> dput()
     get_n_remove(grad, k1 * (k1 + 1) / 2) |> d_upper_to_full() |> dput()
     d_upper_to_full(grad) |> dput()
     */

    constexpr double assoc[]{.67, .25},
                  fix_surv[]{.15, .33},
                     delay[]{.5, 1},
                   true_fn{-1.15679913510741},
                 vcov_vary[]{0.203, -0.294, -0.062, -0.294, 0.46, 0.302, -0.062, 0.302, 1.617},
                 vcov_surv[]{0.555},
                   d_assoc[]{-0.0184346377802509, -0.00919147692590657},
                   d_gamma[]{-0.881901733102239, 0.841649437411157},
               d_vcov_vary[]{0.050779395782231, 0.0254164927279688, 0.00948376621481518, 0.0254164927279688, -0.00166047224672699, -0.000619579387180349, 0.00948376621481518, -0.000619579387180349, -0.000231186169974813},
                 d_vcov_surv{0.113124530606547};

    survival::node_weight const gl_dat
      {gl_nodes, gl_wewights, static_cast<vajoint_uint>(n_gl)};
    ghqCpp::ghq_data const ghq_dat{ghq_nodes, ghq_weights, n_ghq};
    ghqCpp::simple_mem_stack<double> mem;

    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;

    bases_fix.emplace_back(new joint_bases::orth_poly(1, false, true));

    bases_rng.emplace_back(new joint_bases::orth_poly(1, true));
    bases_rng.emplace_back(new joint_bases::orth_poly(1, false));

    {
      subset_params params;
      params.add_marker({ 1L, 2L, 2L});
      params.add_marker({ 2L, 2L, 1L});
      params.add_surv({ 1L, 1L, {1, 1}, true});

      std::vector<double> x(params.n_params(), 0);
      x[params.association(0)] = assoc[0];
      x[params.association(0) + 1] = assoc[1];
      x[params.fixef_surv(0)] = fix_surv[0];
      x[params.fixef_vary_surv(0)] = fix_surv[1];
      std::copy(vcov_vary, vcov_vary + 9, x.data() + params.vcov_vary());
      std::copy(vcov_surv, vcov_surv + 1, x.data() + params.vcov_surv());

      double dsgn[]{1, 1};
      std::vector<simple_mat<double> > design_mats{{dsgn, 1, 2}},
                              fixef_design_varying{{nullptr, 0, 2}},
                              rng_design_varying{{nullptr, 0, 2}};

      std::vector<std::vector<std::vector<int> > > ders{{{0}, {0}}};
      std::vector<delayed_dat::cluster_info>
        info{{{0, 0, delay[0]}, {0, 1, delay[1]}}};

      survival::delayed_dat cmp_dat
        {bases_fix, bases_rng, design_mats, fixef_design_varying,
         rng_design_varying, params, info, ders};

      {
        double const res{cmp_dat(x.data(), mem, 0, gl_dat, ghq_dat)};
        expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
      }

      constexpr double dum_off{-1.5};
      std::vector<double> gr(params.n_params(), dum_off);

      mem.reset();
      double const res
        {cmp_dat.grad(x.data(), gr.data(), mem, 0, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

      expect_true
        (std::abs(gr[params.fixef_surv(0)] - d_gamma[0] - dum_off)
           < std::abs(d_gamma[0]) * 1e-3);
      expect_true
        (std::abs(gr[params.fixef_vary_surv(0)] - d_gamma[1] - dum_off)
           < std::abs(d_gamma[0]) * 1e-3);

      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[params.association(0) + i] - d_assoc[i] - dum_off)
             < std::abs(d_assoc[i]) * 1e-3);

      for(size_t i = 0; i < 9; ++i)
        expect_true
          (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i] - dum_off)
             < std::abs(d_vcov_vary[i]) * 1e-3);

      expect_true
        (std::abs(gr[params.vcov_surv()] - d_vcov_surv - dum_off)
           < std::abs(d_vcov_surv) * 1e-3);
    }

    /*  works with another redundant survival type with frailty and another
     *  observation of the same type. */
    mem.reset();

    bases_fix.clear();

    bases_fix.emplace_back(new joint_bases::orth_poly(2, true));
    bases_fix.emplace_back(new joint_bases::orth_poly(1, false, true));

    subset_params params;
    params.add_marker({ 1L, 2L, 2L});
    params.add_marker({ 2L, 2L, 1L});
    params.add_surv({ 2L, 3L, {1, 1}, true});
    params.add_surv({ 1L, 1L, {1, 1}, true});

    std::vector<double> x(params.n_params(), 0);
    x[params.association(1)] = assoc[0];
    x[params.association(1) + 1] = assoc[1];
    x[params.fixef_surv(1)] = fix_surv[0];
    x[params.fixef_vary_surv(1)] = fix_surv[1];
    std::copy(vcov_vary, vcov_vary + 9, x.data() + params.vcov_vary());
    x[params.vcov_surv()] = -1;
    x[params.vcov_surv() + 1] = -1;
    x[params.vcov_surv() + 2] = -1;
    x[params.vcov_surv() + 3] = *vcov_surv;

    double dsgn_dum[]{0, 1, 2, 3}, dsgn[]{0, 1, 1};
    std::vector<simple_mat<double> >
      design_mats{{dsgn_dum, 2, 1}, {dsgn, 1, 3}},
      fixef_design_varying{{nullptr, 0, 1}, {nullptr, 0, 3}},
      rng_design_varying{{nullptr, 0, 1}, {nullptr, 0, 3}};

    std::vector<std::vector<std::vector<int> > > ders{{{1}, {-1}}, {{0}, {0}}};
    std::vector<delayed_dat::cluster_info>
      info{{{0, 0, 2}, {1, 0, 4}},
           {{1, 1, delay[0]}, {1, 2, delay[1]}}};

    survival::delayed_dat cmp_dat
      {bases_fix, bases_rng, design_mats, fixef_design_varying,
       rng_design_varying, params, info, ders};

    for(size_t n_evals = 0; n_evals < 4; ++n_evals){
      if(n_evals > 0)
        cmp_dat.set_cached_expansions(gl_dat, mem);

      {
        mem.reset();
        double const res{cmp_dat(x.data(), mem, 1, gl_dat, ghq_dat)};
        expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
      }

      constexpr double dum_off{-1.5};
      std::vector<double> gr(params.n_params(), dum_off);

      mem.reset();
      double const res
        {cmp_dat.grad(x.data(), gr.data(), mem, 1, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

      expect_true
        (std::abs(gr[params.fixef_surv(1)] - d_gamma[0] - dum_off)
           < std::abs(d_gamma[0]) * 1e-3);
      expect_true
        (std::abs(gr[params.fixef_vary_surv(1)] - d_gamma[1] - dum_off)
           < std::abs(d_gamma[0]) * 1e-3);

      for(size_t i = 0; i < 2; ++i)
        expect_true
        (std::abs(gr[params.association(1) + i] - d_assoc[i] - dum_off)
           < std::abs(d_assoc[i]) * 1e-3);

      for(size_t i = 0; i < 9; ++i)
        expect_true
        (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i] - dum_off)
           < std::abs(d_vcov_vary[i]) * 1e-3);

      expect_true
        (std::abs(gr[params.vcov_surv() + 3] - d_vcov_surv - dum_off)
           < std::abs(d_vcov_surv) * 1e-3);
    }
  }

  test_that("works with two survival outcomes of the same type without fraitly and with current value and slope") {
    /*
     library(ghqCpp)

     set.seed(3)
     dput(Sigma <- rWishart(1, 3, diag(1/3, 3)) |> drop() |> round(3))

     gamma <- c(.35, .5)
     delayed <- c(.66, 1.25)
     alpha <- c(.67, 1, .25)

     eta <- \(x) cbind(1, log(x)) %*% gamma |> drop()
     M <- \(x) cbind(alpha[1] * x + alpha[2], alpha[1] * x^2 + 2 * alpha[2] * x, alpha[3])

     expand <- \(x, fun) fun(x * gl_dat$node)

     M_comb <- lapply(delayed, expand, fun = M) |> do.call(what = rbind)
     eta_comb <- lapply(delayed, expand, fun = eta) |> unlist()
     ws <- c(gl_dat$weight %o% delayed)

     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log() |> dput()

     fn <- \(par){
     gamma <- get_n_remove(par, length(gamma))
     alpha <- get_n_remove(par, length(alpha))
     Sigma <- upper_to_full(par)

     environment(eta) <- environment()
     environment(M) <- environment()

     M_comb <- lapply(delayed, expand, fun = M) |> do.call(what = rbind)
     eta_comb <- lapply(delayed, expand, fun = eta) |> unlist()
     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log()
     }

     grad <- numDeriv::grad(fn, c(gamma, alpha, Sigma[upper.tri(Sigma, TRUE)]))

     get_n_remove(grad, length(gamma)) |> dput()
     get_n_remove(grad, length(alpha)) |> dput()
     d_upper_to_full(grad) |> dput()
     */

    constexpr double assoc[]{.67, 1, .25},
                  fix_surv[]{.35, .5},
                     delay[]{.66, 1.25},
                     true_fn{-1.53142252348068},
                 vcov_vary[]{0.18, -0.282, -0.298, -0.282, 0.599, 0.757, -0.298, 0.757, 1.806},
                   d_gamma[]{-0.99038243533505, 0.817424734218724},
                   d_assoc[]{0.0315009345747261, 0.147744995604926, 0.605116537123115},
               d_vcov_vary[]{0.435080048959208, 0.512740604206576, 0.0753189883308562, 0.512740604206576, 0.404105825746256, 0.0999765811449094, 0.0753189883308562, 0.0999765811449094, 0.012404404002973};

    survival::node_weight const gl_dat
      {gl_nodes, gl_wewights, static_cast<vajoint_uint>(n_gl)};
    ghqCpp::ghq_data const ghq_dat{ghq_nodes, ghq_weights, n_ghq};
    ghqCpp::simple_mem_stack<double> mem;

    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;

    bases_fix.emplace_back(new joint_bases::orth_poly(1, true, true));

    bases_rng.emplace_back(new joint_bases::orth_poly(2, false));
    bases_rng.emplace_back(new joint_bases::orth_poly(1, false));

    subset_params params;
    params.add_marker({ 1L, 2L, 2L});
    params.add_marker({ 2L, 2L, 1L});
    params.add_surv({ 0L, 2L, {2, 1}, false});

    {
      std::vector<double> x(params.n_params(), 0);
      std::copy(assoc, assoc + 3, x.data() + params.association(0));

      std::copy(fix_surv, fix_surv + 2, x.data() + params.fixef_vary_surv(0));

      std::copy(vcov_vary, vcov_vary + 9, x.data() + params.vcov_vary());

      std::vector<simple_mat<double> >
        design_mats{{nullptr, 0, 2}},
        fixef_design_varying{{nullptr, 0, 2}},
        rng_design_varying{{nullptr, 0, 2}};

      std::vector<std::vector<std::vector<int> > > ders{{{0, 1}, {1}}};
      std::vector<delayed_dat::cluster_info>
        info{{{0, 0, delay[0]}, {0, 1, delay[1]}}};

      survival::delayed_dat cmp_dat
        {bases_fix, bases_rng, design_mats, fixef_design_varying,
         rng_design_varying, params, info, ders};

      {
        double const res{cmp_dat(x.data(), mem, 0, gl_dat, ghq_dat)};
        expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
      }

      std::vector<double> gr(params.n_params(), 0);
      mem.reset();
      double const res
        {cmp_dat.grad(x.data(), gr.data(), mem, 0, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[params.fixef_vary_surv(0) + i] - d_gamma[i]) <
            std::abs(d_gamma[i]) * 1e-3);
      for(size_t i = 0; i < 3; ++i)
        expect_true
          (std::abs(gr[params.association(0) + i] - d_assoc[i]) <
            std::abs(d_assoc[i]) * 1e-3);
      for(size_t i = 0; i < 9; ++i)
        expect_true
        (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i]) <
          std::abs(d_vcov_vary[i]) * 1e-3);
    }

    // works with a another redundant survival type without frailty
    bases_fix.emplace_back(new joint_bases::orth_poly(3, true, true));

    params.add_surv({ 3L, 4L, {2, 2}, false});

    std::vector<double> x(params.n_params(), 0);
    std::copy(assoc, assoc + 3, x.data() + params.association(0));

    std::copy(fix_surv, fix_surv + 2, x.data() + params.fixef_vary_surv(0));

    std::copy(vcov_vary, vcov_vary + 9, x.data() + params.vcov_vary());

    double dum[]{0,1,3};
    std::vector<simple_mat<double> >
      design_mats{{nullptr, 0, 2}, {dum, 3, 1}},
      fixef_design_varying{{nullptr, 0, 2}, {nullptr, 0, 1}},
      rng_design_varying{{nullptr, 0, 2}, {nullptr, 0, 1}};

    std::vector<std::vector<std::vector<int> > > ders{{{0, 1}, {1}},
                                                      {{-1, 1}, {0, 2}}};
    std::vector<delayed_dat::cluster_info>
      info{{{0, 0, delay[0]}, {0, 1, delay[1]}},
           {{1, 0, 4}}};

    survival::delayed_dat cmp_dat
      {bases_fix, bases_rng, design_mats, fixef_design_varying,
       rng_design_varying, params, info, ders};

    for(size_t n_evals = 0; n_evals < 4; ++n_evals){
      if(n_evals > 0)
        cmp_dat.set_cached_expansions(gl_dat, mem);

      {
        mem.reset();
        double const res{cmp_dat(x.data(), mem, 0, gl_dat, ghq_dat)};
        expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
      }

      std::vector<double> gr(params.n_params(), 0);
      mem.reset();
      double const res
        {cmp_dat.grad(x.data(), gr.data(), mem, 0, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[params.fixef_vary_surv(0) + i] - d_gamma[i]) <
            std::abs(d_gamma[i]) * 1e-3);
      for(size_t i = 0; i < 3; ++i)
        expect_true
          (std::abs(gr[params.association(0) + i] - d_assoc[i]) <
            std::abs(d_assoc[i]) * 1e-3);

      for(size_t i = 0; i < 9; ++i)
        expect_true
          (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i]) <
            std::abs(d_vcov_vary[i]) * 1e-3);
    }
  }

  test_that("works with two survival outcomes of different types with fraitly") {
    /*
     library(ghqCpp)

     set.seed(2)
     k1 <- 3
     k2 <- 2
     S1 <- rWishart(1, k1, diag(1/k1, k1)) |> drop() |> round(3)
     S2 <- rWishart(1, k2, diag(1/k2, k2)) |> drop() |> round(3)
     Sigma <- matrix(0, k1 + k2, k1 + k2)
     Sigma[1:k1, 1:k1] <- S1
     Sigma[1:k2 + k1, 1:k2 + k1] <- S2

     dput(S1)
     dput(S2)

     gamma <- list(c(.15, .33), c(-.5, .75))
     delayed <- c(.5, 1)
     alpha <- list(c(.67, .25), c(-.2, .4))

     eta <- \(x, type) cbind(1, log(x)) %*% gamma[[type]] |> drop()
     M <- \(x, type){
     alpha <- alpha[[type]]
     cbind(alpha[1], alpha[1] * x, alpha[2] * x, type == 1, type == 2)
     }

     expand <- \(x, ..., fun) fun(x * gl_dat$node, ...)

     type <- c(2:1)
     M_comb <- mapply(expand, delayed, type, MoreArgs = list(fun = M),
     SIMPLIFY = FALSE) |> do.call(what = rbind)
     eta_comb <- mapply(expand, delayed, type, MoreArgs = list(fun = eta),
     SIMPLIFY = FALSE) |> unlist()
     ws <- c(gl_dat$weight %o% delayed)

     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log() |> dput()

     fn <- \(par){
     gamma <- get_n_remove(par, lengths(gamma) |> sum())
     gamma <- list(gamma[1:2], gamma[3:4])
     alpha <- get_n_remove(par, lengths(alpha) |> sum())
     alpha <- list(alpha[1:2], alpha[3:4])
     S1 <- upper_to_full(get_n_remove(par, k1 * (k1 + 1) / 2))
     S2 <- upper_to_full(par)

     Sigma <- matrix(0, k1 + k2, k1 + k2)
     Sigma[1:k1, 1:k1] <- S1
     Sigma[1:k2 + k1, 1:k2 + k1] <- S2

     environment(eta) <- environment()
     environment(M) <- environment()

     M_comb <- mapply(expand, delayed, type, MoreArgs = list(fun = M),
     SIMPLIFY = FALSE) |> do.call(what = rbind)
     eta_comb <- mapply(expand, delayed, type, MoreArgs = list(fun = eta),
     SIMPLIFY = FALSE) |> unlist()
     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log()
     }

     grad <- numDeriv::grad(fn, c(
     unlist(gamma), unlist(alpha), S1[upper.tri(S1, TRUE)],
     S2[upper.tri(S2, TRUE)]))

     get_n_remove(grad, lengths(gamma) |> sum()) |> dput()
     get_n_remove(grad, lengths(alpha) |> sum()) |> dput()
     get_n_remove(grad, k1 * (k1 + 1) / 2) |> d_upper_to_full() |> dput()
     d_upper_to_full(grad) |> dput()
     */

    constexpr double assoc1[]{.67, .25},
                     assoc2[]{-.2, .4},
                  fix_surv1[]{.15, .33},
                  fix_surv2[]{-.5, .75},
                      delay[]{.5, 1},
                      true_fn{-1.075995739615},
                  vcov_vary[]{0.203, -0.294, -0.062, -0.294, 0.46, 0.302, -0.062, 0.302, 1.617},
                  vcov_surv[]{0.407, -0.466, -0.466, 0.745},
                    d_gamma[]{-0.778544337094796, 0.592757753884664, -0.17034927450825, 0.215809741944685},
                    d_assoc[]{-0.0193417235318428, -0.0124083793784138, 0.00304473971423309, 0.000829696718818985},
                d_vcov_vary[]{-0.00561347289321644, -0.000138079225401536, 0.00682683520257915, -0.000138079225401536, -0.0111717916150825, -0.000540786814489186, 0.00682683520257915, -0.000540786814489186, -0.000493831998060587},
                d_vcov_surv[]{0.0227134396032736, 0.0501122918817207, 0.0501122918817207, -0.0594930011943155};

    survival::node_weight const gl_dat
    {gl_nodes, gl_wewights, static_cast<vajoint_uint>(n_gl)};
    ghqCpp::ghq_data const ghq_dat{ghq_nodes, ghq_weights, n_ghq};
    ghqCpp::simple_mem_stack<double> mem;

    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;

    bases_fix.emplace_back(new joint_bases::orth_poly(1, false, true));
    bases_fix.emplace_back(new joint_bases::orth_poly(1, true, true));

    bases_rng.emplace_back(new joint_bases::orth_poly(1, true));
    bases_rng.emplace_back(new joint_bases::orth_poly(1, false));

    {
      subset_params params;
      params.add_marker({ 1L, 2L, 2L});
      params.add_marker({ 2L, 2L, 1L});
      params.add_surv({ 1L, 1L, {1, 1}, true});
      params.add_surv({ 0L, 2L, {1, 1}, true});

      std::vector<double> x(params.n_params(), 0);
      std::copy(assoc1, assoc1 + 2, x.data() + params.association(0));
      std::copy(assoc2, assoc2 + 2, x.data() + params.association(1));

      x[params.fixef_surv(0)] = fix_surv1[0];
      x[params.fixef_vary_surv(0)] = fix_surv1[1];

      std::copy(fix_surv2, fix_surv2 + 2, x.data() + params.fixef_vary_surv(1));

      std::copy(vcov_vary, vcov_vary + 9, x.data() + params.vcov_vary());
      std::copy(vcov_surv, vcov_surv + 4, x.data() + params.vcov_surv());

      double dsgn1[]{1};
      std::vector<simple_mat<double> >
        design_mats{{dsgn1, 1, 1}, {nullptr, 0, 1}},
        fixef_design_varying{{nullptr, 0, 1}, {nullptr, 0, 1}},
        rng_design_varying{{nullptr, 0, 1}, {nullptr, 0, 1}};

      std::vector<std::vector<std::vector<int> > > ders
        {{{0}, {0}}, {{0}, {0}}};
      std::vector<delayed_dat::cluster_info>
        info{{{1, 0, delay[0]}, {0, 0, delay[1]}}};

      survival::delayed_dat cmp_dat
        {bases_fix, bases_rng, design_mats, fixef_design_varying,
         rng_design_varying, params, info, ders};

      {
        double const res{cmp_dat(x.data(), mem, 0, gl_dat, ghq_dat)};
        expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
      }

      std::vector<double> gr(params.n_params(), 0);
      mem.reset();
      double const res
        {cmp_dat.grad(x.data(), gr.data(), mem, 0, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[i + params.association(0)] - d_assoc[i])
             < std::abs(d_assoc[i]) * 1e-3);
      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[i + params.association(1)] - d_assoc[i + 2])
             < std::abs(d_assoc[i + 2]) * 1e-3);

      expect_true
        (std::abs(gr[params.fixef_surv(0)] - d_gamma[0])
          < std::abs(d_gamma[0]) * 1e-3);
      expect_true
        (std::abs(gr[params.fixef_vary_surv(0)] - d_gamma[1])
           < std::abs(d_gamma[1]) * 1e-3);

      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[params.fixef_vary_surv(1) + i] - d_gamma[i + 2])
             < std::abs(d_gamma[i + 2]) * 1e-3);
      for(size_t i = 0; i < 9; ++i)
        expect_true
          (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i])
             < std::abs(d_vcov_vary[i]) * 1e-3);
      for(size_t i = 0; i < 4; ++i)
        expect_true
          (std::abs(gr[params.vcov_surv() + i] - d_vcov_surv[i])
             < std::abs(d_vcov_surv[i]) * 1e-3);
    }

    // works with another redundant survival type with frailty and another
    // survival outcome of one of the types used
    bases_fix.clear();
    bases_fix.emplace_back(new joint_bases::orth_poly(3, true, true));
    bases_fix.emplace_back(new joint_bases::orth_poly(1, false, true));
    bases_fix.emplace_back(new joint_bases::orth_poly(1, true, true));

    subset_params params;
    params.add_marker({ 1L, 2L, 2L});
    params.add_marker({ 2L, 2L, 1L});
    params.add_surv({ 2L, 4L, {2, 1}, true});
    params.add_surv({ 1L, 1L, {1, 1}, true});
    params.add_surv({ 0L, 2L, {1, 1}, true});

    std::vector<double> x(params.n_params(), 0);
    std::copy(assoc1, assoc1 + 2, x.data() + params.association(1));
    std::copy(assoc2, assoc2 + 2, x.data() + params.association(2));

    x[params.fixef_surv(1)] = fix_surv1[0];
    x[params.fixef_vary_surv(1)] = fix_surv1[1];

    std::copy(fix_surv2, fix_surv2 + 2, x.data() + params.fixef_vary_surv(2));

    std::copy(vcov_vary, vcov_vary + 9, x.data() + params.vcov_vary());
    for(unsigned i = 0; i < 2; ++i)
      std::copy(vcov_surv + 2 * i, vcov_surv + 2 * (i + 1),
                x.data() + params.vcov_surv() + 1 + (i + 1) * 3);

    double dsgn0[]{0, 1}, dsgn1[]{1};
    std::vector<simple_mat<double> >
      design_mats{{dsgn0, 2, 1}, {dsgn1, 1, 1}, {nullptr, 0, 2}},
      fixef_design_varying{{nullptr, 0, 1}, {nullptr, 0, 1}, {nullptr, 0, 2}},
      rng_design_varying{{nullptr, 0, 1}, {nullptr, 0, 1}, {nullptr, 0, 2}};

    std::vector<std::vector<std::vector<int> > > ders
      {{{0, 1}, {0}}, {{0}, {0}}, {{0}, {0}}};
    std::vector<delayed_dat::cluster_info>
      info{{{0, 0, 2}, {2, 0, 3}},
           {{2, 1, delay[0]}, {1, 0, delay[1]}}};

    survival::delayed_dat cmp_dat
      {bases_fix, bases_rng, design_mats, fixef_design_varying,
       rng_design_varying, params, info, ders};

    {
      mem.reset();
      double const res{cmp_dat(x.data(), mem, 1, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
    }

    std::vector<double> gr(params.n_params(), 0);
    mem.reset();
    double const res
      {cmp_dat.grad(x.data(), gr.data(), mem, 1, gl_dat, ghq_dat)};
    expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

    for(size_t i = 0; i < 2; ++i)
      expect_true
      (std::abs(gr[i + params.association(1)] - d_assoc[i])
         < std::abs(d_assoc[i]) * 1e-3);
    for(size_t i = 0; i < 2; ++i)
      expect_true
      (std::abs(gr[i + params.association(2)] - d_assoc[i + 2])
         < std::abs(d_assoc[i + 2]) * 1e-3);

    expect_true
      (std::abs(gr[params.fixef_surv(1)] - d_gamma[0])
         < std::abs(d_gamma[0]) * 1e-3);
    expect_true
      (std::abs(gr[params.fixef_vary_surv(1)] - d_gamma[1])
         < std::abs(d_gamma[1]) * 1e-3);

    for(size_t i = 0; i < 2; ++i)
      expect_true
      (std::abs(gr[params.fixef_vary_surv(2) + i] - d_gamma[i + 2])
         < std::abs(d_gamma[i + 2]) * 1e-3);
    for(size_t i = 0; i < 9; ++i)
      expect_true
      (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i])
         < std::abs(d_vcov_vary[i]) * 1e-3);
    for(size_t i = 0; i < 2; ++i)
      expect_true
      (std::abs(gr[params.vcov_surv() + i + 4] - d_vcov_surv[i])
         < std::abs(d_vcov_surv[i]) * 1e-3);
    for(size_t i = 0; i < 2; ++i)
      expect_true
      (std::abs(gr[params.vcov_surv() + i + 7] - d_vcov_surv[i + 2])
         < std::abs(d_vcov_surv[i + 2]) * 1e-3);
  }

  test_that("works with two survival outcomes of different types with fraitly and time-varying effects") {
    /*
     library(ghqCpp)

     set.seed(2)
     k1 <- 3
     k2 <- 2
     S1 <- rWishart(1, k1, diag(1/k1, k1)) |> drop() |> round(3)
     S2 <- rWishart(1, k2, diag(1/k2, k2)) |> drop() |> round(3)
     Sigma <- matrix(0, k1 + k2, k1 + k2)
     Sigma[1:k1, 1:k1] <- S1
     Sigma[1:k2 + k1, 1:k2 + k1] <- S2

     dput(S1)
     dput(S2)

     gamma <- list(c(.15, .33, .1), c(-.5, .75, -.1))
     delayed <- c(.5, 1)
     alpha <- list(c(.67, .25), c(-.2, .4))

     eta <- \(x, type, data)
     cbind(1, log(x), log(x) * data$w1) %*% gamma[[type]] |> drop()
     M <- \(x, type, data){
     alpha <- alpha[[type]]
     cbind(alpha[1] * data$w1, alpha[1] * x * data$w1,
     alpha[2] * x * data$w2, type == 1, type == 2)
     }

     expand <- \(x, ..., fun) fun(x * gl_dat$node, ...)

     dats <- list(data.frame(w1 = .5, w2 = -.25), data.frame(w1 = -.2, w2 = .3))
     sapply(dats, \(x) as.matrix(x) |> t() |> dput()) |> invisible()

     type <- c(2:1)
     M_comb <- mapply(
     expand, delayed, data = dats, type, MoreArgs = list(fun = M),
     SIMPLIFY = FALSE) |>
     do.call(what = rbind)
     eta_comb <- mapply(
     expand, delayed, data = dats, type, MoreArgs = list(fun = eta),
     SIMPLIFY = FALSE) |>
     unlist()
     ws <- c(gl_dat$weight %o% delayed)

     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log() |> dput()

     fn <- \(par){
     gamma <- get_n_remove(par, lengths(gamma) |> sum())
     gamma <- list(gamma[1:3], gamma[4:6])
     alpha <- get_n_remove(par, lengths(alpha) |> sum())
     alpha <- list(alpha[1:2], alpha[3:4])
     S1 <- upper_to_full(get_n_remove(par, k1 * (k1 + 1) / 2))
     S2 <- upper_to_full(par)

     Sigma <- matrix(0, k1 + k2, k1 + k2)
     Sigma[1:k1, 1:k1] <- S1
     Sigma[1:k2 + k1, 1:k2 + k1] <- S2

     environment(eta) <- environment()
     environment(M) <- environment()

     M_comb <- mapply(
     expand, delayed, data = dats, type, MoreArgs = list(fun = M),
     SIMPLIFY = FALSE) |>
     do.call(what = rbind)
     eta_comb <- mapply(
     expand, delayed, data = dats, type, MoreArgs = list(fun = eta),
     SIMPLIFY = FALSE) |>
     unlist()
     expected_survival_term(
     eta = eta_comb, ws = ws, M = M_comb, Sigma = Sigma,
     weights = ghq_dat$w, nodes = ghq_dat$x) |> log()
     }

     grad <- numDeriv::grad(fn, c(
     unlist(gamma), unlist(alpha), S1[upper.tri(S1, TRUE)],
     S2[upper.tri(S2, TRUE)]))

     get_n_remove(grad, lengths(gamma) |> sum()) |> dput()
     get_n_remove(grad, lengths(alpha) |> sum()) |> dput()
     get_n_remove(grad, k1 * (k1 + 1) / 2) |> d_upper_to_full() |> dput()
     d_upper_to_full(grad) |> dput()
     */

    constexpr double assoc1[]{.67, .25},
                     assoc2[]{-.2, .4},
                  fix_surv1[]{.15, .33, .1},
                  fix_surv2[]{-.5, .75, -.1},
                      delay[]{.5, 1},
                      true_fn{-1.03563541595015},
                  vcov_vary[]{0.203, -0.294, -0.062, -0.294, 0.46, 0.302, -0.062, 0.302, 1.617},
                  vcov_surv[]{0.407, -0.466, -0.466, 0.745},
                    d_gamma[]{-0.76951753741138, 0.558143268669321, 0.279071636795202, -0.167967812073718, 0.211258091891914, -0.0422516220849665},
                    d_assoc[]{-0.00348213568965474, -0.000824202026114233, -7.83108030313273e-06, -0.00106970741301091},
                d_vcov_vary[]{0.00186670151030051, 0.000930969784939519, 0.000369758779798111, 0.000930969784939519, -0.00266338529207848, 0.000830529971156172, 0.000369758779798111, 0.000830529971156172, -0.00033695869555832},
                d_vcov_surv[]{0.00627134641737212, 0.0469434055314739, 0.0469434055314739, -0.0596102839820149};

    constexpr size_t n_fix[]{2, 3},
                     n_rng[]{2, 1},
                     vcov_vary_dim{n_rng[0] + n_rng[1]},
                     vcov_surv_dim{2};

    survival::node_weight const gl_dat
    {gl_nodes, gl_wewights, static_cast<vajoint_uint>(n_gl)};
    ghqCpp::ghq_data const ghq_dat{ghq_nodes, ghq_weights, n_ghq};
    ghqCpp::simple_mem_stack<double> mem;

    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;

    {
      joint_bases::bases_vector expansions;
      expansions.emplace_back(new joint_bases::orth_poly(1, false, true));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly>
           (1, false, true));
      bases_fix.emplace_back(new joint_bases::stacked_basis(expansions));

      expansions.clear();
      expansions.emplace_back(new joint_bases::orth_poly(1, true, true));
      expansions.emplace_back
        (new joint_bases::weighted_basis<joint_bases::orth_poly>
           (1, false, true));
      bases_fix.emplace_back(new joint_bases::stacked_basis(expansions));
    }

    bases_rng.emplace_back
      (new joint_bases::weighted_basis<joint_bases::orth_poly>(1, true));
    bases_rng.emplace_back
      (new joint_bases::weighted_basis<joint_bases::orth_poly>(1, false));

    {
      subset_params params;
      params.add_marker({ 1L, 2L, bases_rng[0]->n_basis()});
      params.add_marker({ 2L, 2L, bases_rng[1]->n_basis()});
      params.add_surv({ 1L, n_fix[0], {1, 1}, true});
      params.add_surv({ 0L, n_fix[1], {1, 1}, true});

      std::vector<double> x(params.n_params(), 0);
      std::copy(assoc1, assoc1 + 2, x.data() + params.association(0));
      std::copy(assoc2, assoc2 + 2, x.data() + params.association(1));

      x[params.fixef_surv(0)] = fix_surv1[0];
      std::copy
        (fix_surv1 + 1, fix_surv1 + 1 + n_fix[0],
         x.data() + params.fixef_vary_surv(0));
      std::copy
        (fix_surv2, fix_surv2 + n_fix[1], x.data() + params.fixef_vary_surv(1));

      std::copy(vcov_vary, vcov_vary + vcov_vary_dim * vcov_vary_dim,
                x.data() + params.vcov_vary());
      std::copy(vcov_surv, vcov_surv + vcov_surv_dim * vcov_surv_dim,
                x.data() + params.vcov_surv());

      double dsgn1[]{1},
             dsgn_fix_varying1[]{.5},
             dsgn_fix_varying2[]{-.2},
             dsgn_rng_varying1[]{0.5, -0.25},
             dsgn_rng_varying2[]{-0.2, 0.3};
      std::vector<simple_mat<double> >
        design_mats{{dsgn1, 1, 1}, {nullptr, 0, 1}},
        fixef_design_varying{{dsgn_fix_varying1, 1, 1},
                             {dsgn_fix_varying2, 1, 1}},
        rng_design_varying{{dsgn_rng_varying1, 2, 1},
                           {dsgn_rng_varying2, 2, 1}};

      std::vector<std::vector<std::vector<int> > > ders
        {{{0}, {0}}, {{0}, {0}}};
      std::vector<delayed_dat::cluster_info>
        info{{{1, 0, delay[0]}, {0, 0, delay[1]}}};

      survival::delayed_dat cmp_dat
        {bases_fix, bases_rng, design_mats, fixef_design_varying,
         rng_design_varying, params, info, ders};

      {
        double const res{cmp_dat(x.data(), mem, 0, gl_dat, ghq_dat)};
        expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);
      }

      std::vector<double> gr(params.n_params(), 0);
      mem.reset();
      double const res
        {cmp_dat.grad(x.data(), gr.data(), mem, 0, gl_dat, ghq_dat)};
      expect_true(std::abs(res - true_fn) < std::abs(true_fn) * 1e-6);

      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[i + params.association(0)] - d_assoc[i])
             < std::abs(d_assoc[i]) * 1e-3);
      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[i + params.association(1)] - d_assoc[i + 2])
             < std::abs(d_assoc[i + 2]) * 1e-3);

      expect_true
        (std::abs(gr[params.fixef_surv(0)] - d_gamma[0])
          < std::abs(d_gamma[0]) * 1e-3);
      for(size_t i = 0; i < 2; ++i)
        expect_true
          (std::abs(gr[params.fixef_vary_surv(0) + i] - d_gamma[1 + i])
             < std::abs(d_gamma[1 + i]) * 1e-3);

      for(size_t i = 0; i < 3; ++i)
        expect_true
          (std::abs(gr[params.fixef_vary_surv(1) + i] - d_gamma[i + 3])
             < std::abs(d_gamma[i + 3]) * 1e-3);
      for(size_t i = 0; i < 9; ++i)
        expect_true
          (std::abs(gr[params.vcov_vary() + i] - d_vcov_vary[i])
             < std::abs(d_vcov_vary[i]) * 1e-3);
      for(size_t i = 0; i < 4; ++i)
        expect_true
          (std::abs(gr[params.vcov_surv() + i] - d_vcov_surv[i])
             < std::abs(d_vcov_surv[i]) * 1e-3);
    }
  }
}
