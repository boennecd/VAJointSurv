#include "ghq.h"
#include "ghq-lp-utils.h"

namespace ghqCpp {

adaptive_problem::mode_problem::mode_problem
  (ghq_problem const &problem, simple_mem_stack<double> &mem):
  problem{problem}, mem{mem} { }

double adaptive_problem::mode_problem::func(double const *val){
  double out{};
  for(PSQN::psqn_uint i = 0; i < size(); ++i)
    out += val[i] * val[i];
  out /= 2;
  out -= problem.log_integrand(val, mem);
  return out;
}

double adaptive_problem::mode_problem::grad
  (double const * __restrict__ val, double * __restrict__ gr){
  double const out{-problem.log_integrand_grad(val, gr, mem)};
  std::for_each(gr, gr + size(), [](double &res){ res *= -1; });

  double extra_term{};
  for(PSQN::psqn_uint i = 0; i < size(); ++i){
    extra_term += val[i] * val[i];
    gr[i] += val[i];
  }
  extra_term /= 2;

  return out + extra_term;
}

adaptive_problem::adaptive_problem
  (ghq_problem const &problem, simple_mem_stack<double> &mem,
   double const rel_eps, PSQN::psqn_uint const max_it, double const c1,
   double const c2, double const gr_tol):
  problem{problem} {
    // attempt to find the mode
    mode_problem my_mode_problem(problem, mem);
    mu.zeros(n_vars());

    // TODO: I can avoid the allocation in PSQN::bfgs with minor changes in the
    //       package
    auto res = PSQN::bfgs
      (my_mode_problem, mu.memptr(), rel_eps, max_it, c1, c2, 0L, gr_tol);

    bool succeeded = res.info == PSQN::info_code::converged;
    if(succeeded){
      // we compute the Hessian
      arma::mat hess(mem.get(2 * n_vars() * n_vars()), n_vars(), n_vars(),
                     false),
            hess_inv(hess.end(), n_vars(), n_vars(), false);
      problem.log_integrand_hess(mu.memptr(), hess.memptr(), mem);
      hess.for_each([](double &res){ res *= -1; });
      for(size_t i = 0; i < n_vars(); ++i)
        hess(i, i) += 1;

      if((succeeded = arma::inv_sympd(hess_inv, hess))){
        succeeded = arma::chol(C, hess_inv);

        sq_C_deter = 1;
        for(arma::uword i = 0; i < C.n_cols; ++i)
          sq_C_deter *= C(i, i);
      }
    }

    if(!succeeded){
      // perform the computation with a non-adaptive version
      mu.zeros(n_vars());
      C.zeros(n_vars(), n_vars());
      C.diag() += 1;
      sq_C_deter = 1;
    }
    mem.reset_to_mark();
  }

void adaptive_problem::eval
  (double const *points, size_t const n_points, double * __restrict__ outs,
   simple_mem_stack<double> &mem) const {
  /// transform the points
  double * const __restrict__
    points_trans{mem.get(n_vars() * n_points + n_points)};
  double * const __restrict__ fac{points_trans + n_vars() * n_points};

  // do the matrix product points.C
  std::copy(points, points + n_points * n_vars(), points_trans);
  {
    int const m = n_points, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_N, &c_N, &m, &n, &alpha, C.memptr(), &n,
       points_trans, &m, 1, 1, 1, 1);
  }

  // add the mode
  for(size_t j = 0; j < n_vars(); ++j)
    std::for_each
      (points_trans + j * n_points, points_trans + (j + 1) * n_points,
       [&](double &lhs){ lhs += mu[j]; });

  // evaluate the inner part
  auto mem_marker = mem.set_mark_raii(); // problem may turn back mem
  problem.eval(points_trans, n_points, outs, mem);

  // add the additional weight
  std::fill(fac, fac + n_points, 0);
  for(size_t j = 0; j < n_vars(); ++j){
    size_t const offset{j * n_points};
    for(size_t i = 0; i < n_points; ++i)
      fac[i] +=
        points[i + offset] * points[i + offset]
        - points_trans[i + offset] * points_trans[i + offset];
  }

  std::for_each
    (fac, fac + n_points,
     [&](double &res) { res = std::exp(res / 2) * sq_C_deter; });

  for(size_t j = 0; j < n_out(); ++j)
    for(size_t i = 0; i < n_points; ++i)
      outs[i + j * n_points] *= fac[i];
}

// recursive functions needed for quadrature implementation
namespace {
void ghq_fill_fixed
  (size_t const lvl, double * const points, double * const weights,
   size_t const n_points, ghq_data const &dat){
  // how many times should we repeat each node?
  size_t const n_nodes{dat.n_nodes};
  size_t n_rep{1};
  for(size_t i = 1; i < lvl; ++i)
    n_rep *= n_nodes;

  // fill the weights and points
  double *p{points}, *w{weights};
  for(size_t j = 0; j < n_points;)
      for(size_t n = 0; n < n_nodes; ++n)
        for(size_t i = 0; i < n_rep; ++i, ++j){
          *p++  = dat.nodes[n];
          *w++ *= dat.weights[n];
        }

  if(lvl > 1)
    ghq_fill_fixed(lvl - 1, points + n_points, weights, n_points, dat);
}

void ghq_inner
  (double * __restrict__ res, size_t const n_res, double * const outs,
   size_t const lvl, size_t const idx_fix, size_t const n_points,
   size_t const n_vars, double * const points, double const * weights,
   ghq_problem const &problem, ghq_data const &dat,
   simple_mem_stack<double> &mem){
  if(lvl == idx_fix){
    // evaluate the integrand and add the result
    problem.eval(points, n_points, outs, mem);
    mem.reset_to_mark();

    for(size_t i = 0; i < n_res; ++i)
      for(size_t j = 0; j < n_points; ++j)
        res[i] += weights[j] * outs[j + i * n_points];

    return;
  }

  // we have to go through all the configurations recursively
  double * const __restrict__ weights_scaled{mem.get(n_points)};
  auto mem_marker = mem.set_mark_raii();

  size_t const n_nodes{dat.n_nodes};
  for(size_t j  = 0; j < n_nodes; ++j){
    double * const p{points + (n_vars - lvl) * n_points};
    for(size_t i = 0; i < n_points; ++i){
      weights_scaled[i] = dat.weights[j] * weights[i];
      p[i] = dat.nodes[j];
    }

    // run the next level
    ghq_inner(res, n_res, outs, lvl - 1, idx_fix, n_points, n_vars, points,
              weights_scaled, problem, dat, mem);
  }
}
} // namespace

void ghq
  (double * __restrict__ res, ghq_data const &ghq_data_in,
   ghq_problem const &problem, simple_mem_stack<double> &mem,
   size_t const target_size){
  size_t const n_nodes{ghq_data_in.n_nodes},
               n_vars{problem.n_vars()},
               n_out{problem.n_out()};

  // checks
  if(n_out < 1)
    return;
  else if(n_nodes < 1)
    throw std::invalid_argument("n_nodes < 1");
  else if(n_vars < 1)
    throw std::invalid_argument("n_vars < 1");

  // determine the maximum number of points we will use and the "fixed" level
  size_t idx_fix{1};
  size_t n_points{n_nodes};
  for(; n_points * n_nodes < target_size && idx_fix < n_vars; ++idx_fix)
    n_points *= n_nodes;

  // get the memory we need
  double * const points
    {mem.get(2 * n_nodes + n_points * (1 + n_vars + n_out))},
         * const outs{points + n_points * n_vars},
         * const weights{outs + n_points * n_out},
         * const ghq_nodes{weights + n_points},
         * const ghq_weigths{ghq_nodes + n_nodes};

  auto mem_marker = mem.set_mark_raii();

  // initialize the objects before the computation
  std::fill(weights, weights + n_points, 1);
  std::fill(res, res + n_out, 0);

  for(size_t i = 0; i < n_nodes; ++i){
    ghq_nodes[i] = ghq_data_in.nodes[i] * 1.4142135623731;  // sqrt(2)
    ghq_weigths[i] = ghq_data_in.weights[i] * 0.564189583547756; // 1 / sqrt(pi)
  }

  ghq_data const ghq_data_use{ghq_nodes, ghq_weigths, n_nodes};

  // the points matrix has a "fixed part" that never changes and set the
  // corresponding weights
  ghq_fill_fixed
    (idx_fix, points + n_points * (n_vars - idx_fix), weights, n_points,
     ghq_data_use);

  ghq_inner(res, n_out, outs, n_vars, idx_fix, n_points, n_vars, points,
            weights, problem, ghq_data_use, mem);
}

combined_problem::combined_problem
  (std::vector<ghq_problem const *> const &problems):
  problems{problems} {
    if(problems.size() > 0){
      size_t const n_vars_first{problems[0]->n_vars()};
      for(ghq_problem const * p : problems){
        if(p->n_vars() != n_vars_first)
          throw std::invalid_argument("p->n_vars() != n_vars_first");
        else if(p->n_out() < 1)
          throw std::invalid_argument("p->n_out() < 1");
      }
    }
  }

void combined_problem::eval
  (double const *points, size_t const n_points, double * __restrict__ outs,
   simple_mem_stack<double> &mem) const {
  double * const __restrict__ scales{mem.get(n_points * (1 + n_out_inner))},
         * const __restrict__ outs_inner{scales + n_points};
  auto mem_marker = mem.set_mark_raii();

  // call eval on each of the problems while setting the value of the integrand
  double * const integrands{outs};
  outs += n_points; // outs are now the derivatives
  std::fill(integrands, integrands + n_points, 1);
  {
    double * outs_inner_p{outs_inner};
    size_t pi{};
    for(auto p : problems){
      p->eval(points, n_points, outs_inner_p, mem);

      for(size_t i = 0; i < n_points; ++i)
        integrands[i] *= outs_inner_p[i];
      outs_inner_p += n_outs[pi] * n_points;
      ++pi;
    }
  }

  // compute the derivatives
  double const * outs_inner_p{outs_inner};
  for(size_t const n_outs_p : n_outs){
    if(n_outs_p > 1){
      // compute the scales to use for the derivatives
      for(size_t i = 0; i < n_points; ++i, ++outs_inner_p)
        scales[i] = integrands[i] > 0 ? integrands[i] / *outs_inner_p : 0;

      // set the derivatives
      for(size_t j = 0; j < n_outs_p - 1; ++j)
        for(size_t i = 0; i < n_points; ++i, ++outs_inner_p, ++outs)
          *outs = *outs_inner_p * scales[i];

    } else
      outs_inner_p += n_points;
  }
}

double combined_problem::log_integrand
  (double const *point, simple_mem_stack<double> &mem) const {
  double out{};
  for(auto p : problems)
    out += p->log_integrand(point, mem);
  return out;
}

double combined_problem::log_integrand_grad
  (double const *point, double * __restrict__ grad,
   simple_mem_stack<double> &mem) const {
  double * const grad_inner{mem.get(n_vars())};
  auto mem_marker = mem.set_mark_raii();

  std::fill(grad, grad + n_vars(), 0);
  double out{};
  for(auto p : problems){
    out += p->log_integrand_grad(point, grad_inner, mem);
    for(size_t i = 0; i < n_vars(); ++i)
      grad[i] += grad_inner[i];
  }
  return out;
}

void combined_problem::log_integrand_hess
  (double const *point, double *hess,
   simple_mem_stack<double> &mem) const {
  size_t const n_vars_sq{n_vars() * n_vars()};
  double * const hess_inner{mem.get(n_vars_sq)};
  auto mem_marker = mem.set_mark_raii();

  std::fill(hess, hess + n_vars_sq, 0);
  for(auto p : problems){
    p->log_integrand_hess(point, hess_inner, mem);
    for(size_t i = 0; i < n_vars_sq; ++i)
      hess[i] += hess_inner[i];
  }
}

outer_prod_problem::outer_prod_problem(size_t const n_vars):
  v_n_vars{n_vars} { }

void outer_prod_problem::eval
  (double const *points, size_t const n_points, double * __restrict__ outs,
   simple_mem_stack<double> &mem) const {
  // have till fill in the ones
  std::fill(outs, outs + n_points, 1);
  outs += n_points;

  // compute the outer products
  size_t out_offset{};
  for(size_t j = 0; j < n_vars(); ++j)
    for(size_t k = 0; k <= j; ++k, ++out_offset)
        for(size_t i = 0; i < n_points; ++i)
          outs[i + out_offset * n_points] =
            points[i + k * n_points] * points[i + j * n_points];
}

double outer_prod_problem::log_integrand
  (double const*, simple_mem_stack<double>&) const {
  return 0;
}

double outer_prod_problem::log_integrand_grad
  (double const*, double * __restrict__ gr, simple_mem_stack<double>&) const {
  std::fill(gr, gr + n_vars(), 0);
  return 0;
}

void outer_prod_problem::log_integrand_hess
  (double const*, double * __restrict__ hess, simple_mem_stack<double>&) const {
  std::fill(hess, hess + n_vars() * n_vars(), 0);
}

void outer_prod_problem::d_Sig
  (double * __restrict__ res, double const *out, double const integral,
   arma::mat const &Sigma) const {
  arma::mat outer_int(n_vars(), n_vars());
  for(arma::uword j = 0; j < outer_int.n_cols; ++j, ++out){
    for(arma::uword i = 0; i < j; ++i, ++out){
      outer_int(i, j) = *out / 2;
      outer_int(j, i) = *out / 2;
    }
    outer_int(j, j) = (*out - integral) / 2;
  }

  arma::mat const sig_chol = arma::chol(Sigma);
  arma::mat lhs(res, n_vars(), n_vars(), false);

  lhs = arma::solve
    (arma::trimatu(sig_chol),
     arma::solve(arma::trimatu(sig_chol), outer_int).t());
}

} // namespace ghqCpp
