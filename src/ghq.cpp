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

    double *bfgs_mem{mem.get(PSQN::bfgs_n_wmem(n_vars()))};
    auto bfgs_marker = mem.set_mark_raii();
    auto res = PSQN::bfgs
      (my_mode_problem, mu.memptr(), bfgs_mem, rel_eps, max_it, c1, c2, 0L,
       gr_tol);

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

void adaptive_problem::post_process
  (double *res, simple_mem_stack<double> &mem) const {
  problem.post_process(res, mem);
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

    char const trans{'T'};
    int const i_n_points = n_points, i_n_res = n_res, incxy{1};
    double const alpha{1};

    F77_CALL(dgemv)
      (&trans, &i_n_points, &i_n_res, &alpha, outs, &i_n_points, weights,
       &incxy, &alpha, res, &incxy, size_t(1));

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

  problem.post_process(res, mem);
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

void combined_problem::post_process
(double *res, simple_mem_stack<double> &mem) const {
  double const integral{res[0]};
  res += 1;

  for(auto p : problems){
    size_t const p_n_out{p->n_out()};
    if(p_n_out < 2)
      continue;

    double * const cp{mem.get(p_n_out)};
    auto cp_marker = mem.set_mark_raii();
    cp[0] = integral;
    std::copy(res, res + p_n_out - 1, cp + 1);
    p->post_process(cp, mem);

    if(cp[0] != integral)
      throw std::runtime_error("post_process changed the first element");
    std::copy(cp + 1, cp + p_n_out, res);
    res += p_n_out - 1;
  }
}

template<bool comp_grad>
rescale_problem<comp_grad>::rescale_problem
  (arma::mat const &Sigma, ghq_problem const &inner_problem):
  Sigma_chol{arma::chol(Sigma)}, inner_problem{inner_problem} {
    if(n_out_inner < 1)
      throw std::invalid_argument("n_out_inner < 1");
    if(inner_problem.n_vars() != n_vars())
      throw std::invalid_argument("inner_problem.n_vars() != n_vars()");
  }

template<bool comp_grad>
void rescale_problem<comp_grad>::eval
  (double const *points, size_t const n_points, double * __restrict__ outs,
   simple_mem_stack<double> &mem) const {
  double * const __restrict__ us{mem.get(n_points * n_vars())};
  auto mem_marker = mem.set_mark_raii();

  // do the matrix product points.chol(Sigma)
  std::copy(points, points + n_points * n_vars(), us);
  {
    int const m = n_points, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_N, &c_N, &m, &n, &alpha, Sigma_chol.memptr(), &n,
       us, &m, 1, 1, 1, 1);
  }

  inner_problem.eval(us, n_points, outs, mem);

  if constexpr(comp_grad){
    // compute the weighted outer products
    double const * const integrands{outs};
    outs += n_points * n_out_inner;
    size_t out_offset{};
    for(size_t j = 0; j < n_vars(); ++j)
      for(size_t k = 0; k <= j; ++k, ++out_offset)
        for(size_t i = 0; i < n_points; ++i)
          outs[i + out_offset * n_points] =
            integrands[i] * points[i + k * n_points] * points[i + j * n_points];
  }
}

template<bool comp_grad>
double * rescale_problem<comp_grad>::rescale
  (double const *point, simple_mem_stack<double> &mem) const {
  double * const __restrict__ u{mem.get(n_vars())};

  std::copy(point, point + n_vars(), u);
  {
    int const m = 1, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_N, &c_N, &m, &n, &alpha, Sigma_chol.memptr(), &n,
       u, &m, 1, 1, 1, 1);
  }
  return u;
}

template<bool comp_grad>
double rescale_problem<comp_grad>::log_integrand
  (double const *point, simple_mem_stack<double> &mem) const {
  auto u = rescale(point, mem);
  auto mem_marker = mem.set_mark_raii();
  return inner_problem.log_integrand(u, mem);
}

template<bool comp_grad>
double rescale_problem<comp_grad>::log_integrand_grad
  (double const *point, double * __restrict__ gr,
   simple_mem_stack<double> &mem) const {
  auto u = rescale(point, mem);
  auto mem_marker = mem.set_mark_raii();
  double const out{inner_problem.log_integrand_grad(u, gr, mem)};

  // compute gr <- Sigma_chol.gr
  {
    int const m = 1, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_L{'L'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_L, &c_U, &c_N, &c_N, &n, &m, &alpha, Sigma_chol.memptr(), &n,
       gr, &n, 1, 1, 1, 1);
  }

  return out;
}

template<bool comp_grad>
void rescale_problem<comp_grad>::log_integrand_hess
  (double const *point, double * __restrict__ hess,
   simple_mem_stack<double>&mem) const {
  auto u = rescale(point, mem);
  auto mem_marker = mem.set_mark_raii();
  inner_problem.log_integrand_hess(u, hess, mem);

  // compute hess <- Sigma_chol.hess.Sigma_chol^T
  {
    int const m = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'}, c_T{'T'}, c_L{'L'};
    F77_CALL(dtrmm)
      (&c_L, &c_U, &c_N, &c_N, &m, &m, &alpha, Sigma_chol.memptr(), &m,
       hess, &m, 1, 1, 1, 1);
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_T, &c_N, &m, &m, &alpha, Sigma_chol.memptr(), &m,
       hess, &m, 1, 1, 1, 1);
  }
}

template<bool comp_grad>
void rescale_problem<comp_grad>::post_process
  (double *res, simple_mem_stack<double> &mem) const {
  inner_problem.post_process(res, mem);
  if constexpr(!comp_grad)
    return;

  double const integral{res[0]};
  res += inner_problem.n_out();

  arma::mat outer_int(n_vars(), n_vars());
  {
    double * res_ij{res};
    for(arma::uword j = 0; j < outer_int.n_cols; ++j, ++res_ij){
      for(arma::uword i = 0; i < j; ++i, ++res_ij){
        outer_int(i, j) = *res_ij / 2;
        outer_int(j, i) = *res_ij / 2;
      }
      outer_int(j, j) = (*res_ij - integral) / 2;
    }
  }

  arma::mat lhs(res, n_vars(), n_vars(), false, true);
  lhs = arma::solve
    (arma::trimatu(Sigma_chol),
     arma::solve(arma::trimatu(Sigma_chol), outer_int).t());
}

template class rescale_problem<false>;
template class rescale_problem<true>;

template<bool comp_grad>
rescale_shift_problem<comp_grad>::rescale_shift_problem
  (arma::mat const &Sigma, arma::vec const &m, ghq_problem const &inner_problem):
  m{m}, Sigma_chol{arma::chol(Sigma)}, inner_problem{inner_problem} {
    if(n_out_inner < 1)
      throw std::invalid_argument("n_out_inner < 1");
    if(inner_problem.n_vars() != n_vars())
      throw std::invalid_argument("inner_problem.n_vars() != n_vars()");
    if(m.n_elem != Sigma_chol.n_rows)
      throw std::invalid_argument("m.n_elem != Sigma_chol.n_rows");
  }

template<bool comp_grad>
void rescale_shift_problem<comp_grad>::eval
  (double const *points, size_t const n_points, double * __restrict__ outs,
   simple_mem_stack<double> &mem) const {
  double * const __restrict__ us{mem.get(n_points * n_vars())};
  auto mem_marker = mem.set_mark_raii();

  // do the matrix product points.chol(Sigma)
  std::copy(points, points + n_points * n_vars(), us);
  {
    int const m = n_points, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_N, &c_N, &m, &n, &alpha, Sigma_chol.memptr(), &n,
       us, &m, 1, 1, 1, 1);
  }

  for(size_t j = 0; j < n_vars(); ++j)
    for(size_t i = 0; i < n_points; ++i)
      us[i + j * n_points] += m[j];

  inner_problem.eval(us, n_points, outs, mem);

  if constexpr(comp_grad){
    // compute the weighted points
    double const * const integrands{outs};
    outs += n_points * n_out_inner;
    for(size_t j = 0; j < n_vars(); ++j)
      for(size_t i = 0; i < n_points; ++i)
        outs[i + j * n_points] = integrands[i] *  points[i + j * n_points];

    // compute the weighted outer products
    outs += n_points * n_vars();
    size_t out_offset{};
    for(size_t j = 0; j < n_vars(); ++j)
      for(size_t k = 0; k <= j; ++k, ++out_offset)
        for(size_t i = 0; i < n_points; ++i)
          outs[i + out_offset * n_points] =
            integrands[i] * points[i + k * n_points] * points[i + j * n_points];
  }
}

template<bool comp_grad>
double * rescale_shift_problem<comp_grad>::rescale_center
  (double const *point, simple_mem_stack<double> &mem) const {
  double * const __restrict__ u{mem.get(n_vars())};

  std::copy(point, point + n_vars(), u);
  {
    int const m = 1, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_N, &c_N, &m, &n, &alpha, Sigma_chol.memptr(), &n,
       u, &m, 1, 1, 1, 1);
  }

  for(size_t i = 0; i < n_vars(); ++i)
    u[i] += m[i];

  return u;
}

template<bool comp_grad>
double rescale_shift_problem<comp_grad>::log_integrand
  (double const *point, simple_mem_stack<double> &mem) const {
  auto u = rescale_center(point, mem);
  auto mem_marker = mem.set_mark_raii();
  return inner_problem.log_integrand(u, mem);
}

template<bool comp_grad>
double rescale_shift_problem<comp_grad>::log_integrand_grad
  (double const *point, double * __restrict__ gr,
   simple_mem_stack<double> &mem) const {
  auto u = rescale_center(point, mem);
  auto mem_marker = mem.set_mark_raii();
  double const out{inner_problem.log_integrand_grad(u, gr, mem)};

  // compute gr <- Sigma_chol.gr
  {
    int const m = 1, n = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_L{'L'}, c_U{'U'}, c_N{'N'};
    F77_CALL(dtrmm)
      (&c_L, &c_U, &c_N, &c_N, &n, &m, &alpha, Sigma_chol.memptr(), &n,
       gr, &n, 1, 1, 1, 1);
  }

  return out;
}

template<bool comp_grad>
void rescale_shift_problem<comp_grad>::log_integrand_hess
  (double const *point, double * __restrict__ hess,
   simple_mem_stack<double>&mem) const {
  auto u = rescale_center(point, mem);
  auto mem_marker = mem.set_mark_raii();
  inner_problem.log_integrand_hess(u, hess, mem);

  // compute hess <- Sigma_chol.hess.Sigma_chol^T
  {
    int const m = n_vars();
    constexpr double const alpha{1};
    constexpr char const c_R{'R'}, c_U{'U'}, c_N{'N'}, c_T{'T'}, c_L{'L'};
    F77_CALL(dtrmm)
      (&c_L, &c_U, &c_N, &c_N, &m, &m, &alpha, Sigma_chol.memptr(), &m,
       hess, &m, 1, 1, 1, 1);
    F77_CALL(dtrmm)
      (&c_R, &c_U, &c_T, &c_N, &m, &m, &alpha, Sigma_chol.memptr(), &m,
       hess, &m, 1, 1, 1, 1);
  }
}

template<bool comp_grad>
void rescale_shift_problem<comp_grad>::post_process
  (double *res, simple_mem_stack<double> &mem) const {
  inner_problem.post_process(res, mem);
  if constexpr(!comp_grad)
    return;

  double const integral{res[0]};
  res += inner_problem.n_out();

  {
    arma::vec rhs(res, n_vars()), lhs(res, n_vars(), false, true);
    lhs = arma::solve(arma::trimatu(Sigma_chol), rhs);
  }

  res += n_vars();

  arma::mat outer_int(n_vars(), n_vars());
  {
    double * res_ij{res};
    for(arma::uword j = 0; j < outer_int.n_cols; ++j, ++res_ij){
      for(arma::uword i = 0; i < j; ++i, ++res_ij){
        outer_int(i, j) = *res_ij / 2;
        outer_int(j, i) = *res_ij / 2;
      }
      outer_int(j, j) = (*res_ij - integral) / 2;
    }
  }
  arma::mat lhs(res, n_vars(), n_vars(), false, true);

  lhs = arma::solve
    (arma::trimatu(Sigma_chol),
     arma::solve(arma::trimatu(Sigma_chol), outer_int).t());
}

template class rescale_shift_problem<false>;
template class rescale_shift_problem<true>;

} // namespace ghqCpp
