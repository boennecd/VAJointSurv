#ifndef GHQ_H
#define GHQ_H

#include "simple_mem_stack.h"
#include "arma-wrap.h"
#include <algorithm>
#include <stdexcept>
#include <psqn-bfgs.h>

namespace ghqCpp {

/**
 * virtual base class for a Gauss-Hermite quadrature problem. It specifies
 * the dimension of the random effects and the number of outputs. The integral
 * we approximate are of the form
 *
 *   int phi(x)g(x)dx
 *
 * where phi(x) is a given dimensional standard multivariate normal
 * distribution. The class also computes g(x).
 *
 * To perform adaptive Gauss-Hermite quadrature for one of the elements of g,
 * say g_i, the class also has member functions to compute log g_i(x), the
 * gradient of it, and the Hessian.
 */
struct ghq_problem {
  /// the number of variables
  virtual size_t n_vars() const = 0;
  /// the number of output
  virtual size_t n_out() const = 0;

  /**
   * n_points x n_vars() points at which to evaluate g is passed in points. The
   * output of evaluating g at each point is then written to n_points x n_out()
   * matrix outs. The mem object can be used for working memory.
   */
  virtual void eval
    (double const *points, size_t const n_points, double * __restrict__ outs,
     simple_mem_stack<double> &mem) const = 0;

  /// evaluates log g_i for the element chosen for the adaptive method
  virtual double log_integrand
    (double const *point, simple_mem_stack<double> &mem) const {
    throw std::runtime_error("not implemented");
    return 0;
  }

  /**
   * evaluates log g_i and the gradient for the element chosen for the adaptive
   * method
   */
  virtual double log_integrand_grad
    (double const *point, double * __restrict__ grad,
     simple_mem_stack<double> &mem) const {
    throw std::runtime_error("not implemented");
    return 0;
  }

  /// evaluates the Hessian of log g_i
  virtual void log_integrand_hess
    (double const *point, double *hess,
     simple_mem_stack<double> &mem) const {
    throw std::runtime_error("not implemented");
  }

  virtual ~ghq_problem() = default;
};

/// Gauss-Hermite quadrature nodes and weights
struct ghq_data {
  double const * nodes, * weights;
  size_t n_nodes;
};

/**
 * performs Gauss-Hermite quadrature. The target_size is the maximum number of
 * integrands to simultaneously process. However, at least the number of
 * quadrature nodes is simultaneously processed.
 *
 * The res pointer needs to have enough memory for problem.n_out().
 */
void ghq
  (double * __restrict__ res, ghq_data const &ghq_data_in,
   ghq_problem const &problem, simple_mem_stack<double> &mem,
   size_t const target_size = 128);

/// overload where the user does not pre-allocated memory
inline std::vector<double> ghq
  (ghq_data const &ghq_data_in, ghq_problem const &problem,
   simple_mem_stack<double> &mem, size_t const target_size = 128){
  std::vector<double> out(problem.n_out());
  ghq(out.data(), ghq_data_in, problem, mem, target_size);
  return out;
}

/**
 * Takes the integral
 *
 *   int phi(x)g(x)dx
 *
 * and makes is adaptive by working with
 *
 *   int phi(x, mu, Psi)phi(x, mu, Psi)^(-1)phi(x)g(x)dx
 *
 * where mu is the mode of the log of the integrand and Psi is the inverse of
 * the negative Hessian of the log of the integrand. This is then transformed to
 *
 *   |Psi|^(1/2)int phi(x)[phi(x)^(-1)phi(mu + C.x)g(mu + C.x)]dx
 *                        |-----------------------------------|
 *                                        h(x)
 *
 * where C.C^T = Psi is the Cholesky decomposition of Psi.
 */
class adaptive_problem final : public ghq_problem  {
  ghq_problem const &problem;
  size_t const v_n_vars{problem.n_vars()},
               v_n_out{problem.n_out()};

  /**
   * the Cholesky decomposition of the -H^(-1) where H is the Hessian at the
   * mode
   */
  arma::mat C;
  /// the mode
  arma::vec mu;
  /// the square root of the determinant of C^T.C
  double sq_C_deter;

public:
  /// class used to find the mode
  class mode_problem final : public PSQN::problem {
    ghq_problem const &problem;
    simple_mem_stack<double> &mem;
    PSQN::psqn_uint const v_n_vars = problem.n_vars();

  public:
    mode_problem(ghq_problem const &problem, simple_mem_stack<double> &mem);

    PSQN::psqn_uint size() const { return v_n_vars; }
    double func(double const *val);

    double grad(double const * __restrict__ val,
                double       * __restrict__ gr);
  };

  adaptive_problem
    (ghq_problem const &problem, simple_mem_stack<double> &mem,
     double const rel_eps = 1e-6, PSQN::psqn_uint const max_it = 1000L,
     double const c1 = .0001, double const c2 = .9, double const gr_tol = -1);

  size_t n_vars() const { return v_n_vars; }
  size_t n_out() const { return v_n_out; }

  void eval
  (double const *points, size_t const n_points, double * __restrict__ outs,
   simple_mem_stack<double> &mem) const;
};

/**
 * The class takes multiple problems g_1(x), g_2(x), ..., g_l(x) and assumes
 * that
 *
 *  1. the first entries, g_(11)(x), g_(21)(x), ..., g_(l1)(x), are for the
 *     integral
 *
 *       A = int phi(x) prod_(i = 1)^l g_(i1)(x) dx
 *
 *  2. the remaining entries are derivatives in the form of
 *
 *       g_i' = d/dz_i g_i(x; z_i).
 *
 * for some fixed vector z_i. The function then returns the estimator of A in
 * the first element and the derivatives of A for each z_i. The latter
 * derivatives are stored in the order the problems are passed and can be
 * computed with
 *
 *   int phi(x) g_j'(x) / g_(j1)(x)prod_(i = 1)^l g_(i1)(x) dx
 */
class combined_problem final : public ghq_problem  {
  std::vector<ghq_problem const *> problems;

  // to avoid the virtual dispatch
  std::vector<size_t> const n_outs
  {
    ([&]{
      std::vector<size_t> out;
      out.reserve(problems.size());
      for(auto p : problems)
        out.emplace_back(p->n_out());
      return out;
    })()
  };

  size_t const v_n_vars{ problems.size() == 0 ? 0 : problems[0]->n_vars() };
  size_t const n_out_inner
    {std::accumulate(
        n_outs.begin(), n_outs.end(), static_cast<size_t>(0))};
  size_t const v_n_out{n_out_inner - problems.size() + 1};

public:
  combined_problem(std::vector<ghq_problem const *> const &problems);

  size_t n_vars() const { return v_n_vars; }
  size_t n_out() const { return v_n_out; }

  void eval
    (double const *points, size_t const n_points, double * __restrict__ outs,
     simple_mem_stack<double> &mem) const;

  double log_integrand
    (double const *point, simple_mem_stack<double> &mem) const;

  double log_integrand_grad
    (double const *point, double * __restrict__ grad,
     simple_mem_stack<double> &mem) const;

  void log_integrand_hess
    (double const *point, double *hess,
     simple_mem_stack<double> &mem) const;
};

/**
 * computes g(x) = (1, x.x^T) where only the upper triangle is stored for the
 * outer product. This is useful when used with combined_problem for problems
 * of the form
 *
 *   int phi(x; 0, Sigma)g_1(x)...g_l(x) dx
 *     = int phi(x)g_1(C.x)...g_l(C.x) dx
 *
 * where C.C^T = Sigma is the Cholesky decomposition. In the framework defined
 * here, we then redfine the g functions as h_i(x) = g_i(C.x).
 *
 * The derivatives w.r.t. Sigma in this case is
 *
 *   1/2 Sigma^(-1)
 *     [int (x.x^T - Sigma) phi(x; 0, Sigma)g_1(x)...g_l(x) dx]Sigma^(-1)
 *    = 1/2 C^(-T)
 *      [int (x.x^T - I) phi(x; 0, Sigma)h_1(x)...h_l(x) dx]C^(-1)
 *
 * so we need to compute
 *
 *   int x.x^T phi(x; 0, Sigma)h_1(x)...h_l(x) dx
 *
 * to compute the gradient w.r.t. Sigma.
 */
class outer_prod_problem final : public ghq_problem  {
  size_t const v_n_vars,
               v_n_out{1 + (v_n_vars * (v_n_vars + 1)) / 2};

public:
  outer_prod_problem(size_t const n_vars);

  size_t n_vars() const { return v_n_vars; }
  size_t n_out() const { return v_n_out; }

  void eval
    (double const *points, size_t const n_points, double * __restrict__ outs,
     simple_mem_stack<double> &mem) const;

  double log_integrand
    (double const *point, simple_mem_stack<double> &mem) const;

  double log_integrand_grad
    (double const *point, double * __restrict__ grad,
     simple_mem_stack<double> &mem) const;

  void log_integrand_hess
    (double const *point, double *hess,
     simple_mem_stack<double> &mem) const;

  /**
   * computes the derivatives w.r.t. Sigma given the last
   * (n_vars() * (n_vars() + 1)) / 2 elements of the final output, the value of
   *  the entire integral, and the covariance matrix Sigma.
   *
   *  The result is stored in a n_vars() x n_vars() matrix.
   */
  void d_Sig
    (double * __restrict__ res, double const *out, double const integral,
     arma::mat const &Sigma) const;
};

} // namespace ghqCpp

#endif
