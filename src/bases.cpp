#include "bases.h"
#include <algorithm> // lower_bound
#include <cmath> // isnan
#include <stdexcept> // invalid_argument

inline void check_splines
  (const arma::vec &boundary_knots, const arma::vec &interior_knots,
   const vajoint_uint order) {
  if(order<1)
    throw std::invalid_argument("order<1");
  if(boundary_knots.size() != 2L)
    throw std::invalid_argument("boundary_knots should have length 2");
  if(interior_knots.size()>0 && boundary_knots(0)>min(interior_knots))
    throw std::invalid_argument("boundary_knots(0)>min(interior_knots)");
  if(interior_knots.size()>0 && boundary_knots(1)<max(interior_knots))
    throw std::invalid_argument("boundary_knots(1)<max(interior_knots)");
  // TODO: check if interior_knots are in ascending order?
}

inline void throw_invalid_out(
    std::string const &cl, vajoint_uint const dim, vajoint_uint const dim_ex){
  std::stringstream msg;
  msg << cl << ": invalid 'out' (dim is " << dim << "; expected "
      << dim_ex << ')';
  throw std::invalid_argument(msg.str());
}

namespace joint_bases {

SplineBasis::SplineBasis(const vec &knots, const vajoint_uint order,
                         bool const use_log, bool const with_integral):
  basisMixin(use_log),
  order(order),
  knots{
    ([](vec ks){
      // sort the knots
      std::sort(ks.begin(), ks.end());
      return ks;
    })(knots)
  },
  integral_basis {
    ([&](){
      if(!with_integral)
        return std::unique_ptr<SplineBasis>();

      // append an additional knot
      arma::vec new_ks(knots.n_elem + 1);
      std::copy(knots.begin(), knots.end(), new_ks.begin());
      if(knots.size() > 0)
        new_ks[knots.n_elem] = knots[knots.n_elem - 1];

      return std::make_unique<SplineBasis>(new_ks, order + 1, use_log, false);
    })()
  } {
  if (order<1)
    throw std::invalid_argument("order<1");
  // TODO: should make a check on the number of knots?
}

inline arma::vec SplineBasis_knots
(vec const &boundary_knots, vec const &interior_knots,
 vajoint_uint const order) {
  check_splines(boundary_knots, interior_knots, order);

  vajoint_uint const nknots = interior_knots.size() + 2 * order;
  vec knots(nknots);
  for(vajoint_uint i = 0; i < order; i++) {
    knots[i] = boundary_knots[0];
    knots[nknots - i - 1] = boundary_knots[1];
  }
  if (interior_knots.size() > 0)
    for(vajoint_uint i = 0; i < interior_knots.size(); i++)
      knots[i + order] = interior_knots[i];

  return knots;
}

bs::bs(const vec &bk, const vec &ik, const bool inter, const vajoint_uint ord,
       bool const use_log):
  SplineBasis(SplineBasis_knots(bk, ik, ord), ord, use_log),
  boundary_knots{bk[0], bk[1]},
  intercept(inter),
  df((int)intercept + order - 1 + ik.size()) {
  check_splines(bk, ik, order);
}

ns::ns(const vec &bk, const vec &interior_knots,
       const bool intercept, const vajoint_uint order,
       const bool use_log):
  basisMixin(use_log),
  s_basis{SplineBasis_knots(bk, interior_knots, order), order},
  boundary_knots{bk[0], bk[1]},
  intercept(intercept)
  {
    {
      arma::vec tmp{boundary_knots[0], boundary_knots[1]};
      mat const_basis = s_basis.basis
        (tmp, wmem::mem_stack().get(s_basis.n_wmem()), 2);
      if (!intercept)
        const_basis = const_basis.cols(1, const_basis.n_cols - 1);
      inplace_trans(const_basis);

      int const m = const_basis.n_rows;
      constexpr int n{2};
      qr_A.resize(const_basis.n_elem);
      std::copy(const_basis.begin(), const_basis.end(), qr_A.begin());
      qr_tau.resize(std::min(m, n));

      int lwork{-1}, info{};
      qr_jpvt[0] = 0;
      qr_jpvt[1] = 0;
      {
        double opt_work{};
        F77_CALL(dgeqp3)
          (&m, &n, qr_A.data(), &m, qr_jpvt.data(), qr_tau.data(), &opt_work,
           &lwork, &info);
        lwork = opt_work;
      }

      F77_CALL(dgeqp3)
        (&m, &n, qr_A.data(), &m, qr_jpvt.data(), qr_tau.data(),
         wmem::mem_stack().get(lwork), &lwork, &info);

      if(info < 0)
        throw std::invalid_argument("ns: QR decomposition failed");
    }

    tl0 = trans
      (s_basis(boundary_knots[0], wmem::mem_stack().get(s_basis.n_wmem()),
               nullptr, 0.));
    tl1 = trans
      (s_basis(boundary_knots[0], wmem::mem_stack().get(s_basis.n_wmem()),
               nullptr, 1.));
    tr0 = trans
      (s_basis(boundary_knots[1], wmem::mem_stack().get(s_basis.n_wmem()),
               nullptr, 0.));
    tr1 = trans
      (s_basis(boundary_knots[1], wmem::mem_stack().get(s_basis.n_wmem()),
               nullptr, 1.));
  }

iSpline::iSpline(const vec &boundary_knots, const vec &interior_knots,
                 const bool intercept, const vajoint_uint order):
  intercept(intercept), order(order),
  bspline(boundary_knots, interior_knots, false, order + 1) { }


mSpline::mSpline(const vec &boundary_knots, const vec &interior_knots,
                 const bool intercept, const vajoint_uint order) :
  bspline(boundary_knots, interior_knots, true, order),
  intercept(intercept) { }

orth_poly::orth_poly(vajoint_uint const degree, bool const intercept,
                     bool const use_log):
basisMixin(use_log),
alpha(), norm2(), raw(true), intercept(intercept),
n_basis_v(degree + intercept) { }

orth_poly::orth_poly(vec const &alpha, vec const &norm2,
                     bool const intercept, bool const use_log):
basisMixin(use_log),
alpha(alpha), norm2(norm2), raw(false), intercept(intercept),
n_basis_v(norm2.size() - 2 + intercept),
orth_map(((alpha.size() + 1) * (alpha.size() + 2)) / 2) {
  for(vajoint_uint i = 0; i < norm2.size(); ++i)
    if(norm2[i] <= 0.)
      throw std::invalid_argument("invalid norm2");
  if(alpha.n_elem + 2L != norm2.n_elem)
    throw std::invalid_argument("invalid alpha");

  if(!raw){
    vajoint_uint nc{alpha.size() + 1};
    auto g = orth_map.begin();
    *g++ = 1;

    if(nc > 1){
      *g++ = -alpha[0] * *(g - 1);
      *g++ = 1;
    }
    if(nc > 2){
      for(vajoint_uint i = 2; i < nc; ++i){
        auto g_prev = g - i;
        auto g_prev_m1 = g_prev;
        auto g_old = g_prev - i + 1;
        vajoint_uint j{};
        double const sig_ratio{norm2[i] / norm2[i - 1]};
        for(; j <= i - 2; ++j, ++g){
          *g = -sig_ratio * *g_old++ - alpha[i - 1] * *g_prev++;
          if(j > 0)
            *g += *g_prev_m1++;
        }

        *g++ += *g_prev_m1++ - alpha[i - 1] * *g_prev;
        *g++ += *g_prev_m1;
      }
    }

    g = orth_map.begin() + 1;
    for(vajoint_uint i = 1; i < nc; ++i){
      double const denom{std::sqrt(norm2[i + 1])};
      for(vajoint_uint j = 0; j <= i; ++j)
        *g++ /= denom;
    }
  }
}

orth_poly orth_poly::poly_basis(vec x, uword const degree, mat &X){
  vajoint_uint const n = x.n_elem,
                    nc = degree + 1;
  double const x_bar = mean(x);
  x -= x_bar;
  mat XX(n, nc);
  XX.col(0).ones();
  for(vajoint_uint d = 1L; d < nc; d++){
    double       * xx_new = XX.colptr(d);
    double const * xx_old = XX.colptr(d - 1);
    for(vajoint_uint i = 0; i < n; ++i, ++xx_new, ++xx_old)
      *xx_new = *xx_old * x[i];
  }

  mat R;
  if(!qr_econ(X, R, XX))
    /* TODO: can be done smarter by calling LAPACK or LINPACK directly */
    throw std::runtime_error(
        "orth_poly::poly_basis(): QR decomposition failed");

  for(vajoint_uint c = 0; c < nc; ++c)
    X.col(c) *= R.at(c, c);

  vec norm2(nc + 1),
  alpha(nc - 1);
  norm2[0] = 1;
  for(vajoint_uint c = 0; c < nc; ++c){
    double z_sq(0),
    x_z_sq(0);
    double const *X_i = X.colptr(c);
    for(vajoint_uint i = 0; i < n; ++i, ++X_i){
      double const z_sq_i = *X_i * *X_i;
      z_sq += z_sq_i;
      if(c < degree)
        x_z_sq += x[i] * z_sq_i;
    }
    norm2[c + 1] = z_sq;
    if(c < degree)
      alpha[c] = x_z_sq / z_sq + x_bar;
  }

  orth_poly out(alpha, norm2, true);
  for(vajoint_uint j = 1; j < nc; ++j)
    for(vajoint_uint i = 0; i < n; ++i)
      X.at(i, j) /= out.sqrt_norm2.at(j + 1);
  return out;
}

} // namespace joint_bases
