#ifndef SPLINES_H
#define SPLINES_H

#include "arma-wrap.h"
#include <limits> // quiet_NaN
#include <stdexcept> // invalid_argument

namespace joint_bases {

constexpr unsigned default_order = 4,
                   default_ders  = 0;
constexpr bool const default_intercept = false;

using namespace arma;

class basisMixin {
public:
  virtual unsigned get_n_basis() const = 0;

  virtual void operator()(
      vec &out, double const x, const unsigned ders = default_ders) const = 0;
  inline vec operator()(
      double const x, unsigned const ders = default_ders) const {
    vec out(get_n_basis());
    operator()(out, x, ders);
    return out;
  }
  inline void operator()(
    double *out, double const x, unsigned const ders = default_ders) const {
    arma::vec v(out, get_n_basis(), false);
    operator()(v, x, ders);
  }

  mat basis(
      const vec &x, const unsigned ders = default_ders,
      const double centre = std::numeric_limits<double>::quiet_NaN())
    const  {
#ifdef DO_CHECKS
    if (ders < 0)
      throw std::invalid_argument("ders<0");
#endif
    unsigned const n_basis(get_n_basis()),
    n_x    (x.n_elem);
    rowvec centering =
      (std::isnan(centre) || ders > 0 ?
      zeros(n_basis) : operator()(centre, 0)).t();

    mat out(n_x, n_basis);
    vec wrk(n_basis);
    for (unsigned i = 0; i < n_x; i++){
      operator()(wrk, x[i], ders);
      out.row(i) = wrk.t() - centering;
    }

    return out;
  }

  virtual ~basisMixin() = default;
};

class SplineBasis : public basisMixin {
public:
  unsigned const order = default_order, /* order of the spline */
                 ordm1 = order - 1;     /* order - 1 (3 for cubic splines) */
  vec const knots;	               /* knot vector */
  unsigned const nknots = knots.n_elem, /* number of knots
                                           except for the boundary case */
                  ncoef =               /* number of coefficients */
                nknots > order ? nknots - order : 0L;

private:
  unsigned mutable curs,		/* current position in knots vector */
               boundary;		/* must have knots[(curs) <= x < knots(curs+1) */
  vec mutable ldel = vec(ordm1), /* differences from knots on the left */
              rdel = vec(ordm1), /* differences from knots on the right */
              a    = vec(order), /* scratch array */
              wrk  = vec(order); /* working memory */

public:
  SplineBasis(const unsigned order);
  SplineBasis(const vec knots, const unsigned order = default_order);

  unsigned get_n_basis() const {
    return ncoef;
  }

  using basisMixin::operator();
  void operator()(
      vec &out, double const x, const unsigned ders = default_ders) const {
    out.zeros();
#ifdef DO_CHECKS
    if(out.n_elem != SplineBasis::get_n_basis())
      throw_invalid_out(
        "splineBasis", out.n_elem, SplineBasis::get_n_basis());
#endif

    set_cursor(x);
    unsigned io;
    if (curs < order || (io = curs - order) > nknots) {
      /* Do nothing. x is already zero by default
       for (size_t j = 0; j < (size_t)order; j++) {
       out(j+io) = double(0); // R_NaN;
       }*/
    } else if (ders > 0) { /* slow method for derivatives */
      for(unsigned i = 0; i < (size_t)order; i++) {
        for(unsigned j = 0; j < (size_t)order; j++)
          a(j) = 0;
        a(i) = 1;
        out(i + io) = slow_evaluate(x, ders);
      }
    } else { /* fast method for value */
      basis_funcs(wrk, x);
      for (unsigned i = 0; i < wrk.n_elem; i++)
        out(i + io) = wrk(i);
    }
  }

  virtual ~SplineBasis() = default;

private:
  inline unsigned set_cursor(const double x) const {
    /* don't assume x's are sorted */
    curs = 0; /* Wall */
    boundary = 0;
    for (unsigned i = 0; i < nknots; i++) {
      if (knots[i] >= x)
        curs = i;
      if (knots[i] > x)
        break;
    }
    if (curs > ncoef) {
      if (x == knots[ncoef]){
        boundary = 1;
        curs = ncoef;
      }
    }
    return curs;
  }

  inline void diff_table(const double x, const unsigned ndiff) const {
    for (unsigned i = 0; i < ndiff; i++) {
      rdel[i] = knots[curs + i] - x;
      ldel[i] = x - knots[curs - (i + 1)];
    }
  }
  inline double slow_evaluate(const double x, unsigned nder) const {
    unsigned apt, lpt, rpt, inner;
    unsigned outer = ordm1;
    if (boundary && nder == ordm1) /* value is arbitrary */
      return 0;
    while(nder--) {  // FIXME: divides by zero
      for(inner = outer, apt = 0, lpt = curs - outer; inner--; apt++, lpt++)
        a[apt] = outer * (a[apt + 1] - a[apt]) /
          (knots[lpt + outer] - knots[lpt]);
      outer--;
    }
    diff_table(x, outer);
    while(outer--)
      for(apt = 0, lpt = outer, rpt = 0, inner = outer + 1;
          inner--; lpt--, rpt++, apt++)
        // FIXME: divides by zero
        a[apt] = (a[apt + 1] * ldel[lpt] + a[apt] * rdel[rpt]) /
          (rdel[rpt] + ldel[lpt]);
    return a[0];
  }
  /* fast evaluation of basis functions */
  inline void basis_funcs(vec &b, const double x) const {
    diff_table(x, ordm1);
    b[0] = 1;
    for (unsigned j = 1; j <= ordm1; j++) {
      double saved(0);
      for(unsigned r = 0; r < j; r++) { // do not divide by zero
        double const den = rdel[r] + ldel[j - 1 - r];
        if(den != 0) {
          double const term = b[r]/den;
          b[r] = saved + rdel[r] * term;
          saved = ldel[j - 1 - r] * term;
        } else {
          if(r != 0 || rdel[r] != 0)
            b[r] = saved;
          saved = 0;
        }
      }
      b[j] = saved;
    }
  }
};

class bs final : public SplineBasis {
public:
  vec const boundary_knots, interior_knots;
  bool const intercept;
  unsigned const df;

private:
  /* working memory */
  vec mutable wrk = arma::vec(SplineBasis::get_n_basis()),
             wrks = arma::vec(wrk.begin(), get_n_basis(), false);

public:
  bs(const vec &bk, const vec &ik,
     const bool inter = default_intercept,
     const unsigned ord = default_order);

  unsigned get_n_basis() const {
    return SplineBasis::get_n_basis() - (!intercept);
  }

  using SplineBasis::operator();
  void operator()(
      vec &out, double const x, const unsigned ders = default_ders) const {
#ifdef DO_CHECKS
    if(out.n_elem != bs::get_n_basis())
      throw_invalid_out("bs", out.n_elem, bs::get_n_basis());
#endif
    if(x < boundary_knots[0] || x > boundary_knots[1]) {
      double const k_pivot =
        x < boundary_knots[0] ?
        0.75 * boundary_knots[0] + 0.25 * knots[order] :
      0.75 * boundary_knots[1] + 0.25 * knots[knots.n_elem - order - 2],
      delta = x - k_pivot;

      auto add_term = [&](unsigned const d, double const f = 1){
        bs::operator()(wrks, k_pivot, d);
        out += f * wrks;
      };

      out.zeros();
      if (ders == 0) {
        add_term(0);
        add_term(1, delta);
        add_term(2, delta * delta/2.);
        add_term(3, delta * delta * delta /6.);

      } else if (ders == 1) {
        add_term(1);
        add_term(2, delta);
        add_term(3, delta * delta / 2.);

      } else if (ders == 2) {
        add_term(2);
        add_term(3, delta);

      } else if (ders==3) {
        add_term(3);

      }

      return;
    }

    if(intercept)
      SplineBasis::operator()(out, x, ders);
    else {
      SplineBasis::operator()(wrk, x, ders);
      for(unsigned i = 1; i < wrk.n_elem; ++i)
        out[i - 1L] = wrk[i];
    }
  }
};

class ns final : public basisMixin {
public:
  bs const bspline; // composition cf. inheritance
  bool const intercept;
  mat const q_matrix = ([&](){
    // calculate the Q matrix
    mat const_basis = bspline.basis(bspline.boundary_knots, 2);
    if (!intercept)
      const_basis = const_basis.cols(1, const_basis.n_cols - 1);
    mat qd, rd;
    if(!qr(qd, rd, const_basis.t()))
#ifdef DO_CHECKS
      throw std::invalid_argument("ns: QR decomposition failed");
#else
      { }
#endif
    inplace_trans(qd);
    return qd;
  })();
  vec const tl0, tl1, tr0, tr1;

  ns(const vec &boundary_knots, const vec &interior_knots,
     const bool intercept = default_intercept,
     const unsigned order = default_order);

  unsigned get_n_basis() const {
    return q_matrix.n_rows - 2;
  }

  using basisMixin::operator();
  void operator()(
      vec &out, double const x, const unsigned ders = default_ders) const {
#ifdef DO_CHECKS
    if(out.n_elem != ns::get_n_basis())
      throw_invalid_out("ns", out.n_elem, ns::get_n_basis());
#endif
    if(x < bspline.boundary_knots[0]) {
      if (ders==0){
        out  = tl1;
        out *= x - bspline.boundary_knots[0];
        out += tl0;

      } else if (ders == 1)
        out = tl1;
      else
        out.zeros();

      return;

    } else if (x > bspline.boundary_knots[1]) {
      if (ders==0){
        out  = tr1;
        out *= x - bspline.boundary_knots[1];
        out += tr0;

      } else if (ders==1)
        out = tr1;
      else
        out.zeros();

      return;
    }

    out = trans(bspline(x, ders));
  }

private:
  inline vec trans(const vec &x) const {
    vec out = q_matrix * (intercept ? x : x(span(1, x.n_elem - 1)));
    return out(span(2, out.size() - 1));
  }
}; // class ns

class iSpline final : public basisMixin {
public:
  bool const intercept;
  unsigned const order;
  bs const bspline; // composition cf. inheritance

private:
  vec mutable wrk = arma::vec(bspline.get_n_basis());

public:
  iSpline(const vec &boundary_knots, const vec &interior_knots,
          const bool intercept = default_intercept,
          const unsigned order = default_order);

  unsigned get_n_basis() const {
    return bspline.get_n_basis() - (!intercept);
  }

  using basisMixin::operator();
  void operator()(
      vec &out, double const x, const unsigned der = default_ders) const  {
#ifdef DO_CHECKS
    if(out.n_elem != iSpline::get_n_basis())
      throw_invalid_out("iSpline", out.n_elem, iSpline::get_n_basis());
#endif
    if(x < 0){
      out.zeros();
      return;

    }
    else if(x <= 1){
      vec &b = wrk;
      bspline(b, x, der);
      unsigned const js = bspline.interior_knots.size() > 0 ?
        static_cast<unsigned>(std::lower_bound(
          bspline.knots.begin(),
          /* TODO: should this not be end and not -1? */
          bspline.knots.end() - 1L, x) -
          bspline.knots.begin()) :
        order + 1;
      for(unsigned j = b.size(); j-- >0;)
        if(j > js)
          b[j] = 0.0;
        else if(j != b.size() - 1)
          b[j] += b[j+1];
      if(der==0)
        for(unsigned j = b.size() - 1; j-- > 0;)
          if(j + order + 1 < js)
            b[j] = 1.0;

      if(intercept)
        out = b;
      else
        out = b.subvec(1, b.n_elem - 1);
      return;

    }
    else if(der > 0)
      out.zeros();
    else
      out.fill(1);
  }
}; // class iSpline

class mSpline final : public basisMixin {
public:
  bs const bspline;
  bool const intercept;

private:
  vec mutable wrk = vec(bspline.get_n_basis());

public:
  mSpline(const vec &boundary_knots, const vec &interior_knots,
          const bool intercept = default_intercept,
          const unsigned order = default_order);

  unsigned get_n_basis() const {
    return bspline.get_n_basis() - (!intercept);
  }

  using basisMixin::operator();
  void operator()(
      vec &out, double const x, const unsigned der = default_ders) const {
#ifdef DO_CHECKS
    if(out.n_elem != mSpline::get_n_basis())
      throw_invalid_out("mSpline", out.n_elem, mSpline::get_n_basis());
#endif
    bspline(wrk, x, der);
    for (unsigned j = 0; j < bspline.get_n_basis(); j++) {
      double denom = bspline.knots(j + bspline.order) - bspline.knots(j);
      wrk(j) *= denom > 0.0 ? bspline.order / denom : 0.0;
    }

    if(intercept)
      out = wrk;
    else
      out = wrk.subvec(1, wrk.size() - 1);
  }
}; // class mSpline

class orth_poly final : public basisMixin {
public:
  vec const alpha,
            norm2,
       sqrt_norm2 = arma::sqrt(norm2);

  orth_poly(vec const &alpha, vec const &norm2);

  /**
   behaves like predict(<poly object>, newdata) except the scale is
   different. See ::get_poly_basis().
   */
  using basisMixin::operator();
  void operator()(vec &out, double const x, const unsigned ders = default_ders)
    const {
#ifdef DO_CHECKS
    if(out.n_elem != get_n_basis())
      throw std::invalid_argument("orth_poly::operator(): invalid out");
#endif
    if(ders > 0)
      throw std::runtime_error("ders > 0 is not implemented");
    out[0] = 1.;

    if(alpha.n_elem > 0L){
      out[1] = x - alpha[0];
      for(unsigned c = 1; c < alpha.n_elem; c++)
        out[c + 1L] =
          (x - alpha[c]) * out[c] - norm2[c + 1L] / norm2[c] * out[c - 1L];
    }

    out /= sqrt_norm2.subvec(1L, sqrt_norm2.n_elem - 1L);
  }

  /**
   behaves like poly(x, degree) though the output is not scaled to have unit
   norm buth rather a norm that scales like the number of samples and there
   is an intercept. I.e. similar to
      x <- rnorm(10)
      B <- poly(x, degree = 4)
      B <- cbind(1, B * sqrt(NROW(B)))
      B
      crossprod(B)
      # same as
      B <- poly(x, degree = 4)
      attr(B, "coefs")$norm2 <- attr(B, "coefs")$norm2 / NROW(B)
      cbind(1, predict(B, x))
   The orthogonal polynomial is returned by reference.
   */
  static orth_poly get_poly_basis(vec, unsigned const, mat&);

  unsigned get_n_basis() const {
    return norm2.n_elem - 1L;
  }
}; // orth_poly

} // namespace joint_bases

#endif // SPLINES_H
