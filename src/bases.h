#ifndef SPLINES_H
#define SPLINES_H

#include "arma-wrap.h"
#include <limits> // quiet_NaN
#include <stdexcept> // invalid_argument
#include "VA-joint-config.h"
#include <memory.h>
#include "wmem.h"

namespace joint_bases {

constexpr vajoint_uint default_order{4},
                       default_ders{0};
constexpr bool default_intercept{false};

using namespace arma;

/// base class for basis expansions
class basisMixin {
public:
  /// the required working memory
  virtual size_t n_wmem() const = 0;

  /// the number of basis functions
  virtual vajoint_uint n_basis() const = 0;

  /**
   * fills a vector with the (possibly derivatives) of the basis expansions
   * evaluated at x */
  void operator()
    (vec &out, double *wk_mem, double const x,
     const vajoint_uint ders = default_ders) const {
    (*this)(out.memptr(), wk_mem, x, ders);
  }
  void operator()
    (vec &out, double const x, const vajoint_uint ders = default_ders) const {
    (*this)(out.memptr(), wmem::get_double_mem(n_wmem()), x, ders);
  }
  /// returns an allocated vector
  vec operator()
    (double const x, double * wk_mem,
     vajoint_uint const ders = default_ders) const {
    vec out(n_basis());
    (*this)(out.begin(), wk_mem, x, ders);
    return out;
  }
  vec operator()
    (double const x, vajoint_uint const ders = default_ders) const {
    return (*this)(x, wmem::get_double_mem(n_wmem()), ders);
  }
  /// same as the other operator() calls but filling the out
  virtual void operator()
    (double *out, double *wk_mem, double const x,
     vajoint_uint const ders = default_ders) const = 0;
  void operator()
    (double *out, double const x,
     vajoint_uint const ders = default_ders) const {
    (*this)(out, wmem::get_double_mem(n_wmem()), x, ders);
  }

  mat basis
    (const vec &x, double *wk_mem, const vajoint_uint ders = default_ders,
     const double centre = std::numeric_limits<double>::quiet_NaN())
    const  {
    vajoint_uint const n_basis_v(n_basis()),
    n_x    (x.n_elem);
    rowvec centering =
      (std::isnan(centre) || ders > 0 ?
      zeros(n_basis_v) : operator()(centre, wk_mem, 0)).t();

    mat out(n_x, n_basis_v);
    vec wrk(n_basis_v);
    for (vajoint_uint i = 0; i < n_x; i++){
      (*this)(wrk, wk_mem, x[i], ders);
      out.row(i) = wrk.t() - centering;
    }

    return out;
  }

  virtual ~basisMixin() = default;

  virtual std::unique_ptr<basisMixin> clone() const = 0;
};

class SplineBasis : public basisMixin {
public:
  vajoint_uint const order = default_order, /* order of the spline */
                     ordm1 = order - 1;     /* order - 1 (3 for cubic splines) */
  vec const knots;	               /* knot vector */
  vajoint_uint const nknots = knots.n_elem, /* number of knots
                                           except for the boundary case */
                      ncoef =               /* number of coefficients */
                         nknots > order ? nknots - order : 0L;

  SplineBasis(const vajoint_uint order);
  SplineBasis(const vec knots, const vajoint_uint order = default_order);

  vajoint_uint n_basis() const override {
    return ncoef;
  }

  size_t n_wmem() const override {
    return 2 * ordm1 + 2 * order;
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     const vajoint_uint ders = default_ders)
    const override {
    // setup the object we need
    double * const ldel{wk_mem},       /* differences from knots on the left */
           * const rdel{ldel + ordm1}, /* differences from knots on the right */
           * const a   {rdel + ordm1}, /* scratch array */
           * const wrk {a + order};    /* working memory */
    vajoint_uint curs{}, boundary{};

    // the function we need to perform the computation
    auto set_cursor = [&](const double x){
      /* don't assume x's are sorted */
      curs = 0; /* Wall */
      boundary = 0;
      for (vajoint_uint i = 0; i < nknots; i++) {
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
    };

    auto diff_table = [&](const double x, const vajoint_uint ndiff){
      for (vajoint_uint i = 0; i < ndiff; i++) {
        rdel[i] = knots[curs + i] - x;
        ldel[i] = x - knots[curs - (i + 1)];
      }
    };

    auto slow_evaluate = [&](const double x, vajoint_uint nder){
      vajoint_uint apt, lpt, rpt, inner;
      vajoint_uint outer = ordm1;
      if (boundary && nder == ordm1) /* value is arbitrary */
        return 0.;
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
    };

    // fast evaluation of basis functions
    auto basis_funcs = [&](double *b, const double x){
      diff_table(x, ordm1);
      b[0] = 1;
      for (vajoint_uint j = 1; j <= ordm1; j++) {
        double saved(0);
        for(vajoint_uint r = 0; r < j; r++) { // do not divide by zero
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
    };

    // evaluate the basis
    std::fill(out, out + SplineBasis::n_basis(), 0);

    set_cursor(x);
    vajoint_uint io;
    if (curs < order || (io = curs - order) > nknots) {
      /* Do nothing. x is already zero by default
       for (size_t j = 0; j < (size_t)order; j++) {
       out(j+io) = double(0); // R_NaN;
       }*/
    } else if (ders > 0) { /* slow method for derivatives */
      for(vajoint_uint i = 0; i < (size_t)order; i++) {
        for(vajoint_uint j = 0; j < (size_t)order; j++)
          a[j] = 0;
        a[i] = 1;
        out[i + io] = slow_evaluate(x, ders);
      }
    } else { /* fast method for value */
      basis_funcs(wrk, x);
      for (vajoint_uint i = 0; i < order; i++)
        out[i + io] = wrk[i];
    }
  }

  virtual ~SplineBasis() = default;

  std::unique_ptr<basisMixin> clone() const override {
    return std::make_unique<SplineBasis>(*this);
  }
};

class bs final : public SplineBasis {
public:
  vec const boundary_knots, interior_knots;
  bool const intercept;
  vajoint_uint const df;

public:
  size_t n_wmem() const {
    return
      2 * std::max(SplineBasis::n_basis(), bs::n_basis()) +
        SplineBasis::n_wmem();
  }

  bs(const vec &bk, const vec &ik,
     const bool inter = default_intercept,
     const vajoint_uint ord = default_order);

  vajoint_uint n_basis() const {
    return SplineBasis::n_basis() - (!intercept);
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<bs>(*this);
  }

  using SplineBasis::operator();
  void operator()
      (double *out, double *wk_mem, double const x,
       const vajoint_uint ders = default_ders) const {
    double * const my_wk_mem{wk_mem};
    wk_mem += std::max(SplineBasis::n_basis(), bs::n_basis());

    if(x < boundary_knots[0] || x > boundary_knots[1]) {
      double const k_pivot =
        x < boundary_knots[0]
          ? 0.75 * boundary_knots[0] + 0.25 * knots[order]
          : 0.75 * boundary_knots[1] + 0.25 * knots[knots.n_elem - order - 2],
                     delta = x - k_pivot;

      auto add_term = [&](vajoint_uint const d, double const f = 1){
        bs::operator()(my_wk_mem, wk_mem, k_pivot, d);
        for(vajoint_uint i = 0; i < bs::n_basis(); ++i)
          out[i] += f * my_wk_mem[i];
      };

      std::fill(out, out + bs::n_basis(), 0);
      if (ders == 0) {
        add_term(0);
        add_term(1, delta);
        add_term(2, delta * delta / 2);
        add_term(3, delta * delta * delta / 6);

      } else if (ders == 1) {
        add_term(1);
        add_term(2, delta);
        add_term(3, delta * delta / 2.);

      } else if (ders == 2) {
        add_term(2);
        add_term(3, delta);

      } else if (ders==3)
        add_term(3);

      return;
    }

    if(intercept)
      SplineBasis::operator()(out, wk_mem, x, ders);
    else {
      SplineBasis::operator()(my_wk_mem, wk_mem, x, ders);
      for(vajoint_uint i = 1; i < SplineBasis::n_basis(); ++i)
        out[i - 1L] = my_wk_mem[i];
    }
  }
};

class ns final : public basisMixin {
public:
  bs const bspline; // composition cf. inheritance
  bool const intercept;
  mat const q_matrix = ([&](){
    // calculate the Q matrix
    mat const_basis = bspline.basis
      (bspline.boundary_knots, wmem::get_double_mem(bspline.n_wmem()), 2);
    if (!intercept)
      const_basis = const_basis.cols(1, const_basis.n_cols - 1);
    mat qd, rd;
    if(!qr(qd, rd, const_basis.t()))
      throw std::invalid_argument("ns: QR decomposition failed");
    inplace_trans(qd);
    return qd;
  })();
  vec const tl0, tl1, tr0, tr1;

  ns(const vec &boundary_knots, const vec &interior_knots,
     const bool intercept = default_intercept,
     const vajoint_uint order = default_order);

  size_t n_wmem() const {
    return bspline.n_wmem();
  }

  vajoint_uint n_basis() const {
    return q_matrix.n_rows - 2;
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<ns>(*this);
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     vajoint_uint const ders = default_ders) const {
    if(x < bspline.boundary_knots[0]) {
      if (ders==0){
        for(vajoint_uint i = 0; i < ns::n_basis(); ++i){
          out[i] = tl1[i];
          out[i] *= x - bspline.boundary_knots[0];
          out[i] += tl0[i];
        }

      } else if (ders == 1)
        std::copy(tl1.begin(), tl1.end(), out);
      else
        std::fill(out, out + ns::n_basis(), 0);

      return;

    } else if (x > bspline.boundary_knots[1]) {
      if (ders==0){
        for(vajoint_uint i = 0; i < ns::n_basis(); ++i){
          out[i] = tr1[i];
          out[i] *= x - bspline.boundary_knots[1];
          out[i] += tr0[i];
        }

      } else if (ders==1)
        std::copy(tr1.begin(), tr1.end(), out);
      else
        std::fill(out, out + ns::n_basis(), 0);

      return;
    }

    // TODO: very inefficient
    arma::vec res{trans(bspline(x, ders))};
    std::copy(res.begin(), res.end(), out);
  }

private:
  vec trans(const vec &x) const {
    // TODO: very inefficient
    vec out = q_matrix * (intercept ? x : x(span(1, x.n_elem - 1)));
    return out(span(2, out.size() - 1));
  }
}; // class ns

class iSpline final : public basisMixin {
public:
  bool const intercept;
  vajoint_uint const order;
  bs const bspline; // composition cf. inheritance

public:
  iSpline(const vec &boundary_knots, const vec &interior_knots,
          const bool intercept = default_intercept,
          const vajoint_uint order = default_order);

  size_t n_wmem() const {
    return bspline.n_wmem() + bspline.n_basis();
  }

  vajoint_uint n_basis() const {
    return bspline.n_basis() - (!intercept);
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<iSpline>(*this);
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     vajoint_uint const ders = default_ders) const  {
    double * const b{wk_mem};
    vajoint_uint const n_b{bspline.n_basis()};
    wk_mem += n_b;

    if(x < 0){
      std::fill(out, out + n_basis(), 0);
      return;

    }
    else if(x <= 1){
      bspline(b, wk_mem, x, ders);
      vajoint_uint const js = bspline.interior_knots.size() > 0 ?
        static_cast<vajoint_uint>(std::lower_bound(
          bspline.knots.begin(),
          /* TODO: should this not be end and not -1? */
          bspline.knots.end() - 1L, x) -
          bspline.knots.begin()) :
        order + 1;
      for(vajoint_uint j = n_b; j-- >0;)
        if(j > js)
          b[j] = 0.0;
        else if(j != n_b - 1)
          b[j] += b[j+1];
      if(ders==0)
        for(vajoint_uint j = n_b - 1; j-- > 0;)
          if(j + order + 1 < js)
            b[j] = 1.0;

      if(intercept)
        std::copy(b, b + bspline.n_basis(), out);
      else
        std::copy(b + 1, b + bspline.n_basis(), out);
      return;

    }
    else if(ders > 0)
      std::fill(out, out + n_basis(), 0);
    else
      std::fill(out, out + n_basis(), 1);
  }
}; // class iSpline

class mSpline final : public basisMixin {
public:
  bs const bspline;
  bool const intercept;

public:
  mSpline(const vec &boundary_knots, const vec &interior_knots,
          const bool intercept = default_intercept,
          const vajoint_uint order = default_order);

  size_t n_wmem() const {
    return bspline.n_wmem() + bspline.n_basis();
  }

  vajoint_uint n_basis() const {
    return bspline.n_basis() - (!intercept);
  }

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<mSpline>(*this);
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     vajoint_uint const ders = default_ders) const {
    double * const wrk{wk_mem};
    wk_mem += bspline.n_basis();

    bspline(wrk, wk_mem, x, ders);
    for (vajoint_uint j = 0; j < bspline.n_basis(); ++j) {
      double denom = bspline.knots(j + bspline.order) - bspline.knots(j);
      wrk[j] *= denom > 0.0 ? bspline.order / denom : 0.0;
    }

    if(intercept)
      std::copy(wrk, wrk + bspline.n_basis(), out);
    else
      std::copy(wrk + 1, wrk + bspline.n_basis(), out);
  }
}; // class mSpline

class orth_poly final : public basisMixin {
  // coefficients for the orthogonal polynomial
  vec alpha,
      norm2,
      sqrt_norm2{arma::sqrt(norm2)};
  // flags for whether a raw polynomial is used and whether there is a intercept
  bool raw,
       intercept;
  // the number of basis function plus the possible intercept
  vajoint_uint n_basis_v;

public:
  // constructor for raw == true
  orth_poly
    (vajoint_uint const degree, bool const intercept = default_intercept);

  // constructor for raw == false
  orth_poly(vec const &alpha, vec const &norm2,
            bool const intercept = default_intercept);

  std::unique_ptr<basisMixin> clone() const {
    return std::make_unique<orth_poly>(*this);
  }

  size_t n_wmem() const {
    return 0;
  }

  /**
   * behaves like predict(<poly object>, newdata) except there might be an
   * intercept */
  using basisMixin::operator();
  void operator()(double *out, double *wk_mem, double const x,
                  vajoint_uint const ders = default_ders)
    const {
    if(ders > 0)
      throw std::runtime_error("ders > 0 is not implemented");

    if(raw){
      if(intercept){
        out[0] = 1.;
        for(vajoint_uint c = 1; c < n_basis_v; c++)
          out[c] = out[c - 1] * x;

      } else {
        double val{1};
        for(vajoint_uint c = 0; c < n_basis_v; c++)
          out[c] = val *= x;
      }

      return;
    }

    out[0] = 1.;

    if(alpha.n_elem > 0L){
      out[1] = x - alpha[0];
      for(vajoint_uint c = 1; c < alpha.n_elem; c++)
        out[c + 1L] =
          (x - alpha[c]) * out[c] - norm2[c + 1L] / norm2[c] * out[c - 1L];
    }

    for(vajoint_uint j = 1; j < alpha.n_elem + 1; ++j)
      out[j] /= sqrt_norm2.at(j + 1);
  }

  /**
   * behaves like poly(x, degree). The orthogonal polynomial is returned by
   * reference.
   */
  static orth_poly poly_basis(vec, vajoint_uint const, mat&);

  vajoint_uint n_basis() const {
    return n_basis_v;
  }
}; // orth_poly


using bases_vector = std::vector<std::unique_ptr<basisMixin> >;

/// simple function clone bases
inline bases_vector clone_bases(const bases_vector &bases){
  std::vector<std::unique_ptr<basisMixin> > out;
  out.reserve(bases.size());
  for(auto &x : bases)
    out.emplace_back(x->clone());
  return out;
}

} // namespace joint_bases

#endif // SPLINES_H
