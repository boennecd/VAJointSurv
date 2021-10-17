#ifndef SPLINES_H
#define SPLINES_H

#include "arma-wrap.h"
#include <limits> // quiet_NaN
#include <stdexcept> // invalid_argument
#include "VA-joint-config.h"
#include <memory.h>
#include "wmem.h"
#include "lp-joint.h"
#include <algorithm>

namespace joint_bases {

constexpr vajoint_uint default_order{4};
constexpr int default_ders{0};
constexpr bool default_intercept{false};

using namespace arma;

/// base class for basis expansions
class basisMixin {
public:
  /// the required working memory
  virtual size_t n_wmem() const = 0;

  /// the number of basis functions
  virtual vajoint_uint n_basis() const = 0;

  /// lower limit of integrals
  double lower_limit{};

  /**
   * fills a vector with the (possibly derivatives) of the basis expansions
   * evaluated at x */
  void operator()
    (vec &out, double *wk_mem, double const x,
     const int ders = default_ders) const {
    (*this)(out.memptr(), wk_mem, x, ders);
  }
  /// returns an allocated vector
  vec operator()
    (double const x, double * wk_mem,
     int const ders = default_ders) const {
    vec out(n_basis());
    (*this)(out.begin(), wk_mem, x, ders);
    return out;
  }
  /// same as the other operator() calls but filling the out
  virtual void operator()
    (double *out, double *wk_mem, double const x,
     int const ders = default_ders) const = 0;

  mat basis
    (const vec &x, double *wk_mem, const int ders = default_ders,
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

  SplineBasis(const vec &knots, const vajoint_uint order = default_order,
              bool const with_intercept = true);

  SplineBasis(SplineBasis const &other):
    SplineBasis(other.knots, other.order,
                static_cast<bool>(other.integral_basis)) { }

  vajoint_uint n_basis() const override {
    return ncoef;
  }

  size_t n_wmem() const override {
    return n_wmem_v;
  }

  using basisMixin::operator();
  void operator()
    (double *out, double *wk_mem, double const x,
     const int ders = default_ders)
    const override {
    if(ders < 0){
      if(ders < -1)
        throw std::runtime_error("not implemented for ders < -1");
      // use formulas from Monotone Regression Splines in Action
      //    https://doi.org/10.1214/ss/1177012761
      // TODO: can be implemented smarter...
      double * const basis_mem = wk_mem;
      wk_mem += integral_basis->n_basis();

      double dorder{static_cast<double>(order)};

      // computes the indefinte integral at the upper or lower limit. The
      // function must first be called at the upper limit
      auto add_int = [&](double const lim, bool const is_upper){
        // evaluate the basis which is one order greater
        (*integral_basis)(basis_mem, wk_mem, lim, ders + 1);

        // find the index j such that knots[j] <= lim < knots[j + 1]
        // use -1 if x < knots[0]
        auto const idx_knot_start =
          ([&]{
            auto knot_j = std::upper_bound(knots.begin(), knots.end(), lim);
            return std::distance(knots.begin(), knot_j) - 1;
          })();

        // x is too small for these basis function to be active
        vajoint_uint const idx_no_support
            {static_cast<vajoint_uint>(idx_knot_start + 1)};
        if(is_upper)
          std::fill(out + idx_no_support, out + ncoef, 0);
        // x is large enough that we have integrated over the full support of
        // these basis functions
        vajoint_uint i{};
        vajoint_uint const is_capped =
          idx_knot_start + 1 >= order ? idx_knot_start + 1 - order : 0;
        for(; i < is_capped; ++i)
          if(is_upper)
            out[i]  = (knots[i + order] - knots[i]) / dorder;
          else
            out[i] -= (knots[i + order] - knots[i]) / dorder;

        // the residual is somewhere in between
        for(; i < idx_no_support; ++i){
          double su{};
          for(vajoint_uint j = i; j < idx_no_support; ++j)
            // TODO: redundant computations
            su += basis_mem[j];
          if(is_upper)
            out[i]  = su * (knots[i + order] - knots[i]) / dorder;
          else
            out[i] -= su * (knots[i + order] - knots[i]) / dorder;
        }
      };

      add_int(x, true); // the upper limit
      if(lower_limit > knots[0])
        add_int(lower_limit, false); // the lower limit

      return;
    }

    if(ders == 0){
      comp_basis(x, out, order, wk_mem);
      return;
    }

    // setup the object we need
    double * const ldel{wk_mem},       /* differences from knots on the left */
           * const rdel{ldel + ordm1}, /* differences from knots on the right */
           * const a   {rdel + ordm1}; /* scratch array */
    vajoint_uint curs{}, boundary{};

    // the function we need to perform the computation
    auto set_cursor = [&](const double x){
      curs = 0; /* Wall */
      boundary = 0;
      while(curs < nknots && knots[curs] <= x)
        ++curs;
      if(curs > ncoef && x == knots[ncoef]){
        boundary = 1;
        curs = ncoef;
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

    // evaluate the basis
    // TODO: need to replace
    std::fill(out, out + SplineBasis::n_basis(), 0);
    set_cursor(x);
    vajoint_uint const io{curs < order ? 0 : curs - order};

    vajoint_uint const uders{static_cast<vajoint_uint>(ders)};
    for(vajoint_uint i = 0; i < order; i++) {
      for(vajoint_uint j = 0; j < order; j++)
        a[j] = 0;
      a[i] = 1;
      out[i + io] = slow_evaluate(x, uders);
    }
  }

  void comp_basis(double const x, double *out, unsigned const ord_p1,
                  double * wk_mem) const {
    // De Boor's algorithm
    //    https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
    unsigned const ord{ord_p1 - 1};

    // find the knot such that knot[i] <= x < knot[i + 1]
    auto it_inter = std::upper_bound(knots.begin(), knots.end(), x);
    unsigned const n_knots{knots.size()},
                n_knots_m1{n_knots - 1},
                   dim_out{n_knots - ord_p1};

    // deal with the matching boundaries
    if(it_inter == knots.end())
      while(it_inter != knots.begin() && *std::prev(it_inter) == x)
        --it_inter;

    std::fill(out, out + dim_out, 0);
    if(it_inter == knots.begin() ||
       (it_inter == knots.end() && *std::prev(it_inter) < x))
      return;

    --it_inter;
    unsigned const idx_greater = std::distance(knots.begin(), it_inter);
    double * const D{wk_mem};
    std::fill(D, D + ord_p1 * ord_p1, 0);

    for(unsigned j = idx_greater + 1 >= ord ? 0 : ord - idx_greater - 1;
        j < ord_p1; ++j)
      D[j + j * ord_p1] = 1;

    for(unsigned r = 0; r < ord; ++r){
      for(unsigned j = std::min(ord, r + n_knots_m1 - idx_greater);
          j-- > r && j + idx_greater + 1 >= ord; ){
        double const k1{knots[j + 1 + idx_greater - r]},
                     k2{knots[j     + idx_greater - ord + 1]};

        double const alpha{k1 != k2 ? (x - k2) / (k1 - k2) : 0},
          diff_alpha{1 - alpha};

        unsigned const end_update{std::min(j + 2, ord_p1)};
        unsigned i{end_update > r + 1 ? end_update - r - 2 : 0};

        for(; i < end_update; ++i){
          D[i + (j + 1) * ord_p1] *= alpha;
          D[i + (j + 1) * ord_p1] += diff_alpha * D[i + j * ord_p1];
        }
      }
    }

    if(idx_greater < ord){
      unsigned const shift_D{ord - idx_greater};
      std::copy(D + shift_D + ord * ord_p1, D + ord_p1 + ord * ord_p1, out);
    } else {
      unsigned const shift_out{idx_greater - ord},
      n_cp{std::min(dim_out - shift_out, ord_p1)};
      std::copy(D + ord * ord_p1, D + n_cp + ord * ord_p1, out + shift_out);
    }
  }

  virtual ~SplineBasis() = default;

  std::unique_ptr<basisMixin> clone() const override {
    return std::make_unique<SplineBasis>(*this);
  }

private:
  std::unique_ptr<SplineBasis> integral_basis;
  size_t const n_wmem_v
    {
    integral_basis
      ? integral_basis->n_wmem() + integral_basis->n_basis()
        : std::max(2 * ordm1 + 2 * order, order * order)
    };
};

class bs final : public SplineBasis {
public:
  vec const boundary_knots, interior_knots;
  bool const intercept;
  vajoint_uint const df;
  size_t n_wmem_v
    {2 * std::max(SplineBasis::n_basis(), bs::n_basis()) +
      SplineBasis::n_wmem()};

public:
  size_t n_wmem() const {
    return n_wmem_v;
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
       const int ders = default_ders) const {
    double * const my_wk_mem{wk_mem};
    wk_mem += std::max(SplineBasis::n_basis(), bs::n_basis());

    if(x < boundary_knots[0] || x > boundary_knots[1]) {
      double const k_pivot =
        x < boundary_knots[0]
          ? 0.75 * boundary_knots[0] + 0.25 * knots[order]
          : 0.75 * boundary_knots[1] + 0.25 * knots[knots.n_elem - order - 2],
                     delta = x - k_pivot;

      auto add_term = [&](int const d, double const f = 1){
        bs::operator()(my_wk_mem, wk_mem, k_pivot, d);
        for(vajoint_uint i = 0; i < bs::n_basis(); ++i)
          out[i] += f * my_wk_mem[i];
      };

      std::fill(out, out + bs::n_basis(), 0);
      if(ders < 0 )
        throw std::invalid_argument("integration outside of bounds is not implemented with bs");
      else if(ders == 0){
        add_term(0);
        add_term(1, delta);
        add_term(2, delta * delta / 2);
        add_term(3, delta * delta * delta / 6);

      } else if(ders == 1){
        add_term(1);
        add_term(2, delta);
        add_term(3, delta * delta / 2.);

      } else if(ders == 2){
        add_term(2);
        add_term(3, delta);
      } else if(ders == 3)
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
    return bspline.n_wmem() + q_matrix.n_rows + bspline.n_basis();
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
     int const ders = default_ders) const {
    if(x < bspline.boundary_knots[0]) {
      if(ders==0){
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

    double * const lhs = wk_mem;
    wk_mem += q_matrix.n_rows;
    double * const b = wk_mem;
    wk_mem += bspline.n_basis();
    bspline(b, wk_mem, x, ders);

    std::fill(lhs, lhs + q_matrix.n_rows, 0);
    lp_joint::mat_vec
      (lhs, q_matrix.begin(), b + (!intercept), q_matrix.n_rows,
       q_matrix.n_cols);

    std::copy(lhs + 2, lhs + q_matrix.n_rows, out);
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
     int const ders = default_ders) const  {
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
     int const ders = default_ders) const {
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
  // the matrix to map from the raw polynomial to the orthogonal polynomial
  // see https://stats.stackexchange.com/a/472289/81865
  std::vector<double> orth_map;

  // evaluates the polynomial with raw == TRUE
  static void eval_raw
    (double *out, double const x, bool const inter, int const ders,
     vajoint_uint const degree, double const lb) {
    vajoint_uint const dim{degree + inter};

    if(ders == 0){
      double val{inter ? 1 : x};
      for(vajoint_uint c = 0; c < dim; ++c, val *= x)
        out[c] = val;
      return;

    } else if(ders < 0){
      // compute the starting value
      double val_upper{x},
             val_lower{lb};
      vajoint_uint const uders{static_cast<vajoint_uint>(-ders)};
      for(vajoint_uint i = 2; i <= uders; ++i){
        val_upper *= x  / static_cast<double>(i);
        val_lower *= lb / static_cast<double>(i);
      }

      if(!inter){
        val_upper *= x  / (uders + 1);
        val_lower *= lb / (uders + 1);
      }

      for(vajoint_uint c = 0; c < dim; c++){
        out[c] = val_upper - val_lower;
        val_upper *= x  / static_cast<double>(c + uders + 1 + !inter);
        val_lower *= lb / static_cast<double>(c + uders + 1 + !inter);
        if(c + 1 + !inter >= uders){
          val_upper *= c + 1. + !inter;
          val_lower *= c + 1. + !inter;
        }
      }

    } else { // ders > 0
      vajoint_uint const uders{static_cast<vajoint_uint>(ders)};

      if(inter){
        std::fill(out, out + uders, 0);
        double val{1};
        for(vajoint_uint c = uders; c < dim; c++){
          vajoint_uint mult{c};
          for(vajoint_uint cc = c; --cc > c - uders;)
            mult *= cc;
          out[c] = mult * val;
          val *= x;
        }

      } else {
        std::fill(out, out + uders - 1, 0);
        double val{1};
        for(vajoint_uint c = uders - 1; c < dim; c++){
          vajoint_uint mult{c + 1};
          for(vajoint_uint cc = c + 1; --cc > c - uders + 1;)
            mult *= cc;
          out[c] = mult * val;
          val *= x;
        }
      }
    }
  }

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
    return intercept ? n_basis_v : n_basis_v + 1;
  }

  /**
   * behaves like predict(<poly object>, newdata) except there might be an
   * intercept */
  using basisMixin::operator();
  void operator()(double *out, double *wk_mem, double const x,
                  int const ders = default_ders)
    const {

    if(raw){
      eval_raw(out, x, intercept, ders, n_basis_v - intercept, lower_limit);
      return;
    }

    if(ders == 0){
      out[0] = 1.;
      double old{1};

      if(alpha.n_elem > 0L){
        out[0 + intercept] = x - alpha[0];
        for(vajoint_uint c = 1; c < alpha.n_elem; c++){
          out[c + intercept] =
            (x - alpha[c]) * out[c - 1 + intercept] - norm2[c + 1] /
            norm2[c] * old;
          old = out[c - 1 + intercept];
        }
      }

      for(vajoint_uint j = 1; j < alpha.n_elem + 1; ++j)
        out[j - 1 + intercept] /= sqrt_norm2.at(j + 1);

      return;
    }

    // compute the raw polynomial and multiply on the matrix
    // TODO: can likely be done in more stable way?
    eval_raw(wk_mem, x, true, ders, n_basis_v - intercept, lower_limit);

    std::fill(out, out + n_basis_v, 0);
    // handle the intercept term
    auto g = orth_map.begin() + !intercept;
    if(intercept)
      out[0] = *g++ * wk_mem[0];

    // handle the other terms
    for(vajoint_uint j = 1; j < alpha.size() + 1; ++j)
      for(vajoint_uint i = 0; i <= j; ++i)
        out[j - !intercept] += wk_mem[i] * *g++;
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
