#include "arma-eigen-wrap.h"
#include "marker-term.h"
#include "psqn.h"
#include "psqn-reporter.h"
#include "log-cholesky.h"
#include "kl-term.h"
#include "cfaad/AAD.h"
#include <numeric>
#include <limits>
#include "survival-term.h"
#include <array>
#include "prof-vajoint.h"
#include "ghq-delayed-entry.h"

using Rcpp::List;
using Rcpp::NumericMatrix;
using Rcpp::NumericVector;

namespace {
/// function to perform max(a1, a2, a3, ...)
template<class T>
T many_max(T a){
  return a;
}

template<class T, class ... Args>
T many_max(T a, Args ... args){
 return std::max<T>(many_max(args...), a);
}

/**
 * takes a functor that returns a double and returns NaN if the functor returns
 * throws any error.
 */
template<class T>
double nan_if_fail(T x){
  try {
    return x();
  } catch (...) {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

/// tapes used for the Number class
std::vector<cfaad::Tape> number_tapes;
/// sets the number of tapes
void set_num_tapes(vajoint_uint const n_threads){
  number_tapes.resize(n_threads);
}
/// sets the tape for this thread
void set_my_tape(){
#ifdef _OPENMP
  if(get_thread_num())
    cfaad::Number::tape = &number_tapes[get_thread_num()];
#endif
}
/// the currently set quadrature rule to use
survival::node_weight const * cur_quad_rule;
/// the currently set Gauss-Hermite quadrature rule to use
ghqCpp::ghq_data const * cur_gh_quad_rule;

survival::node_weight node_weight_from_list(List dat){
  NumericVector nodes = dat["node"],
              weigths = dat["weight"];
  if(nodes.size() != weigths.size())
    throw std::runtime_error("nodes.size() != weigths.size()");

  return { &nodes[0], &weigths[0], static_cast<vajoint_uint>(nodes.size()) };
}

ghqCpp::ghq_data gh_node_weight_from_list(List dat){
  NumericVector nodes = dat["node"],
              weigths = dat["weight"];
  if(nodes.size() != weigths.size())
    throw std::runtime_error("nodes.size() != weigths.size()");

  return { &nodes[0], &weigths[0], static_cast<vajoint_uint>(nodes.size()) };
}

} // namespace

using cfaad::Number;

/// creates a poly term
template<class basisT>
std::unique_ptr<joint_bases::basisMixin> poly_term_from_list(List dat){
  if(!Rf_inherits(dat, "poly_term"))
    throw std::runtime_error("wrong class of term was passed");

  List coefs = dat["coefs"];
  arma::vec alpha{Rcpp::as<arma::vec>(coefs["alpha"])},
            norm2{Rcpp::as<arma::vec>(coefs["norm2"])};
  bool const raw{Rcpp::as<bool>(dat["raw"])},
       intercept{Rcpp::as<bool>(dat["intercept"])},
         use_log{Rcpp::as<bool>(dat["use_log"])};

  return raw
    ? std::make_unique<basisT>(alpha.size(), intercept, use_log)
    : std::make_unique<basisT>(alpha, norm2, intercept, use_log);
}

/// creates a bs term
template<class basisT>
std::unique_ptr<joint_bases::basisMixin> bs_term_from_list(List dat){
  if(!Rf_inherits(dat, "bs_term"))
    throw std::runtime_error("wrong class of term was passed");

  arma::vec i_knots{Rcpp::as<arma::vec>(dat["knots"])},
            b_knots{Rcpp::as<arma::vec>(dat["Boundary.knots"])};

  bool const intercept{Rcpp::as<bool>(dat["intercept"])},
               use_log{Rcpp::as<bool>(dat["use_log"])};
  vajoint_uint const degree{Rcpp::as<vajoint_uint>(dat["degree"])};

  return std::make_unique<basisT>
    (b_knots, i_knots, intercept, degree + 1, use_log);
}

/// creates a ns term
template<class basisT>
std::unique_ptr<joint_bases::basisMixin> ns_term_from_list(List dat){
  if(!Rf_inherits(dat, "ns_term"))
    throw std::runtime_error("wrong class of term was passed");

  arma::vec i_knots{Rcpp::as<arma::vec>(dat["knots"])},
            b_knots{Rcpp::as<arma::vec>(dat["Boundary.knots"])};

  bool const intercept{Rcpp::as<bool>(dat["intercept"])},
               use_log{Rcpp::as<bool>(dat["use_log"])};
  vajoint_uint const degree{Rcpp::as<vajoint_uint>(dat["degree"])};

  return std::make_unique<basisT>
    (b_knots, i_knots, intercept, degree + 1, use_log);
}

/// creates a stacked term
std::unique_ptr<joint_bases::basisMixin> basis_from_list(List dat);

template<class basisT>
std::unique_ptr<joint_bases::basisMixin> stacked_term_from_list(List dat){
  if(!Rf_inherits(dat, "stacked_term"))
    throw std::runtime_error("wrong class of term was passed");

  List terms = dat["terms"];

  if(terms.size() < 1)
    throw std::invalid_argument("stacked_term without terms");

  joint_bases::bases_vector bases;
  bases.reserve(terms.size());
  for(auto sexp_term : terms){
    bases.emplace_back(basis_from_list(List(sexp_term)));
  }

  return std::make_unique<basisT>(bases);
}

/// returns a pointer to an expansions given an R List with data.
std::unique_ptr<joint_bases::basisMixin> basis_from_list(List dat){
  if(Rf_inherits(dat, "weighted_term")){
    List term = dat["term"];
    if(Rf_inherits(term, "weighted_term")){
      throw std::invalid_argument("weighted_term of weighted_term is not supported");

    } else if(Rf_inherits(term, "poly_term")){
      return
        poly_term_from_list
        <joint_bases::weighted_basis<joint_bases::orth_poly> >(term);

    } else if(Rf_inherits(term, "bs_term")){
      return
        bs_term_from_list
        <joint_bases::weighted_basis<joint_bases::bs> >(term);

    } else if(Rf_inherits(term, "ns_term")){
      return
        ns_term_from_list
        <joint_bases::weighted_basis<joint_bases::ns> >(term);

    } else if(Rf_inherits(term, "stacked_term")){
      return
        stacked_term_from_list
        <joint_bases::weighted_basis<joint_bases::stacked_basis> >(term);

    }

  } else if(Rf_inherits(dat, "poly_term")){
    return poly_term_from_list<joint_bases::orth_poly>(dat);

  } else if(Rf_inherits(dat, "bs_term")){
    return bs_term_from_list<joint_bases::bs>(dat);

  } else if(Rf_inherits(dat, "ns_term")){
    return ns_term_from_list<joint_bases::ns>(dat);

  } else if(Rf_inherits(dat, "stacked_term")){
    return stacked_term_from_list<joint_bases::stacked_basis>(dat);

  }

  throw std::invalid_argument("expansions is not implemented");
  return std::make_unique<joint_bases::orth_poly>(2, false);
}

/// returns a smart pointer to basisMixin object
// [[Rcpp::export(rng = false)]]
SEXP expansion_object(List dat){
  auto res = basis_from_list(dat);
  return Rcpp::XPtr<joint_bases::basisMixin>{res.release()};
}

/// evaluates an expansion at x with using an XPtr
// [[Rcpp::export(rng = false)]]
NumericMatrix eval_expansion
  (SEXP ptr, NumericVector const x, NumericMatrix const weights, int const ders,
   double lower_limit){
  Rcpp::XPtr<joint_bases::basisMixin> basis(ptr);
  if(basis->n_weights()!=static_cast<unsigned>(weights.nrow()))
    throw std::invalid_argument("Number of weights does not match");

  NumericMatrix out(basis->n_basis(), x.size());

  std::unique_ptr<double[]> wmem(new double[basis->n_wmem()]);
  basis->set_lower_limit(lower_limit);
  for(R_len_t i = 0; i < x.size(); ++i)
    (*basis)(&out.column(i)[0], wmem.get(), x[i], &weights(0,i), ders);

  return out;
}

/// needs to be forward declared for lower_bound_caller
class lower_bound_term;

/// class to evaluate the elements function
class lower_bound_caller {
  subset_params *par_idx;
  marker::marker_dat *m_dat;
  kl_term *kl_dat;
  friend lower_bound_term;
  /**
   * the parameter vector with the full matrices. Filled during setup. It is
   * only the model parameters and not the variational parameters.
   */
  std::vector<double> par_vec;
  /// set to true if setup() failed
  bool setup_failed;

public:
  static bool optimze_survival;

  lower_bound_caller(std::vector<lower_bound_term const*> &);

  void setup(double const *val, bool const comp_grad);
  double eval_func(lower_bound_term const &obj, double const * val);
  double eval_grad(lower_bound_term const &obj, double const * val,
                   double *gr);
};

bool lower_bound_caller::optimze_survival = true;

class lower_bound_term {
  subset_params const &par_idx;
  marker::marker_dat const &m_dat;
  survival::survival_dat const &s_dat;
  kl_term const &kl_dat;
  survival::delayed_dat const &d_dat;

  std::vector<vajoint_uint> marker_indices;
  /// indices of survival outcomes stored as (index, type of outcome)
  std::vector<std::array<vajoint_uint, 2L> > surv_indices;
  size_t n_global, n_private;

  bool has_delayed_entry{false};
  size_t delayed_entry_idx;

  friend lower_bound_caller;
public:
  lower_bound_term
  (subset_params const &par_idx, marker::marker_dat const &m_dat,
   survival::survival_dat const &s_dat, kl_term const &kl_dat,
   survival::delayed_dat const &d_dat):
  par_idx(par_idx), m_dat(m_dat), s_dat(s_dat), kl_dat(kl_dat), d_dat(d_dat),
  n_global(par_idx.n_params<true>()),
  n_private{par_idx.n_va_params<true>()}
  { }

  void add_marker_index(const vajoint_uint idx){
    marker_indices.emplace_back(idx);
  }

  void add_surv_index(const vajoint_uint idx, const vajoint_uint type){
    surv_indices.emplace_back(std::array<vajoint_uint, 2L>{idx, type});
  }

  void add_delayed_entry(size_t const idx){
    has_delayed_entry = true;
    delayed_entry_idx = idx;
  }

  void shrink_to_fit() {
    marker_indices.shrink_to_fit();
    surv_indices.shrink_to_fit();
  }

  size_t global_dim() const {
    return n_global;
  }
  size_t private_dim() const {
    return n_private;
  }

  double comp(double const *p, double *gr,
              lower_bound_caller const &caller, bool const comp_grad) const {
    if(caller.setup_failed)
      return std::numeric_limits<double>::quiet_NaN();
    set_my_tape();

    vajoint_uint const n_rng{par_idx.va_mean_end() - par_idx.va_mean()};
    wmem::rewind();

    if(!comp_grad){
      // find the required working memory
      // TODO: do this once
      vajoint_uint const n_wmem =
        many_max<vajoint_uint>(log_chol::pd_mat::n_wmem(n_rng),
                               m_dat.n_wmem(),
                               kl_dat.n_wmem(),
                               s_dat.n_wmem()[0] + s_dat.n_wmem()[1]);

      double * const inter_mem{wmem::get_double_mem(n_wmem)};
      double * const par_vec{wmem::get_double_mem(par_idx.n_params_w_va())};
      auto mem_mark = wmem::mem_stack().set_mark_raii();

      // copy the global parameters
      std::copy(caller.par_vec.begin(), caller.par_vec.end(), par_vec);

      // insert the variational parameters
      std::copy(p + par_idx.va_mean<true>(), p + par_idx.va_mean_end<true>(),
                par_vec + par_idx.va_mean<false>());
      log_chol::pd_mat::get(p + par_idx.va_vcov<true>(), n_rng,
                            par_vec + par_idx.va_vcov<false>(), inter_mem);

      // compute the lower bound terms and return
      double res = kl_dat.eval(par_vec, inter_mem);
      for(vajoint_uint idx : marker_indices)
        res += m_dat(par_vec, inter_mem, idx);
      if(lower_bound_caller::optimze_survival){
        for(auto &idx : surv_indices)
          res += s_dat(par_vec, inter_mem, idx[0], idx[1],
                       inter_mem + s_dat.n_wmem()[0], *cur_quad_rule);

        if(has_delayed_entry){
          ghqCpp::simple_mem_stack<double> &my_stack = wmem::mem_stack();
          res += d_dat(par_vec, my_stack, delayed_entry_idx,
                       *cur_quad_rule, *cur_gh_quad_rule);
        }
      }

      return res;
    }

    // find the required working memory
    // TODO: compute this once
    vajoint_uint const n_wmem_num
      {many_max<vajoint_uint>(m_dat.n_wmem(), s_dat.n_wmem()[0])};
    vajoint_uint const n_wmem_dub
      {many_max<vajoint_uint>
        (log_chol:: pd_mat::n_wmem(n_rng),
         log_chol::dpd_mat::n_wmem(n_rng),
         log_chol::dpd_mat::n_wmem(par_idx.marker_info().size()),
         log_chol::dpd_mat::n_wmem(par_idx.n_shared_surv()),
         kl_dat.n_wmem(),
         s_dat.n_wmem()[1])};

    Number * const inter_mem_num{wmem::get_Number_mem(n_wmem_num)};
    double * const inter_mem_dub{wmem::get_double_mem(n_wmem_dub)};

    // setup the parameter vectors
    double * const par_vec_dub
      {wmem::get_double_mem(par_idx.n_params_w_va<false>())},
           * const par_vec_gr
      {wmem::get_double_mem(par_idx.n_params_w_va<false>())};
    Number * const par_vec_num
      {wmem::get_Number_mem(par_idx.n_params_w_va<false>())};
    auto mem_mark = wmem::mem_stack().set_mark_raii();

    Number::tape->rewind();
    std::fill(par_vec_gr, par_vec_gr + par_idx.n_params_w_va<false>(), 0);

    // copy the global parameters
    std::copy(caller.par_vec.begin(), caller.par_vec.end(), par_vec_dub);
    cfaad::convertCollection
      (caller.par_vec.begin(), caller.par_vec.end(), par_vec_num);

    // insert the variational parameters
    std::copy(p + par_idx.va_mean<true>(), p + par_idx.va_mean_end<true>(),
              par_vec_dub + par_idx.va_mean<false>());
    cfaad::convertCollection
      (p + par_idx.va_mean<true>(), p + par_idx.va_mean_end<true>(),
       par_vec_num + par_idx.va_mean<false>());

    log_chol::pd_mat::get
      (p + par_idx.va_vcov<true>(), n_rng,
       par_vec_dub + par_idx.va_vcov<false>(), inter_mem_dub);
    cfaad::convertCollection
      (par_vec_dub + par_idx.va_vcov<false>(),
       par_vec_dub + par_idx.va_vcov<false>() + n_rng * n_rng,
       par_vec_num + par_idx.va_vcov<false>());

    // compute the lower bound terms.
    Number res{0};
    for(vajoint_uint idx : marker_indices)
      res += m_dat(par_vec_num, inter_mem_num, idx);

    res += kl_dat.grad(par_vec_gr, par_vec_dub, inter_mem_dub);

    if(lower_bound_caller::optimze_survival){
      for(auto &idx : surv_indices)
        res += s_dat(par_vec_num, inter_mem_num, idx[0], idx[1],
                     inter_mem_dub, *cur_quad_rule);

      if(has_delayed_entry){
        ghqCpp::simple_mem_stack<double> &my_stack = wmem::mem_stack();
        res += d_dat.grad(par_vec_dub, par_vec_gr, my_stack, delayed_entry_idx,
                          *cur_quad_rule, *cur_gh_quad_rule);
      }
    }

    // compute the gradient and return
    res.propagateToStart();

    for(vajoint_uint i = 0; i < par_idx.n_params_w_va<false>(); ++i)
      par_vec_gr[i] += par_vec_num[i].adjoint();

    // the covariance matrix parameters
    std::copy(par_vec_gr, par_vec_gr + par_idx.vcov_start<false>(), gr);
    std::copy(par_vec_gr + par_idx.va_mean<false>(),
              par_vec_gr + par_idx.va_mean_end<false>(),
              gr + par_idx.va_mean<true>());

    // the covariance matrix parameters
    std::fill(gr + par_idx.vcov_start<true>(),
              gr + par_idx.vcov_end<true>(), 0);
    std::fill(gr + par_idx.va_vcov<true>(),
              gr + par_idx.va_vcov_end<true>(), 0);

    log_chol::dpd_mat::get
      (p + par_idx.vcov_marker<true>(), par_idx.marker_info().size(),
       gr + par_idx.vcov_marker<true>(),
       par_vec_gr + par_idx.vcov_marker<false>(), inter_mem_dub);

    log_chol::dpd_mat::get
      (p + par_idx.vcov_surv<true>(), par_idx.n_shared_surv(),
       gr + par_idx.vcov_surv<true>(),
       par_vec_gr + par_idx.vcov_surv<false>(), inter_mem_dub);

    log_chol::dpd_mat::get
      (p + par_idx.vcov_vary<true>(), par_idx.n_shared(),
       gr + par_idx.vcov_vary<true>(),
       par_vec_gr + par_idx.vcov_vary<false>(), inter_mem_dub);

    log_chol::dpd_mat::get
      (p + par_idx.va_vcov<true>(), n_rng,
       gr + par_idx.va_vcov<true>(),
       par_vec_gr + par_idx.va_vcov<false>(), inter_mem_dub);

    return res.value();
  }

  double func(double const *point, lower_bound_caller const &caller) const {
    return nan_if_fail([&]{ return comp(point, nullptr, caller, false); });
  }

  double grad
  (double const * point, double * gr, lower_bound_caller const &caller) const {
    return nan_if_fail([&]{ return comp(point, gr, caller, true); });
  }

  bool thread_safe() const {
    return true;
  }
};

void lower_bound_caller::setup(double const *val, bool const comp_grad){
  setup_failed = false;

  try {
    // setup the global parameter vector
    par_vec.resize(par_idx->n_params<false>());
    vajoint_uint const n_wmem
    {many_max<vajoint_uint>
      (log_chol::pd_mat::n_wmem(par_idx->marker_info().size()),
       log_chol::pd_mat::n_wmem(par_idx->n_shared_surv()),
       log_chol::pd_mat::n_wmem(par_idx->n_shared()),
       m_dat->n_wmem(),
       kl_dat->n_wmem())
    };

    double * const wmem{wmem::get_double_mem(n_wmem)};
    log_chol::pd_mat::get
      (val + par_idx->vcov_marker<true>(), par_idx->marker_info().size(),
       par_vec.data() + par_idx->vcov_marker<false>(), wmem);
    log_chol::pd_mat::get
      (val + par_idx->vcov_surv<true>(), par_idx->n_shared_surv(),
       par_vec.data() + par_idx->vcov_surv<false>(), wmem);
    log_chol::pd_mat::get
      (val + par_idx->vcov_vary<true>(), par_idx->n_shared(),
       par_vec.data() + par_idx->vcov_vary<false>(), wmem);

    std::copy(val, val + par_idx->vcov_start<true>(), par_vec.data());

    // setup the market data, kl term, and the survival data
    m_dat->setup(par_vec.data(), wmem);
    kl_dat->setup(par_vec.data(), wmem,
                  optimze_survival ? lb_terms::all : lb_terms::markers);
  } catch(...){
    setup_failed = true;
  }
}

double lower_bound_caller::eval_func
  (lower_bound_term const &obj, double const * val){
  return obj.func(val, *this);
}
double lower_bound_caller::eval_grad
  (lower_bound_term const &obj, double const * val, double *gr){
  return obj.grad(val, gr, *this);
}

lower_bound_caller::lower_bound_caller
  (std::vector<lower_bound_term const*> &terms):
  // TODO: we need the non-const reference because we need to call setup
  par_idx
  {terms.size() == 0
    ? nullptr : const_cast<subset_params*>(&terms[0]->par_idx) },
  m_dat
  {terms.size() == 0
    ? nullptr : const_cast<marker::marker_dat*>(&terms[0]->m_dat)},
  kl_dat
  {terms.size() == 0
    ? nullptr : const_cast<kl_term*>(&terms[0]->kl_dat)},
  par_vec(par_idx->n_params<false>()) { }

/// psqn class to perform the optimization
using lb_optim = PSQN::optimizer
  <lower_bound_term, PSQN::R_reporter, PSQN::R_interrupter,
   lower_bound_caller>;

/**
 * class that holds the psqn class along with other objects that are referenced
 * needed to compute the function, gradient, and to find the maximum lower
 * bound.
 */
class problem_data {
  subset_params par_idx;
  marker::marker_dat m_dat;
  survival::survival_dat s_dat;
  kl_term kl_dat;
  survival::delayed_dat d_dat;
  std::unique_ptr<lb_optim> optim_obj;

public:
  problem_data(List markers, List survival_terms,
               unsigned const max_threads, List delayed_terms) {
    // handle the markers
    std::vector<marker::setup_marker_dat_helper> input_dat;
    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;

    input_dat.reserve(markers.size());
    bases_fix.reserve(markers.size());
    bases_rng.reserve(markers.size());

    for(auto m : markers){
      List marker = m;
      bases_fix.emplace_back(basis_from_list(marker["time_fixef"]));
      bases_rng.emplace_back(basis_from_list(marker["time_rng"]));

      NumericMatrix X{Rcpp::as<NumericMatrix>(marker["X"])},
                    fixef_design_varying
                      {Rcpp::as<NumericMatrix>(marker["fixef_design_varying"])},
                    rng_design_varying
                      {Rcpp::as<NumericMatrix>(marker["rng_design_varying"])};
      Rcpp::IntegerVector id = Rcpp::as<Rcpp::IntegerVector>(marker["id"]);
      NumericVector y = Rcpp::as<NumericVector>(marker["y"]),
                       time = Rcpp::as<NumericVector>(marker["time"]);

      vajoint_uint const n_fixef = X.nrow(),
                         n_obs   = X.ncol();

      input_dat.emplace_back(
        &X[0], n_fixef, n_obs, &id[0], &time[0], &y[0],
        &fixef_design_varying[0], fixef_design_varying.nrow(),
        &rng_design_varying[0], rng_design_varying.nrow());
      par_idx.add_marker
        ({n_fixef, bases_fix.back()->n_basis(), bases_rng.back()->n_basis()});
    }

    // handle the survival terms
    joint_bases::bases_vector bases_fix_surv;
    std::vector<survival::obs_input> surv_input;
    std::vector<simple_mat<double> > s_fixef_design,
                                     s_fixef_design_varying,
                                     s_rng_design_varying;
    std::vector<Rcpp::IntegerVector> s_id_vecs;

    surv_input.reserve(survival_terms.size());
    bases_fix_surv.reserve(survival_terms.size());
    s_fixef_design.reserve(survival_terms.size());
    s_fixef_design_varying.reserve(survival_terms.size());
    s_rng_design_varying.reserve(survival_terms.size());
    s_id_vecs.reserve(survival_terms.size());

    std::vector<std::vector<std::vector<int> > > ders;
    std::vector<unsigned> n_associations(markers.size());
    for(auto s : survival_terms){
      List surv = s;
      bases_fix_surv.emplace_back(basis_from_list(surv["time_fixef"]));

      NumericMatrix Z{Rcpp::as<NumericMatrix>(surv["Z"])},
                    fixef_design_varying
                      {Rcpp::as<NumericMatrix>(surv["fixef_design_varying"])},
                    rng_design_varying
                      {Rcpp::as<NumericMatrix>(surv["rng_design_varying"])};
      s_id_vecs.emplace_back(Rcpp::as<Rcpp::IntegerVector>(surv["id"]));
      NumericMatrix y{Rcpp::as<NumericMatrix>(surv["y"])};

      bool with_frailty{Rcpp::as<bool>(surv["with_frailty"])};

      // handle the integral/derivative argument to the basis expansions of the
      // markers
      ders.emplace_back();
      List s_ders = surv["ders"];
      for(SEXP ele : s_ders){
        Rcpp::IntegerVector ele_vec{ele};
        ders.back().emplace_back();
        ders.back().back().reserve(ele_vec.size());
        for(int val : ele_vec)
          ders.back().back().emplace_back(val);
      }
      if(ders.back().size() != static_cast<size_t>(markers.size()))
        throw std::runtime_error("Number of association parameters do not match the number of markers");
      for(unsigned i = 0; i < ders.back().size(); ++i)
        n_associations[i] = ders.back()[i].size();

      vajoint_uint const n_fixef = Z.nrow(),
                         n_obs   = Z.ncol();

      surv_input.emplace_back
        (survival::obs_input{n_obs, &y[0], &y[y.nrow()], &y[2 * y.nrow()]});
      s_fixef_design.emplace_back(&Z[0], n_fixef, n_obs);
      s_fixef_design_varying.emplace_back
        (&fixef_design_varying[0], fixef_design_varying.nrow(),
         fixef_design_varying.ncol());
      s_rng_design_varying.emplace_back
        (&rng_design_varying[0], rng_design_varying.nrow(),
         rng_design_varying.ncol());
      par_idx.add_surv
        ({n_fixef, bases_fix_surv.back()->n_basis(), n_associations,
         with_frailty});
    }

    // handle the terms for the delayed entry
    if(delayed_terms.size() != survival_terms.size())
      throw std::invalid_argument
        ("delayed_terms.size() != survival_terms.size()");

    std::vector<simple_mat<double> > d_fixef_design,
                        d_fixef_design_varying_mats,
                          d_rng_design_varying_mats;
    std::vector<Rcpp::IntegerVector> d_id_vecs;
    std::vector<Rcpp::NumericVector> delay_times;

    d_fixef_design.reserve(delayed_terms.size());
    d_fixef_design_varying_mats.reserve(delayed_terms.size());
    d_rng_design_varying_mats.reserve(delayed_terms.size());
    d_id_vecs.reserve(delayed_terms.size());
    delay_times.reserve(delayed_terms.size());

    for(auto d : delayed_terms){
      List delay = d;

      NumericMatrix Z{Rcpp::as<NumericMatrix>(delay["Z"])},
                    fixef_design_varying
                      {Rcpp::as<NumericMatrix>(delay["fixef_design_varying"])},
                    rng_design_varying
                      {Rcpp::as<NumericMatrix>(delay["rng_design_varying"])};
      d_id_vecs.emplace_back(Rcpp::as<Rcpp::IntegerVector>(delay["id"]));
      delay_times.emplace_back(Rcpp::as<NumericVector>(delay["y"]));

      vajoint_uint const n_fixef = Z.nrow(),
                         n_obs   = Z.ncol();
      d_fixef_design.emplace_back(&Z[0], n_fixef, n_obs);
      d_fixef_design_varying_mats.emplace_back
        (&fixef_design_varying[0], fixef_design_varying.nrow(),
         fixef_design_varying.ncol());
      d_rng_design_varying_mats.emplace_back
        (&rng_design_varying[0], rng_design_varying.nrow(),
         rng_design_varying.ncol());
    }

    // construct the objects to compute the different terms of the lower bound
    auto dat_n_idx =
      marker::get_comp_dat(input_dat, par_idx, bases_fix, bases_rng);
    kl_dat = kl_term(par_idx);
    m_dat = std::move(dat_n_idx.dat);
    s_dat = survival::survival_dat
      (bases_fix_surv, bases_rng, s_fixef_design, s_fixef_design_varying,
       s_rng_design_varying, par_idx, surv_input, ders);

    // check that all ids are sorted
    for(auto b = dat_n_idx.id.begin(); b != dat_n_idx.id.end(); ++b)
      if(b != dat_n_idx.id.begin() && *b < *(b - 1))
        throw std::invalid_argument("Marker ids are not sorted");

    for(size_t i = 0; i < s_id_vecs.size(); ++i)
      for(auto b = s_id_vecs[i].begin(); b != s_id_vecs[i].end(); ++b)
        if(b != s_id_vecs[i].begin() && *b < *(b - 1))
          throw std::invalid_argument
          ("ids for survival type " + std::to_string(i + 1) + " are not sorted");

    for(size_t i = 0; i < d_id_vecs.size(); ++i)
      for(auto b = d_id_vecs[i].begin(); b != d_id_vecs[i].end(); ++b)
        if(b != d_id_vecs[i].begin() && *b < *(b - 1))
          throw std::invalid_argument
          ("ids for delayed entry type " + std::to_string(i + 1) + " are not sorted");

    // have to make a pass through the data to construct the delayed entry
    // object
    std::vector<std::vector<survival::delayed_dat::cluster_obs> >
      delayed_cluster_obs;
    std::vector<int> delayed_cluster_ids;
    for(auto ids : d_id_vecs){
      delayed_cluster_obs.reserve(ids.size());
      delayed_cluster_ids.reserve(ids.size());
    }

    {
      std::vector<vajoint_uint> indices(d_id_vecs.size(), 0);
      auto all_at_end = [&]{
        bool out{true};
        for(size_t i = 0; out && i < indices.size(); ++i)
          out &= indices[i] >= d_id_vecs[i].size();
        return out;
      };

      while(!all_at_end()){
        // find the next smallest id
        int smallest_id{std::numeric_limits<int>::max()};
        for(size_t i = 0; i < indices.size(); ++i)
          if(indices[i] < d_id_vecs[i].size())
            smallest_id = std::min(smallest_id, d_id_vecs[i][indices[i]]);

        // create a cluster with ids that match the smallest id
        delayed_cluster_ids.emplace_back(smallest_id);
        delayed_cluster_obs.emplace_back();
        auto &next_cluster = delayed_cluster_obs.back();

        for(vajoint_uint type = 0; type < indices.size(); ++type){
          while(indices[type] < d_id_vecs[type].size() &&
                d_id_vecs[type][indices[type]] == smallest_id){
            next_cluster.emplace_back(
              survival::delayed_dat::cluster_obs
                {type, indices[type], delay_times[type][indices[type]]});
            ++indices[type];
          }
        }

        next_cluster.shrink_to_fit();
      }
    }
    delayed_cluster_obs.shrink_to_fit();

    d_dat = survival::delayed_dat
      (bases_fix_surv, bases_rng, d_fixef_design, d_fixef_design_varying_mats,
       d_rng_design_varying_mats, par_idx, delayed_cluster_obs, ders);

    // create the object to use for optimization
    std::vector<lower_bound_term> ele_funcs;
    ele_funcs.reserve(dat_n_idx.id.size());

    {
      // add the element functions
      auto id_marker = dat_n_idx.id.begin();
      std::vector<Rcpp::IntegerVector::iterator> s_indices(s_id_vecs.size());
      for(size_t i = 0; i < s_indices.size(); ++i)
        s_indices[i] = s_id_vecs[i].begin();
      size_t delayed_idx{};

      // returns true if all iterators are at the end
      auto all_at_end = [&]{
        if(id_marker != dat_n_idx.id.end())
          return false;
        for(size_t i = 0; i < s_indices.size(); ++i)
          if(s_indices[i] != s_id_vecs[i].end())
            return false;
        return true;
      };

      while(!all_at_end()){
        // find the lowest id as the ids are supposed to be sorted
        int cur_id{std::numeric_limits<int>::max()};
        if(id_marker != dat_n_idx.id.end())
          cur_id = std::min(*id_marker, cur_id);
        for(size_t i = 0; i < s_indices.size(); ++i)
          if(s_indices[i] != s_id_vecs[i].end())
            cur_id = std::min(*s_indices[i], cur_id);

        // add the observation where the id does match
        ele_funcs.emplace_back(par_idx, m_dat, s_dat, kl_dat, d_dat);
        auto &ele_func = ele_funcs.back();
        while(id_marker != dat_n_idx.id.end() && *id_marker == cur_id)
          ele_func.add_marker_index
            (std::distance(dat_n_idx.id.begin(), id_marker++));

        for(size_t i = 0; i < s_indices.size(); ++i)
          while(s_indices[i] != s_id_vecs[i].end() && *s_indices[i] == cur_id)
            ele_func.add_surv_index
              (std::distance(s_id_vecs[i].begin(), s_indices[i]++), i);

        if(delayed_idx < delayed_cluster_ids.size() &&
            delayed_cluster_ids[delayed_idx] < cur_id)
          throw std::invalid_argument(
              "delayed entry wihtout an outcome (id "
                + std::to_string(delayed_cluster_ids[delayed_idx])
                + " " + std::to_string(cur_id) + ")");

        if(delayed_idx < delayed_cluster_ids.size() &&
            delayed_cluster_ids[delayed_idx] == cur_id){
          ele_func.add_delayed_entry(delayed_idx);
          ++delayed_idx;
        }

        ele_func.shrink_to_fit();
      }
    }

    ele_funcs.shrink_to_fit();
    optim_obj.reset(new lb_optim(ele_funcs, max_threads));
  }

  lb_optim & optim(){
    return *optim_obj;
  }

  /// sets the cached expansions for the survival terms
  void set_cached_expansions(survival::node_weight const &nws){
    s_dat.set_cached_expansions(nws);
    d_dat.set_cached_expansions(nws, wmem::mem_stack());
  }

  /***
   * clears the cached expansions for the survival terms. Thus, they are
   * recomputed every time
   */
  void clear_cached_expansions(){
    s_dat.clear_cached_expansions();
    d_dat.clear_cached_expansions();
  }

  lb_optim const & optim() const {
    return *optim_obj;
  }

  subset_params const &params() const {
    return par_idx;
  }

  vajoint_uint n_markers_obs() const {
    return m_dat.n_obs();
  }

  vajoint_uint n_surv_obs(vajoint_uint const type) const {
    return s_dat.n_terms(type);
  }

  vajoint_uint n_surv_types() const {
    return s_dat.n_outcomes();
  }

  void set_n_threads(unsigned const n_threads){
    optim().set_n_threads(n_threads);
    wmem::setup_working_memory(n_threads);
    set_num_tapes(n_threads);
  }
};

inline void check_par_length(problem_data const &obj, NumericVector par){
  if(obj.optim().n_par != static_cast<size_t>(par.size()))
    throw std::invalid_argument("invalid parameter size");
}

inline void set_or_clear_cached_expansions
  (problem_data &dat, survival::node_weight const &nws, bool const do_set){
  if(do_set)
    dat.set_cached_expansions(nws);
  else
    dat.clear_cached_expansions();
}

/// returns a pointer to problem_data object
// [[Rcpp::export(".joint_ms_ptr", rng = false)]]
SEXP joint_ms_ptr
  (List markers, List survival_terms, unsigned const max_threads,
   List delayed_terms){
  profiler pp(".joint_ms_ptr");

  return Rcpp::XPtr<problem_data>
    (new problem_data(markers, survival_terms, max_threads, delayed_terms));
}

/// returns the number of lower bound terms of different types
// [[Rcpp::export(rng = false)]]
List joint_ms_n_terms(SEXP ptr){
  profiler pp("joint_ms_n_terms");

  Rcpp::XPtr<problem_data> obj(ptr);

  Rcpp::IntegerVector surv_count(obj->n_surv_types());
  for(vajoint_uint i = 0; i < obj->n_surv_types(); ++i)
    surv_count[i] = obj->n_surv_obs(i);

  return List::create(
    Rcpp::_("Marker terms") = obj->n_markers_obs(),
    Rcpp::_("Survival terms") = std::move(surv_count),
    Rcpp::_("Number of clusters") = obj->optim().get_ele_funcs().size());
}

/// evaluates the lower bound
// [[Rcpp::export(rng = false)]]
double joint_ms_eval_lb
  (NumericVector val, SEXP ptr, unsigned const n_threads, List quad_rule,
   bool const cache_expansions, List gh_quad_rule){
  profiler pp("joint_ms_eval_lb");

  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  survival::node_weight quad_rule_use{node_weight_from_list(quad_rule)};
  cur_quad_rule = &quad_rule_use;
  ghqCpp::ghq_data gh_quad_rule_use{gh_node_weight_from_list(gh_quad_rule)};
  cur_gh_quad_rule = &gh_quad_rule_use;

  set_or_clear_cached_expansions(*obj, quad_rule_use, cache_expansions);
  obj->set_n_threads(n_threads);
  double const out{obj->optim().eval(&val[0], nullptr, false)};
  wmem::rewind_all();

  return out;
}

/// evaluates the gradient of the lower bound
// [[Rcpp::export(rng = false)]]
NumericVector joint_ms_eval_lb_gr
  (NumericVector val, SEXP ptr, unsigned const n_threads, List quad_rule,
   bool const cache_expansions, List gh_quad_rule){
  profiler pp("joint_ms_eval_lb_gr");

  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  survival::node_weight quad_rule_use{node_weight_from_list(quad_rule)};
  cur_quad_rule = &quad_rule_use;
  ghqCpp::ghq_data gh_quad_rule_use{gh_node_weight_from_list(gh_quad_rule)};
  cur_gh_quad_rule = &gh_quad_rule_use;

  set_or_clear_cached_expansions(*obj, quad_rule_use, cache_expansions);

  NumericVector grad(val.size());
  obj->set_n_threads(n_threads);
  grad.attr("value") = obj->optim().eval(&val[0], &grad[0], true);
  wmem::rewind_all();

  return grad;
}

/// computes the Hessian with numerical differentiation
// [[Rcpp::export(".joint_ms_hess", rng = false)]]
Eigen::SparseMatrix<double> joint_ms_hess
  (NumericVector val, SEXP ptr,  List quad_rule, bool const cache_expansions,
   double const eps, double const scale, double const tol,
   unsigned const order, List gh_quad_rule){
  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  survival::node_weight quad_rule_use{node_weight_from_list(quad_rule)};
  cur_quad_rule = &quad_rule_use;
  ghqCpp::ghq_data gh_quad_rule_use{gh_node_weight_from_list(gh_quad_rule)};
  cur_gh_quad_rule = &gh_quad_rule_use;

  set_or_clear_cached_expansions(*obj, quad_rule_use, cache_expansions);

  return obj->optim().true_hess_sparse(&val[0], eps, scale, tol, order);
}

/// returns the names of the parameters
// [[Rcpp::export(rng = false)]]
List joint_ms_parameter_names(SEXP ptr){
  Rcpp::XPtr<problem_data> obj(ptr);

  auto param_names = obj->params().param_names(true);
  auto va_param_names = obj->params().va_param_names(true);

  Rcpp::CharacterVector param_names_out(param_names.size());
  for(size_t i = 0; i < param_names.size(); ++i)
    param_names_out[i] = param_names[i];

  Rcpp::CharacterVector va_param_names_out(va_param_names.size());
  for(size_t i = 0; i < va_param_names.size(); ++i)
    va_param_names_out[i] = va_param_names[i];

  return List::create(
    Rcpp::_("param_names") = std::move(param_names_out),
    Rcpp::_("VA_param_names") = std::move(va_param_names_out));
}

// TODO: need to test this function
/// returns indices of the parameters
// [[Rcpp::export(rng = false)]]
List joint_ms_parameter_indices(SEXP ptr){
  Rcpp::XPtr<problem_data> obj(ptr);

  auto &params = obj->params();
  List m_out(params.marker_info().size()),
       s_out(params.surv_info().size());

  // handle the fixed effects
  for(size_t i = 0; i < params.marker_info().size(); ++i){
    auto &info = params.marker_info()[i];
    Rcpp::IntegerVector fixef(info.n_fix),
                   fixef_vary(info.n_variying);

    std::iota(fixef.begin(), fixef.end(), info.idx_fix + 1);
    std::iota(fixef_vary.begin(), fixef_vary.end(), info.idx_varying + 1);

    m_out[i] = List::create(
      Rcpp::_("fixef") = std::move(fixef),
      Rcpp::_("fixef_vary") = std::move(fixef_vary));
  }

  for(size_t i = 0; i < params.surv_info().size(); ++i){
    auto &info = params.surv_info()[i];
    unsigned const n_assoc{std::accumulate(info.n_associations.begin(),
                                           info.n_associations.end(), 0U)};

    Rcpp::IntegerVector fixef(info.n_fix),
                   fixef_vary(info.n_variying),
                 associations(n_assoc);

    std::iota(fixef.begin(), fixef.end(), info.idx_fix + 1);
    std::iota(fixef_vary.begin(), fixef_vary.end(), info.idx_varying + 1);
    std::iota
      (associations.begin(), associations.end(), info.idx_association + 1);

    s_out[i] = List::create(
      Rcpp::_("fixef") = std::move(fixef),
      Rcpp::_("fixef_vary") = std::move(fixef_vary),
      Rcpp::_("associations") = std::move(associations));
  }

  // handle the covariance matrices
  vajoint_uint const n_vcov_marker
    {static_cast<vajoint_uint>(dim_tri(params.marker_info().size()))};
  Rcpp::IntegerVector vcov_marker(n_vcov_marker);
  std::iota
    (vcov_marker.begin(), vcov_marker.end(), params.vcov_marker<true>() + 1);

  vajoint_uint const n_vcov_surv
    {static_cast<vajoint_uint>(dim_tri(params.n_shared_surv()))};
  Rcpp::IntegerVector vcov_surv(n_vcov_surv);
  std::iota
    (vcov_surv.begin(), vcov_surv.end(), params.vcov_surv<true>() + 1);

  vajoint_uint const n_vcov_vary
    {static_cast<vajoint_uint>(dim_tri(params.n_shared()))};
  Rcpp::IntegerVector vcov_vary(n_vcov_vary);
  std::iota
    (vcov_vary.begin(), vcov_vary.end(), params.vcov_vary<true>() + 1);

  List vcovs = List::create(
    Rcpp::_("vcov_marker") = std::move(vcov_marker),
    Rcpp::_("vcov_surv") = std::move(vcov_surv),
    Rcpp::_("vcov_vary") = std::move(vcov_vary));

  return List::create(
    Rcpp::_("markers") = std::move(m_out),
    Rcpp::_("survival") = std::move(s_out),
    Rcpp::_("vcovs") = std::move(vcovs),
    Rcpp::_("va_params_start") = params.va_mean<true>() + 1,
    Rcpp::_("n_va_params") = params.n_va_params<true>(),
    Rcpp::_("va_dim") = params.va_mean_end() - params.va_mean());
}

/// returns the number of parameters
// [[Rcpp::export(rng = false)]]
int joint_ms_n_params(SEXP ptr){
  Rcpp::XPtr<problem_data> obj(ptr);
  return static_cast<int>(obj->optim().n_par);
}

/// optimizes the private parameters
// [[Rcpp::export(rng = false)]]
NumericVector opt_priv
  (NumericVector val, SEXP ptr,
   double const rel_eps, unsigned const max_it, unsigned const n_threads,
   double const c1, double const c2, List quad_rule,
   bool const cache_expansions, double const gr_tol, List gh_quad_rule){
  profiler pp("opt_priv");

  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  survival::node_weight quad_rule_use{node_weight_from_list(quad_rule)};
  cur_quad_rule = &quad_rule_use;
  ghqCpp::ghq_data gh_quad_rule_use{gh_node_weight_from_list(gh_quad_rule)};
  cur_gh_quad_rule = &gh_quad_rule_use;

  set_or_clear_cached_expansions(*obj, quad_rule_use, cache_expansions);

  NumericVector par = clone(val);
  obj->set_n_threads(n_threads);
  double const res = obj->optim().
    optim_priv(&par[0], rel_eps, max_it, c1, c2, gr_tol);
  par.attr("value") = res;
  wmem::rewind_all();

  return par;
}

/// optimizes the lower bound using the psqn package
// [[Rcpp::export(rng = false)]]
List joint_ms_opt_lb
  (NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1,
   double const c2, bool const use_bfgs, unsigned const trace,
   double const cg_tol, bool const strong_wolfe, size_t const max_cg,
   unsigned const pre_method, List quad_rule, Rcpp::IntegerVector mask,
   bool const cache_expansions, bool const only_markers,
   double const gr_tol, List gh_quad_rule){
  profiler pp("joint_ms_opt_lb");

  lower_bound_caller::optimze_survival = !only_markers;
  struct reset_opt_surv {
    ~reset_opt_surv() { lower_bound_caller::optimze_survival = true; }
  } reset_opt_surv_obj;

  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  obj->optim().set_masked(mask.begin(), mask.end());
  struct clear_masked {
    problem_data &dat;
    clear_masked(problem_data &dat): dat{dat} { }
    ~clear_masked() { dat.optim().clear_masked(); }
  } clear_m(*obj);

  survival::node_weight quad_rule_use{node_weight_from_list(quad_rule)};
  cur_quad_rule = &quad_rule_use;
  ghqCpp::ghq_data gh_quad_rule_use{gh_node_weight_from_list(gh_quad_rule)};
  cur_gh_quad_rule = &gh_quad_rule_use;

  set_or_clear_cached_expansions(*obj, quad_rule_use, cache_expansions);

  NumericVector par = clone(val);
  obj->set_n_threads(n_threads);
  auto res = obj->optim().optim(&par[0], rel_eps, max_it, c1, c2,
                                use_bfgs, trace, cg_tol, strong_wolfe, max_cg,
                                static_cast<PSQN::precondition>(pre_method),
                                gr_tol);

  NumericVector counts = NumericVector::create(
    res.n_eval, res.n_grad,  res.n_cg);
  counts.names() =
    Rcpp::CharacterVector::create("function", "gradient", "n_cg");
  wmem::rewind_all();

  int const info{static_cast<int>(res.info)};
  return List::create(
    Rcpp::_["par"] = par, Rcpp::_["value"] = res.value,
    Rcpp::_["info"] = info, Rcpp::_["counts"] = counts,
    Rcpp::_["convergence"] =  res.info == PSQN::info_code::converged);
}


/**
 * class to evaluate the negative log likelihood and its gradient for a
 * proportional hazard model
 */
class ph_model {
  std::unique_ptr<joint_bases::basisMixin> expansion;
  simple_mat<double> Z,
                     fixef_design_varying,
                     rng_design_varying,
                     surv;

  survival::expected_cum_hazzard cum_haz;

  std::array<size_t, 2> const n_wmem_v;

  double lb(vajoint_uint const idx) const{
    return *surv.col(idx);
  }
  double ub(vajoint_uint const idx) const{
    return surv.col(idx)[1];
  }
  double event(vajoint_uint const idx) const{
    return surv.col(idx)[2];
  }

  template<class T>
  T eval(T const * param, survival::node_weight const &quad_rule,
           vajoint_uint const start, vajoint_uint const end,
           T * T_mem, double * wk_mem, double const va_var_in) const{

    T const * const fixef{param},
            * const fixef_vary{param + Z.n_rows()};

    T out{0.};
    T association{0};
    T va_mean(0);
    T va_var(va_var_in);
    for(vajoint_uint i = start; i < end; ++i){
      // the hazard term
      if(event(i) > 0){
        out -= cfaad::dotProd(Z.col(i), Z.col(i) + Z.n_rows(), fixef);
        (*expansion)(wk_mem, wk_mem + expansion->n_basis(), ub(i),
                     fixef_design_varying.col(i));
        out -= cfaad::dotProd
          (wk_mem, wk_mem + expansion->n_basis(), fixef_vary);
      }

      // the cumulative hazard term
      out += cum_haz
        (quad_rule, lb(i), ub(i), Z.col(i), fixef_design_varying.col(i),
         rng_design_varying.col(i), fixef, fixef_vary, &association, &va_mean,
         &va_var, T_mem, wk_mem, nullptr);
    }

    return out;
  }

public:
  ph_model
    (joint_bases::basisMixin const * expansion_in, simple_mat<double> const &Z,
     simple_mat<double> const &fixef_design_varying,
     simple_mat<double> const &rng_design_varying,
     simple_mat<double> const &surv, bool const with_frailty):
  expansion{expansion_in->clone()}, Z{Z},
  fixef_design_varying{fixef_design_varying},
  rng_design_varying{rng_design_varying},
  surv{surv},
  cum_haz{*expansion, survival::bases_vector{}, Z.n_rows(),
          std::vector<std::vector<int> >{}, with_frailty},
  n_wmem_v{cum_haz.n_wmem()[0],
           std::max(cum_haz.n_wmem()[1],
                    expansion->n_wmem() + expansion->n_basis())}
  { }

  vajoint_uint n_params() const {
    return Z.n_rows() + expansion->n_basis();
  }

  std::array<size_t, 2> const & n_wmem() const {
    return n_wmem_v;
  }

  double eval(double const *param, survival::node_weight const &quad_rule,
              double const va_var) const {
    double * const w1{wmem::get_double_mem(n_wmem()[0])},
           * const w2{wmem::get_double_mem(n_wmem()[1])};
    double const out(eval(param, quad_rule, 0, Z.n_cols(), w1, w2, va_var));
    wmem::rewind_all();
    return out;
  }

  double gr
  (double const *param, double *gr,
   survival::node_weight const &quad_rule, double const va_var) const {
    Number::tape->clear();
    Number * const param_num{wmem::get_Number_mem(n_params())};
    cfaad::convertCollection(param, param + n_params(), param_num);

    Number * const w1{wmem::get_Number_mem(n_wmem()[0])};
    double * const w2{wmem::get_double_mem(n_wmem()[1])};

    double out{};
    for(vajoint_uint i = 0; i < Z.n_cols(); ){
      Number::tape->rewind();
      cfaad::putOnTape(param_num, param_num + n_params());

      vajoint_uint const inc{std::min<vajoint_uint>(Z.n_cols() - i, 256)};
      Number res = eval(param_num, quad_rule, i, i + inc, w1, w2, va_var);
      i += inc;
      res.propagateToStart();
      out += res.value();

      for(vajoint_uint i = 0; i < n_params(); ++i)
        gr[i] += param_num[i].adjoint();
    }

    Number::tape->clear();
    wmem::rewind_all();
    return out;
  }
};

// [[Rcpp::export(rng = false)]]
List ph_ll
  (List time_fixef, NumericMatrix Z, NumericMatrix surv,
   bool const with_frailty, NumericMatrix fixef_design_varying,
   NumericMatrix rng_design_varying){
  profiler pp("ph_ll");

  auto expansion = basis_from_list(time_fixef);
  simple_mat Z_sm(&Z[0], Z.nrow(), Z.ncol()),
             fixef_design_varying_sm(&fixef_design_varying[0],
                                     fixef_design_varying.nrow(),
                                     fixef_design_varying.ncol()),
             rng_design_varying_sm(&rng_design_varying[0],
                                   rng_design_varying.nrow(),
                                   rng_design_varying.ncol()),
             surv_sm(&surv[0], surv.nrow(), surv.ncol());

  // basic check
  if(surv_sm.n_rows() != 3)
    throw std::invalid_argument("surv.nrow() != 3");
  if(Z_sm.n_cols() != surv_sm.n_cols())
    throw std::invalid_argument("Z_sm.n_cols() != surv_sm.n_cols()");
  if(fixef_design_varying_sm.n_cols() != surv_sm.n_cols())
    throw std::invalid_argument("fixef_design_varying_sm.n_cols() != surv_sm.n_cols()");
  if(rng_design_varying_sm.n_cols() != surv_sm.n_cols())
    throw std::invalid_argument("rng_design_varying_sm.n_cols() != surv_sm.n_cols()");

  Rcpp::XPtr<ph_model> ptr
    (new ph_model
       (expansion.get(), Z_sm, fixef_design_varying_sm, rng_design_varying_sm,
        surv_sm, with_frailty));
  vajoint_uint const n_params = ptr->n_params();

  return List::create(Rcpp::_("n_params") = n_params,
                      Rcpp::_("ptr") = std::move(ptr));
}

// [[Rcpp::export(rng = false)]]
double ph_eval
  (SEXP ptr, NumericVector par, List quad_rule, double const va_var){
  profiler pp("ph_eval");

  Rcpp::XPtr<ph_model> comp_obj(ptr);
  if(par.size() != static_cast<R_len_t>(comp_obj->n_params()))
    throw std::invalid_argument("par.size() != n_params");

  return comp_obj->eval(&par[0], node_weight_from_list(quad_rule), va_var);
}

// [[Rcpp::export(rng = false)]]
NumericVector ph_grad
  (SEXP ptr, NumericVector par, List quad_rule, double const va_var){
  profiler pp("ph_grad");

  Rcpp::XPtr<ph_model> comp_obj(ptr);
  if(par.size() != static_cast<R_len_t>(comp_obj->n_params()))
    throw std::invalid_argument("par.size() != n_params");

  NumericVector out(comp_obj->n_params(), 0);
  out.attr("logLik") =
    comp_obj->gr(&par[0], &out[0], node_weight_from_list(quad_rule), va_var);
  return out;
}
