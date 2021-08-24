#include "arma-wrap.h"
#include "marker-term.h"
#include "psqn.h"
#include "psqn-reporter.h"
#include "log-cholesky.h"
#include "kl-term.h"
#include <AAD.h>
#include <numeric>
#include <limits>

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

} // namespace

using cfaad::Number;

/// returns a pointer to an expansions given an R List with data.
std::unique_ptr<joint_bases::basisMixin> get_basis_from_list(Rcpp::List dat){
  if(Rf_inherits(dat, "poly_term")){
    Rcpp::List coefs = dat["coefs"];
    arma::vec alpha{Rcpp::as<arma::vec>(coefs["alpha"])},
              norm2{Rcpp::as<arma::vec>(coefs["norm2"])};
    bool const raw{Rcpp::as<bool>(dat["raw"])},
         intercept{Rcpp::as<bool>(dat["intercept"])};

    return raw
      ? std::make_unique<joint_bases::orth_poly>(alpha.size(), intercept)
      : std::make_unique<joint_bases::orth_poly>(alpha, norm2, intercept);

  } else if(Rf_inherits(dat, "bs_term")){
    arma::vec i_knots{Rcpp::as<arma::vec>(dat["knots"])},
              b_knots{Rcpp::as<arma::vec>(dat["Boundary.knots"])};

    bool const intercept{Rcpp::as<bool>(dat["intercept"])};
    vajoint_uint const degree{Rcpp::as<vajoint_uint>(dat["degree"])};

    return std::make_unique<joint_bases::bs>
      (b_knots, i_knots, intercept, degree + 1);

  } else if(Rf_inherits(dat, "ns_term")){
    arma::vec i_knots{Rcpp::as<arma::vec>(dat["knots"])},
              b_knots{Rcpp::as<arma::vec>(dat["Boundary.knots"])};

    bool const intercept{Rcpp::as<bool>(dat["intercept"])};
    vajoint_uint const degree{Rcpp::as<vajoint_uint>(dat["degree"])};

    return std::make_unique<joint_bases::ns>
      (b_knots, i_knots, intercept, degree + 1);

  }

  throw std::invalid_argument("expansions is not implemented");
  return std::make_unique<joint_bases::orth_poly>(2, false);
}

/// functions which is exported only to check that the bases work as in R
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix eval_expansion(Rcpp::List dat, Rcpp::NumericVector const x){
  auto basis = get_basis_from_list(dat);
  Rcpp::NumericMatrix out(basis->n_basis(), x.size());

  wmem::clear_all();
  double *wmem = wmem::get_double_mem(basis->n_wmem());
  for(R_len_t i = 0; i < x.size(); ++i)
    (*basis)(&out.column(i)[0], wmem, x[i]);

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
  lower_bound_caller(std::vector<lower_bound_term const*> &);

  void setup(double const *val, bool const comp_grad);
  double eval_func(lower_bound_term const &obj, double const * val);
  double eval_grad(lower_bound_term const &obj, double const * val,
                   double *gr);
};

class lower_bound_term {
  subset_params const &par_idx;
  marker::marker_dat const &m_dat;
  kl_term const &kl_dat;

  std::vector<vajoint_uint> marker_indices;
  size_t n_global, n_private;

  friend lower_bound_caller;
public:
  lower_bound_term
  (subset_params const &par_idx, marker::marker_dat const &m_dat,
   kl_term const &kl_dat):
  par_idx(par_idx), m_dat(m_dat), kl_dat(kl_dat),
  n_global(par_idx.n_params<true>()),
  n_private{par_idx.n_va_params<true>()}
  { }

  void add_marker_index(vajoint_uint idx){
    marker_indices.emplace_back(idx);
  }

  void shrink_to_fit() {
    marker_indices.shrink_to_fit();
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

    vajoint_uint const n_rng{par_idx.va_mean_end() - par_idx.va_mean()};
    wmem::rewind();

    if(!comp_grad){
      // find the required working memory
      vajoint_uint const n_wmem =
        many_max<vajoint_uint>(log_chol::pd_mat::n_wmem(n_rng),
                               m_dat.n_wmem(),
                               kl_dat.n_wmem());

      double * const inter_mem{wmem::get_double_mem(n_wmem)};
      double * const par_vec{wmem::get_double_mem(par_idx.n_params_w_va())};

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

      return res;
    }

    // find the required working memory
    vajoint_uint const n_wmem_num
      {many_max<vajoint_uint>(m_dat.n_wmem())};
    vajoint_uint const n_wmem_dub
      {many_max<vajoint_uint>
        (log_chol:: pd_mat::n_wmem(n_rng),
         log_chol::dpd_mat::n_wmem(n_rng),
         log_chol::dpd_mat::n_wmem(par_idx.marker_info().size()),
         log_chol::dpd_mat::n_wmem(par_idx.surv_info().size()),
         kl_dat.n_wmem())};

    Number * const inter_mem_num{wmem::get_Number_mem(n_wmem_num)};
    double * const inter_mem_dub{wmem::get_double_mem(n_wmem_dub)};

    // setup the parameter vectors
    double * const par_vec_dub
      {wmem::get_double_mem(par_idx.n_params_w_va<false>())},
           * const par_vec_gr
      {wmem::get_double_mem(par_idx.n_params_w_va<false>())};
    Number * const par_vec_num
      {wmem::get_Number_mem(par_idx.n_params_w_va<false>())};

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

    // compute the lower bound terms
    Number res{0};
    for(vajoint_uint idx : marker_indices)
      res += m_dat(par_vec_num, inter_mem_num, idx);

    res += kl_dat.grad(par_vec_gr, par_vec_dub, inter_mem_dub);

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
      (p + par_idx.vcov_surv<true>(), par_idx.surv_info().size(),
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
      (log_chol::pd_mat::n_wmem(par_idx->surv_info().size()),
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
    kl_dat->setup(par_vec.data(), wmem);
    // TODO: setup for the survival terms
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
  kl_term kl_dat;
  std::unique_ptr<lb_optim> optim_obj;

public:
  problem_data(Rcpp::List markers, unsigned const max_threads) {
    // handle the markers
    std::vector<marker::setup_marker_dat_helper> input_dat;
    joint_bases::bases_vector bases_fix;
    joint_bases::bases_vector bases_rng;

    for(auto m : markers){
      Rcpp::List marker = m;
      bases_fix.emplace_back(get_basis_from_list(marker["time_fixef"]));
      bases_rng.emplace_back(get_basis_from_list(marker["time_rng"]));

      Rcpp::NumericMatrix X{marker["X"]};
      Rcpp::IntegerVector id = Rcpp::as<Rcpp::IntegerVector>(marker["id"]);
      Rcpp::NumericVector y = Rcpp::as<Rcpp::NumericVector>(marker["y"]),
                       time = Rcpp::as<Rcpp::NumericVector>(marker["time"]);

      vajoint_uint const n_fixef = X.nrow(),
                         n_obs   = X.ncol();

      input_dat.emplace_back(&X[0], n_fixef, n_obs, &id[0], &time[0], &y[0]);
      par_idx.add_marker
        ({n_fixef, bases_fix.back()->n_basis(), bases_rng.back()->n_basis()});
    }

    // TODO: can first be constructed later when all objects have been added to
    //       par_idx
    auto dat_n_idx =
      marker::get_comp_dat(input_dat, par_idx, bases_fix, bases_rng);
    m_dat = std::move(dat_n_idx.dat);

    kl_dat = kl_term(par_idx);

    // TODO: create the other terms

    // create the object to use for optimization
    std::vector<lower_bound_term> ele_funcs;
    ele_funcs.reserve(dat_n_idx.id.size());

    {
      vajoint_uint m_idx{};
      auto id_marker = dat_n_idx.id.begin();
      while(id_marker != dat_n_idx.id.end()){
        int const cur_id{*id_marker++};
        ele_funcs.emplace_back(par_idx, m_dat, kl_dat);
        ele_funcs.back().add_marker_index(m_idx++);

        for(; id_marker != dat_n_idx.id.end() && *id_marker == cur_id;
            ++id_marker)
          ele_funcs.back().add_marker_index(m_idx++);

        ele_funcs.back().shrink_to_fit();
      }
    }

    ele_funcs.shrink_to_fit();
    optim_obj.reset(new lb_optim(ele_funcs, max_threads));
  }

  lb_optim & optim(){
    return *optim_obj;
  }

  lb_optim const & optim() const {
    return *optim_obj;
  }

  subset_params const &params() const {
    return par_idx;
  }

  void set_n_threads(unsigned const n_threads){
    optim().set_n_threads(n_threads);
    wmem::setup_working_memory(n_threads);
    // TODO: we have not set the Number::tape for all threads!
  }
};

inline void check_par_length(problem_data const &obj, Rcpp::NumericVector par){
  if(obj.optim().n_par != static_cast<size_t>(par.size()))
    throw std::invalid_argument("invalid parameter size");
}

/// returns a pointer to problem_data object
// [[Rcpp::export(".joint_ms_ptr", rng = false)]]
SEXP joint_ms_ptr(Rcpp::List markers, unsigned const max_threads){
  return Rcpp::XPtr<problem_data>(new problem_data(markers, max_threads));
}

/// evaluates the lower bound
// [[Rcpp::export(rng = false)]]
double joint_ms_eval_lb
  (Rcpp::NumericVector val, SEXP ptr, unsigned const n_threads){
  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  obj->set_n_threads(n_threads);
  return obj->optim().eval(&val[0], nullptr, false);
}

/// evaluates the gradient of the lower bound
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector joint_ms_eval_lb_gr
  (Rcpp::NumericVector val, SEXP ptr, unsigned const n_threads){
  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  Rcpp::NumericVector grad(val.size());
  obj->set_n_threads(n_threads);
  grad.attr("value") = obj->optim().eval(&val[0], &grad[0], true);

  return grad;
}

/// returns the names of the parameters
// [[Rcpp::export(rng = false)]]
Rcpp::List joint_ms_parameter_names(SEXP ptr){
  Rcpp::XPtr<problem_data> obj(ptr);

  auto param_names = obj->params().param_names(true);
  auto va_param_names = obj->params().va_param_names(true);

  Rcpp::CharacterVector param_names_out(param_names.size());
  for(size_t i = 0; i < param_names.size(); ++i)
    param_names_out[i] = param_names[i];

  Rcpp::CharacterVector va_param_names_out(va_param_names.size());
  for(size_t i = 0; i < va_param_names.size(); ++i)
    va_param_names_out[i] = va_param_names[i];

  return Rcpp::List::create(
    Rcpp::_("param_names") = std::move(param_names_out),
    Rcpp::_("VA_param_names") = std::move(va_param_names_out));
}

// TODO: need to test this function
/// returns indices of the parameters
// [[Rcpp::export(rng = false)]]
Rcpp::List joint_ms_parameter_indices(SEXP ptr){
  Rcpp::XPtr<problem_data> obj(ptr);

  auto &params = obj->params();
  Rcpp::List m_out(params.marker_info().size()),
             s_out(params.surv_info().size());

  // handle the fixed effects
  for(size_t i = 0; i < params.marker_info().size(); ++i){
    auto &info = params.marker_info()[i];
    Rcpp::IntegerVector fixef(info.n_fix),
                   fixef_vary(info.n_variying);

    std::iota(fixef.begin(), fixef.end(), info.idx_fix + 1);
    std::iota(fixef_vary.begin(), fixef_vary.end(), info.idx_varying + 1);

    m_out[i] = Rcpp::List::create(
      Rcpp::_("fixef") = std::move(fixef),
      Rcpp::_("fixef_vary") = std::move(fixef_vary));
  }

  for(size_t i = 0; i < params.surv_info().size(); ++i){
    auto &info = params.surv_info()[i];

    Rcpp::IntegerVector fixef(info.n_fix),
                   fixef_vary(info.n_variying),
                 associations(params.marker_info().size());

    std::iota(fixef.begin(), fixef.end(), info.idx_fix + 1);
    std::iota(fixef_vary.begin(), fixef_vary.end(), info.idx_varying + 1);
    std::iota
      (associations.begin(), associations.end(), info.idx_association + 1);

    s_out[i] = Rcpp::List::create(
      Rcpp::_("fixef") = std::move(fixef),
      Rcpp::_("fixef_vary") = std::move(fixef_vary),
      Rcpp::_("associations") = std::move(associations));
  }

  // handle the covariance matrices
  vajoint_uint const n_vcov_marker{dim_tri(params.marker_info().size())};
  Rcpp::IntegerVector vcov_marker(n_vcov_marker);
  std::iota
    (vcov_marker.begin(), vcov_marker.end(), params.vcov_marker<true>() + 1);

  vajoint_uint const n_vcov_surv{dim_tri(params.surv_info().size())};
  Rcpp::IntegerVector vcov_surv(n_vcov_surv);
  std::iota
    (vcov_surv.begin(), vcov_surv.end(), params.vcov_surv<true>() + 1);

  vajoint_uint const n_vcov_vary{dim_tri(params.n_shared())};
  Rcpp::IntegerVector vcov_vary(n_vcov_vary);
  std::iota
    (vcov_vary.begin(), vcov_vary.end(), params.vcov_vary<true>() + 1);

  Rcpp::List vcovs = Rcpp::List::create(
    Rcpp::_("vcov_marker") = std::move(vcov_marker),
    Rcpp::_("vcov_surv") = std::move(vcov_surv),
    Rcpp::_("vcov_vary") = std::move(vcov_vary));

  return Rcpp::List::create(
    Rcpp::_("markers") = std::move(m_out),
    Rcpp::_("survival") = std::move(s_out),
    Rcpp::_("vcovs") = std::move(vcovs),
    Rcpp::_("va_params_start") = params.va_mean<true>(),
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
Rcpp::NumericVector opt_priv
  (Rcpp::NumericVector val, SEXP ptr,
   double const rel_eps, unsigned const max_it, unsigned const n_threads,
   double const c1, double const c2){
  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  Rcpp::NumericVector par = clone(val);
  obj->set_n_threads(n_threads);
  double const res = obj->optim().optim_priv(&par[0], rel_eps, max_it, c1, c2);
  par.attr("value") = res;
  return par;
}

/// optimizes the lower bound using the psqn package
// [[Rcpp::export(rng = false)]]
Rcpp::List joint_ms_opt_lb
  (Rcpp::NumericVector val, SEXP ptr, double const rel_eps,
   unsigned const max_it, unsigned const n_threads, double const c1,
   double const c2, bool const use_bfgs, unsigned const trace,
   double const cg_tol, bool const strong_wolfe, size_t const max_cg,
   unsigned const pre_method){
  Rcpp::XPtr<problem_data> obj(ptr);
  check_par_length(*obj, val);

  Rcpp::NumericVector par = clone(val);
  obj->set_n_threads(n_threads);
  auto res = obj->optim().optim(&par[0], rel_eps, max_it, c1, c2,
                                use_bfgs, trace, cg_tol, strong_wolfe, max_cg,
                                static_cast<PSQN::precondition>(pre_method));

  Rcpp::NumericVector counts = Rcpp::NumericVector::create(
    res.n_eval, res.n_grad,  res.n_cg);
  counts.names() =
    Rcpp::CharacterVector::create("function", "gradient", "n_cg");

  int const info = static_cast<int>(res.info);
  return Rcpp::List::create(
    Rcpp::_["par"] = par, Rcpp::_["value"] = res.value,
    Rcpp::_["info"] = info, Rcpp::_["counts"] = counts,
    Rcpp::_["convergence"] =  res.info == PSQN::info_code::converged);
}
