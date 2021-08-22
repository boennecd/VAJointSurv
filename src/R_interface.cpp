#include "arma-wrap.h"
#include "marker-term.h"

// TODO: need to check that expansions are created correctly
/// returns a pointer to an expansions given an R List with data.
// std::unique_ptr<joint_bases::basisMixin> get_basis_from_list
//   (Rcpp::List dat){
//   std::string dat_class; // TODO: get the class
//   if(dat_class == "poly_term"){
//     Rcpp::List coefs{dat["coefs"]};
//     arma::vec alpha{Rcpp::as<arma::vec>(coefs["alpha"])},
//               norm2{Rcpp::as<arma::vec>(coefs["norm2"])};
//     bool const raw{Rcpp::as<bool>(dat["raw"])},
//          intercept{Rcpp::as<bool>(dat["intercept"])};
//
//     if(raw)
//       return std::make_unique<joint_bases::orth_poly>
//         (new joint_bases::orth_poly(alpha.size(), intercept));
//     else
//       return std::make_unique<joint_bases::orth_poly>
//         (new joint_bases::orth_poly(alpha, norm2, intercept));
//
//   } else if(dat_class == "bs_term"){
//     arma::vec i_knots{Rcpp::as<arma::vec>(dat["knots"])},
//               b_knots{Rcpp::as<arma::vec>(dat["Boundary.knots"])};
//
//     bool const intercept{Rcpp::as<bool>(dat["intercept"])};
//     vajoint_uint const degree{Rcpp::as<vajoint_uint>(dat["degree"])};
//
//     return std::make_unique<joint_bases::bs>
//       (new joint_bases::bs(b_knots, i_knots, intercept, degree + 1));
//
//   } else if(dat_class == "ns_term"){
//     arma::vec i_knots{Rcpp::as<arma::vec>(dat["knots"])},
//               b_knots{Rcpp::as<arma::vec>(dat["Boundary.knots"])};
//
//     bool const intercept{Rcpp::as<bool>(dat["intercept"])};
//     vajoint_uint const degree{Rcpp::as<vajoint_uint>(dat["degree"])};
//
//     return std::make_unique<joint_bases::ns>
//       (new joint_bases::ns(b_knots, i_knots, intercept, degree + 1));
//
//   }
//
//   throw std::invalid_argument("expansions is not implemented");
//   return std::make_unique<joint_bases::orth_poly>
//     (new joint_bases::orth_poly(2, false));
// }
//
// struct comp_obj {
//   subset_params par_idx;
//
// public:
//   comp_obj(Rcpp::List markers){
//     // handle the markers
//     std::vector<marker::setup_marker_dat_helper> input_dat;
//     joint_bases::bases_vector bases_fix;
//     joint_bases::bases_vector bases_rng;
//
//     for(auto m : markers){
//       Rcpp::List marker{m};
//       bases_fix.emplace_back(get_basis_from_list(marker["time_fixef"]));
//       bases_rng.emplace_back(get_basis_from_list(marker["time_rng"]));
//
//       Rcpp::NumericMatrix X{marker["X"]};
//       Rcpp::IntegerVector id{marker["id"]};
//       Rcpp::NumericVector y{marker["y"]},
//                        time{marker["time"]};
//
//       vajoint_uint const n_fixef = X.nrow(),
//                          n_obs   = X.ncol();
//
//       input_dat.emplace_back(&X[0], n_fixef, n_obs, &id[0], &time[0], &y[0]);
//       par_idx.add_marker
//         ({n_fixef, bases_fix.back()->n_basis(), bases_rng.back()->n_basis()});
//     }
//   }
//
// };


// Rcpp::XPtr<comp_obj> joint_ms_ptr(Rcpp::List markers){
//
//
// }
