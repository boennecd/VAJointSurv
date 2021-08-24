#include "testthat-wrapper.h"
#include "log-cholesky.h"
#include <algorithm>
#include <iterator>

context("log-cholesky works as expected") {
  test_that("log_chol::pd_mat works as expected") {
    /*
     set.seed(1)
     n <- 4
     X <- drop(rWishart(1, n, diag(n)))
     M <- chol(X)
     diag(M) <- log(diag(M))
     dput(M[upper.tri(M, TRUE)])
     dput(X)
     */
    constexpr vajoint_uint dim = 4;
    constexpr double X[dim * dim] { 1.66173000475602, 1.64026455103121, -1.98512044759873,  -0.00743435321565165, 1.64026455103121, 7.16285794371291, -4.14581062030046, 5.65447492748679, -1.98512044759873, -4.14581062030046, 4.90603816923005, -1.23652487489528, -0.00743435321565165, 5.65447492748679, -1.23652487489528, 6.77907155498759 };
    constexpr double theta[dim_tri(dim)] { 0.253929615612238, 1.2724293214294, 0.856338430490728, -1.53995004190371, -0.928567034713538, 0.25711649595832, -0.00576717274753696, 2.40465338885795, 0.763593461140459, -0.441421449571498 };
    double res[dim * dim];

    log_chol::pd_mat::get(theta, dim, res);
    for(vajoint_uint i = 0; i < dim * dim; ++i)
      expect_true(pass_rel_err(res[i], X[i]));

    // with working memory supplied
    std::unique_ptr<double[]>
      mem(new double[log_chol::pd_mat::n_wmem(dim)]);

    std::fill(std::begin(res), std::end(res), 0);
    log_chol::pd_mat::get(theta, dim, res, mem.get());
    for(vajoint_uint i = 0; i < dim * dim; ++i)
      expect_true(pass_rel_err(res[i], X[i]));

    // clean up
    wmem::clear_all();
  }

  test_that("log_chol::dpd_mat works as expected") {
    /*
     set.seed(1)
     options(width = 1000)
     n <- 4
     X <- drop(rWishart(1, n, diag(n)))
     M <- chol(X)
     diag(M) <- log(diag(M))
     dput(M[upper.tri(M, TRUE)])

     f <- function(x){
     L <- matrix(0, n, n)
     L[upper.tri(L, TRUE)] <- x
     diag(L) <- exp(diag(L))
     sum(exp(crossprod(L)))
     }

     library(numDeriv)
     deriv <- jacobian(f, M[upper.tri(M, TRUE)])
     dput(deriv)
     dput(c(exp(X)))
     tmp <- exp(X)
     tmp[lower.tri(tmp)] <- 0
     dput(tmp)
     */
    constexpr vajoint_uint dim = 4,
                      dim_ltri = dim_tri(dim);
    constexpr double res[dim_ltri]  = { 33.8654055763531, 3294.33138190854, 17543.1228663122, -415.712956102195, -249.433435111634, 452.454270319483, 718.249462770731, 5572.7994952304, 1343.53323421608, 727.32675735327 };

    constexpr double theta[dim_ltri] = { 0.253929615612238, 1.2724293214294, 0.856338430490728, -1.53995004190371, -0.928567034713538, 0.25711649595832, -0.00576717274753696, 2.40465338885795, 0.763593461140459, -0.441421449571498 };
    constexpr double derivs[dim * dim] = { 5.26841735210069, 5.15653349805919, 0.137364067948365, 0.99259321323301, 5.15653349805919, 1290.59411256876, 0.0158305981445064, 285.566500234151, 0.137364067948365, 0.0158305981445064, 135.103097103552, 0.29039161369905, 0.99259321323301, 285.566500234151, 0.29039161369905, 879.252007887043 };
    constexpr double derivs_half[dim * dim] = { 5.26841735210069, 0, 0, 0, 5.15653349805919, 1290.59411256876, 0, 0, 0.137364067948365, 0.0158305981445064, 135.103097103552, 0, 0.99259321323301, 285.566500234151, 0.29039161369905, 879.252007887043 };
    double output[dim_ltri];

    std::fill(std::begin(output), std::end(output), 0);
    log_chol::dpd_mat::get(theta, dim, output, derivs);
    double const eps = std::sqrt(std::numeric_limits<double>::epsilon());
    for(vajoint_uint i = 0; i < dim_ltri; ++i)
      expect_true(pass_rel_err(output[i], res[i]));

    // works correctly if add a value to output
    std::fill(std::begin(output), std::end(output), 100);
    log_chol::dpd_mat::get(theta, dim, output, derivs);
    for(vajoint_uint i = 0; i < dim_ltri; ++i)
      expect_true(pass_rel_err(output[i] - 100, res[i]));

    // works when only the upper triangular matrix has elements
    std::fill(std::begin(output), std::end(output), 0);
    log_chol::dpd_mat::get(theta, dim, output, derivs_half);
    for(vajoint_uint i = 0; i < dim_ltri; ++i)
      expect_true(pass_rel_err(output[i], res[i]));

    // works with pre-allocated memory
    std::unique_ptr<double[]> mem(new double[log_chol::dpd_mat::n_wmem(dim)]);
    std::fill(std::begin(output), std::end(output), 0);
    log_chol::dpd_mat::get(theta, dim, output, derivs, mem.get());
    for(vajoint_uint i = 0; i < dim_ltri; ++i)
      expect_true(pass_rel_err(output[i], res[i]));

    // clean up
    wmem::clear_all();
  }
}
