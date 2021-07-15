#include "lp-joint.h"
#include <limits>
#include <testthat.h>

context("testing lp_joint functions") {
  test_that("quad_form works as expected") {
    /*
     set.seed(1)
     n <- 3
     X <- drop(rWishart(1, n, diag(n)))
     M <- chol(X)
     v <- rnorm(n)
     dput(v)
     dput(c(M))
     dput(v %*% X %*% v)
     */
    constexpr vajoint_uint dim = 3;
    constexpr double M[dim * dim] = { 0.971243824697038, 0, 0, 1.2724293214294, 1.94031007642937, 0, -0.928567034713538, -0.29472044679056, 0.103451447135964 };
    constexpr double x[dim] = { -0.00576717274753696, 2.40465338885795, 0.76359346114046 };
    double mem[dim];

    double const val = lp_joint::quad_form(x, M, dim, mem);
    expect_true(std::abs(val - 25.2257982961411) < 1e-8);
  }

  test_that("submat_trace works as expected") {
    /*
     k <- 3L
     n <- 5L
     set.seed(1)
     A <- drop(rWishart(1, k, diag(k)))
     M <- drop(rWishart(1, n, diag(n)))
     dput(c(A))
     dput(c(M))
     dput(sum(diag(A %*% M[1:k, 1:k])))
     dput(sum(diag(A %*% M[1:k + 2, 1:k + 2])))
     */
    constexpr vajoint_uint k = 3,
                           n = 5;
    constexpr double A[k * k] = { 0.94331456701213, 1.23583912080175, -0.901864998282764, 1.23583912080175, 5.38387957072665, -1.75338497451975, -0.901864998282764, -1.75338497451975, 0.959799081627646 },
                     X[n * n] = { 3.98370460230852, -1.59476013270181, -4.42036821280508, -1.78020499400603, 0.753252269893148, -1.59476013270181, 2.93700109157094, 3.475083142122, 1.37319643063454, -0.0993902072615885, -4.42036821280508, 3.475083142122, 6.53360944781029, 1.71958991015682, -0.201140885765227, -1.78020499400603, 1.37319643063454, 1.71958991015682, 3.45430837572081, -1.32902247007167, 0.753252269893148, -0.0993902072615885, -0.201140885765227, -1.32902247007167, 1.28437823319705 };

    {
      double const val = lp_joint::submat_trace(A, X, k, n, 0);
      expect_true(std::abs(val - 17.6863987933819) < 1e-8);
    }
    {
      double const val = lp_joint::submat_trace(A, X, k, n, 2);
      expect_true(std::abs(val - 35.2672271852575) < 1e-8);
    }
  }
}
