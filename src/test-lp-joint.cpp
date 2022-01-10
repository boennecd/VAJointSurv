#include "testthat-wrapper.h"
#include "lp-joint.h"

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
    expect_true(pass_rel_err(val, 25.2257982961411));
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
      expect_true(pass_rel_err(val, 17.6863987933819));
    }
    {
      double const val = lp_joint::submat_trace(A, X, k, n, 2);
      expect_true(pass_rel_err(val, 35.2672271852575));
    }
  }

  test_that("mat_add works as expected") {
    /*
     k <- 3L
     n <- 5L
     set.seed(1)
     M <- drop(rWishart(1, k, diag(k)))
     A <- drop(rWishart(1, n, diag(n)))
     dput(c(A))
     dput(c(M))

     V <- A
     V[1:k, 1:k] <- V[1:k, 1:k] +  2 * M
     dput(c(V))

     V <- A
     V[1:k + 2, 1:k + 2] <- V[1:k + 2, 1:k + 2] + M
     dput(c(V))

     V <- A
     V[1:k + 2, 1:k + 2] <- V[1:k + 2, 1:k + 2] + 2 * M
     dput(c(V))
     */
    constexpr vajoint_uint k = 3,
                           n = 5;
    constexpr double B[k * k] = { 0.94331456701213, 1.23583912080175, -0.901864998282764, 1.23583912080175, 5.38387957072665, -1.75338497451975, -0.901864998282764, -1.75338497451975, 0.959799081627646 },
                     A[n * n] = { 3.98370460230852, -1.59476013270181, -4.42036821280508, -1.78020499400603, 0.753252269893148, -1.59476013270181, 2.93700109157094, 3.475083142122, 1.37319643063454, -0.0993902072615885, -4.42036821280508, 3.475083142122, 6.53360944781029, 1.71958991015682, -0.201140885765227, -1.78020499400603, 1.37319643063454, 1.71958991015682, 3.45430837572081, -1.32902247007167, 0.753252269893148, -0.0993902072615885, -0.201140885765227, -1.32902247007167, 1.28437823319705 };

    double res[n * n];
    {
      constexpr double val[n * n] = { 5.87033373633278, 0.876918108901697, -6.22409820937061, -1.78020499400603,
                                      0.753252269893148, 0.876918108901697, 13.7047602330243, -0.0316868069175014,
                                      1.37319643063454, -0.0993902072615885, -6.22409820937061, -0.0316868069175014,
                                      8.45320761106558, 1.71958991015682, -0.201140885765227, -1.78020499400603,
                                      1.37319643063454, 1.71958991015682, 3.45430837572081, -1.32902247007167,
                                      0.753252269893148, -0.0993902072615885, -0.201140885765227, -1.32902247007167,
                                      1.28437823319705 };
      std::copy(A, A + n * n, res);
      lp_joint::mat_add(res, B, k, n, 0, 2.);

      for(vajoint_uint i = 0; i < n * n; ++i)
        expect_true(pass_rel_err(res[i], val[i]));

    }
    {
      std::copy(A, A + n * n, res);
      lp_joint::mat_add(res, B, k, n, 2, 1.);

      constexpr double val[n * n] = { 3.98370460230852, -1.59476013270181, -4.42036821280508, -1.78020499400603,
                                      0.753252269893148, -1.59476013270181, 2.93700109157094, 3.475083142122,
                                      1.37319643063454, -0.0993902072615885, -4.42036821280508, 3.475083142122,
                                      7.47692401482242, 2.95542903095857, -1.10300588404799, -1.78020499400603,
                                      1.37319643063454, 2.95542903095857, 8.83818794644746, -3.08240744459143,
                                      0.753252269893148, -0.0993902072615885, -1.10300588404799, -3.08240744459143,
                                      2.2441773148247 };

      for(vajoint_uint i = 0; i < n * n; ++i)
        expect_true(pass_rel_err(res[i], val[i]));

      constexpr double val2[n * n] = { 3.98370460230852, -1.59476013270181, -4.42036821280508, -1.78020499400603,
                                       0.753252269893148, -1.59476013270181, 2.93700109157094, 3.475083142122,
                                       1.37319643063454, -0.0993902072615885, -4.42036821280508, 3.475083142122,
                                       8.42023858183455, 4.19126815176032, -2.00487088233075, -1.78020499400603,
                                       1.37319643063454, 4.19126815176032, 14.2220675171741, -4.83579241911118,
                                       0.753252269893148, -0.0993902072615885, -2.00487088233075, -4.83579241911118,
                                       3.20397639645234 };

      lp_joint::mat_add(res, B, k, n, 2, 1.);
      for(vajoint_uint i = 0; i < n * n; ++i)
        expect_true(pass_rel_err(res[i], val2[i]));

    }
  }


  test_that("mat_vec_prod works as expected") {
    /*
     k <- 3L
     set.seed(2)
     b <- rnorm(k)
     M <- drop(rWishart(1, k, diag(k)))
     dput(c(b))
     dput(c(M))
     dput(drop(M %*% b))
     */
    constexpr vajoint_uint k = 3;
    constexpr double b[k] = { -0.896914546624981, 0.184849184646742, 1.58784533120882 },
                 M[k * k] = { 0.378123233971301, -0.436506592144316, 0.646754211145205, -0.436506592144316,
                              1.68944864054823, -1.56613751629636, 0.646754211145205, -1.56613751629636,
                              1.7643953134308 };

    double res[k];

    std::fill(res, res + k, 0);
    lp_joint::mat_vec_prod(M, b, res, k);

    constexpr double val[k] = { 0.607113537990013, -1.78298182738233, 1.93200435783866 };
    for(vajoint_uint i = 0; i < k; ++i)
      expect_true(pass_rel_err(res[i], val[i]));

    lp_joint::mat_vec_prod(M, b, res, k);
    for(vajoint_uint i = 0; i < k; ++i)
      expect_true(pass_rel_err(res[i], 2 * val[i]));
  }

  test_that("rank_one works as expected") {
    /*
     k <- 3L
     set.seed(2)
     b <- rnorm(k)
     M <- drop(rWishart(1, k, diag(k)))
     dput(c(b))
     dput(c(M))

     dput(drop(M + outer(b, b)))
     */
    constexpr vajoint_uint k = 3;
    constexpr double M[k * k] = { 0.378123233971301, -0.436506592144316, 0.646754211145205, -0.436506592144316,
                                  1.68944864054823, -1.56613751629636, 0.646754211145205, -1.56613751629636,
                                  1.7643953134308 },
                     b[k] = { -0.896914546624981, 0.184849184646742, 1.58784533120882 },
                   val[k * k] = { 1.1825789379188, -0.602300514785747, -0.77740736420655,
                                  -0.602300514785747, 1.72361786161279, -1.27262560147727, -0.77740736420655,
                                  -1.27262560147727, 4.28564810927245 };

    double res[k * k];
    std::copy(M, M + k * k, res);
    lp_joint::rank_one<true>(res, b, k);
    for(vajoint_uint i = 0; i < k * k; ++i)
      expect_true(pass_rel_err(res[i], val[i]));

    lp_joint::rank_one<false>(res, b, k);
    for(vajoint_uint i = 0; i < k * k; ++i)
      expect_true(pass_rel_err(res[i], M[i]));
  }

  test_that("copy_block_upper_tri works as expected"){
    /*
     n <- 3L
     k <- 4L
     set.seed(1)
     A <- round(drop(rWishart(1, n, diag(n))), 2)
     B <- round(drop(rWishart(1, k, diag(k))), 2)
     truth <- matrix(0, n + k, n + k)
     truth[1:n, 1:n] <- A
     truth[1:k + n, 1:k + n] <- B
     dput(A[upper.tri(A, TRUE)])
     dput(B[upper.tri(B, TRUE)])
     dput(truth[upper.tri(truth, TRUE)])
     */

    constexpr size_t na{3}, nb{4}, n_all{na + nb},
                dim_res{(n_all * (n_all + 1)) / 2};
    constexpr double A[]{0.94, 1.24, 5.38, -0.9, -1.75, 0.96},
                     B[]{2.99, -1.38, 2.07, -3.83, 3.12, 6.21, 1.63, 0.23, -1.05, 2.6},
                 truth[]{0.94, 1.24, 5.38, -0.9, -1.75, 0.96, 0, 0, 0, 2.99, 0, 0, 0, -1.38, 2.07, 0, 0, 0, -3.83, 3.12, 6.21, 0, 0, 0, 1.63, 0.23, -1.05, 2.6};

    double res[dim_res];
    lp_joint::copy_block_upper_tri(res, A, B, na, nb);
    for(size_t i = 0; i < dim_res; ++i)
      expect_true(res[i] == truth[i]);
  }
}
