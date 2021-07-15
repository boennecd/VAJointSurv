#include <testthat.h>
#include "VA-parameter.h"

context("subset_params works as expected") {
  test_that("Works with both markers and survival outcomes") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    params.add_survival({ 1L, 4L});
    params.add_survival({ 3L, 2L});

    expect_true(params.get_marker_info().size() == 3);
    expect_true(params.get_survival_info().size() == 2);

    expect_true(params.get_fixef_idx_marker(0) == 0);
    expect_true(params.get_varying_idx_marker(0) == 3);

    expect_true(params.get_fixef_idx_marker(1) == 5);
    expect_true(params.get_varying_idx_marker(1) == 6);

    expect_true(params.get_fixef_idx_marker(2) == 7);
    expect_true(params.get_varying_idx_marker(2) == 10);

    expect_true(params.get_fixef_idx_survival(0) == 13);
    expect_true(params.get_varying_idx_survival(0) == 14);
    expect_true(params.get_idx_association_parameter(0) == 18);

    expect_true(params.get_fixef_idx_survival(1) == 21);
    expect_true(params.get_varying_idx_survival(1) == 24);
    expect_true(params.get_idx_association_parameter(1) == 26);

    expect_true(params.get_idx_error_term() == 29);
    expect_true(params.get_idx_shared_effect() == 38);
    expect_true(params.get_idx_shared_survival() == 74);

    expect_true(params.get_n_params() == 78);
    expect_true(params.get_n_shared() == 6);
    expect_true(params.get_n_shared_survival() == 2);

    expect_true(params.get_idx_va_mean() == 78);
    expect_true(params.get_idx_va_vcov() == 86);

    expect_true(params.get_n_parms_w_va() == 150);

    // with the triangular matrices
    expect_true(params.get_fixef_idx_marker<true>(0) == 0);
    expect_true(params.get_varying_idx_marker<true>(0) == 3);

    expect_true(params.get_fixef_idx_marker<true>(1) == 5);
    expect_true(params.get_varying_idx_marker<true>(1) == 6);

    expect_true(params.get_fixef_idx_marker<true>(2) == 7);
    expect_true(params.get_varying_idx_marker<true>(2) == 10);

    expect_true(params.get_fixef_idx_survival<true>(0) == 13);
    expect_true(params.get_varying_idx_survival<true>(0) == 14);
    expect_true(params.get_idx_association_parameter<true>(0) == 18);

    expect_true(params.get_fixef_idx_survival<true>(1) == 21);
    expect_true(params.get_varying_idx_survival<true>(1) == 24);
    expect_true(params.get_idx_association_parameter<true>(1) == 26);

    expect_true(params.get_idx_error_term<true>() == 29);
    expect_true(params.get_idx_shared_effect<true>() == 35);
    expect_true(params.get_idx_shared_survival<true>() == 56);

    expect_true(params.get_n_params<true>() == 59);

    expect_true(params.get_idx_va_mean<true>() == 59);
    expect_true(params.get_idx_va_vcov<true>() == 67);

    expect_true(params.get_n_parms_w_va<true>() == 103);
  }

  test_that("Works only with markers") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    expect_true(params.get_marker_info().size() == 3);
    expect_true(params.get_survival_info().size() == 0);

    expect_true(params.get_fixef_idx_marker(0) == 0);
    expect_true(params.get_varying_idx_marker(0) == 3);

    expect_true(params.get_fixef_idx_marker(1) == 5);
    expect_true(params.get_varying_idx_marker(1) == 6);

    expect_true(params.get_fixef_idx_marker(2) == 7);
    expect_true(params.get_varying_idx_marker(2) == 10);

    expect_true(params.get_idx_error_term() == 13);
    expect_true(params.get_idx_shared_effect() == 22);

    expect_true(params.get_idx_shared_survival() == 58);
    expect_true(params.get_n_params() == 58);
    expect_true(params.get_n_shared() == 6);
    expect_true(params.get_n_shared_survival() == 0);

    expect_true(params.get_idx_va_mean() == 58);
    expect_true(params.get_idx_va_vcov() == 64);

    expect_true(params.get_n_parms_w_va() == 100);

    // with the triangular matrices
    expect_true(params.get_fixef_idx_marker<true>(0) == 0);
    expect_true(params.get_varying_idx_marker<true>(0) == 3);

    expect_true(params.get_fixef_idx_marker<true>(1) == 5);
    expect_true(params.get_varying_idx_marker<true>(1) == 6);

    expect_true(params.get_fixef_idx_marker<true>(2) == 7);
    expect_true(params.get_varying_idx_marker<true>(2) == 10);

    expect_true(params.get_idx_error_term<true>() == 13);
    expect_true(params.get_idx_shared_effect<true>() == 19);

    expect_true(params.get_idx_shared_survival<true>() == 40);
    expect_true(params.get_n_params<true>() == 40);

    expect_true(params.get_idx_va_mean<true>() == 40);
    expect_true(params.get_idx_va_vcov<true>() == 46);

    expect_true(params.get_n_parms_w_va<true>() == 67);
  }

  test_that("Works only with survival outcomes") {
    subset_params params;

    params.add_survival({ 1L, 4L});
    params.add_survival({ 3L, 2L});

    expect_true(params.get_marker_info().size() == 0);
    expect_true(params.get_survival_info().size() == 2);

    expect_true(params.get_fixef_idx_survival(0) == 0);
    expect_true(params.get_varying_idx_survival(0) == 1);
    expect_true(params.get_idx_association_parameter(0) == 5);

    expect_true(params.get_fixef_idx_survival(1) == 5);
    expect_true(params.get_varying_idx_survival(1) == 8);
    expect_true(params.get_idx_association_parameter(1) == 10);

    expect_true(params.get_idx_error_term() == 10);
    expect_true(params.get_idx_shared_effect() == 10);
    expect_true(params.get_idx_shared_survival() == 10);
    expect_true(params.get_n_params() == 14);
    expect_true(params.get_n_shared() == 0);
    expect_true(params.get_n_shared_survival() == 2);

    expect_true(params.get_idx_va_mean() == 14);
    expect_true(params.get_idx_va_vcov() == 16);

    expect_true(params.get_n_parms_w_va() == 20);

    // with the triangular matrices
    expect_true(params.get_fixef_idx_survival<true>(0) == 0);
    expect_true(params.get_varying_idx_survival<true>(0) == 1);
    expect_true(params.get_idx_association_parameter<true>(0) == 5);

    expect_true(params.get_fixef_idx_survival<true>(1) == 5);
    expect_true(params.get_varying_idx_survival<true>(1) == 8);
    expect_true(params.get_idx_association_parameter<true>(1) == 10);

    expect_true(params.get_idx_error_term<true>() == 10);
    expect_true(params.get_idx_shared_effect<true>() == 10);
    expect_true(params.get_idx_shared_survival<true>() == 10);
    expect_true(params.get_n_params<true>() == 13);

    expect_true(params.get_idx_va_mean<true>() == 13);
    expect_true(params.get_idx_va_vcov<true>() == 15);

    expect_true(params.get_n_parms_w_va<true>() == 18);
  }
}
