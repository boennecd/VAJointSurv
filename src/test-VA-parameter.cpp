#include "testthat-wrapper.h"
#include "VA-parameter.h"

context("subset_params works as expected") {
  test_that("Works with both markers and survival outcomes") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    params.add_surv({ 1L, 4L});
    params.add_surv({ 3L, 2L});

    expect_true(params.marker_info().size() == 3);
    expect_true(params.surv_info().size() == 2);

    expect_true(params.fixef_marker(0) == 0);
    expect_true(params.fixef_marker(1) == 3);
    expect_true(params.fixef_marker(2) == 4);

    expect_true(params.fixef_vary_marker(0) == 7);
    expect_true(params.fixef_vary_marker(1) == 9);
    expect_true(params.fixef_vary_marker(2) == 10);

    expect_true(params.fixef_surv(0) == 13);
    expect_true(params.fixef_vary_surv(0) == 14);
    expect_true(params.association(0) == 18);

    expect_true(params.fixef_surv(1) == 21);
    expect_true(params.fixef_vary_surv(1) == 24);
    expect_true(params.association(1) == 26);

    expect_true(params.vcov_marker() == 29);
    expect_true(params.vcov_vary() == 38);
    expect_true(params.vcov_surv() == 74);

    expect_true(params.n_params() == 78);
    expect_true(params.n_shared() == 6);
    expect_true(params.n_shared_surv() == 2);

    expect_true(params.va_mean() == 78);
    expect_true(params.va_vcov() == 86);

    expect_true(params.n_parms_w_va() == 150);

    // with the triangular matrices
    expect_true(params.fixef_marker<true>(0) == 0);
    expect_true(params.fixef_marker<true>(1) == 3);
    expect_true(params.fixef_marker<true>(2) == 4);

    expect_true(params.fixef_vary_marker<true>(0) == 7);
    expect_true(params.fixef_vary_marker<true>(1) == 9);
    expect_true(params.fixef_vary_marker<true>(2) == 10);

    expect_true(params.fixef_surv<true>(0) == 13);
    expect_true(params.fixef_vary_surv<true>(0) == 14);
    expect_true(params.association<true>(0) == 18);

    expect_true(params.fixef_surv<true>(1) == 21);
    expect_true(params.fixef_vary_surv<true>(1) == 24);
    expect_true(params.association<true>(1) == 26);

    expect_true(params.vcov_marker<true>() == 29);
    expect_true(params.vcov_vary<true>() == 35);
    expect_true(params.vcov_surv<true>() == 56);

    expect_true(params.n_params<true>() == 59);

    expect_true(params.va_mean<true>() == 59);
    expect_true(params.va_vcov<true>() == 67);

    expect_true(params.n_parms_w_va<true>() == 103);
  }

  test_that("Works only with markers") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    expect_true(params.fixef_marker(0) == 0);
    expect_true(params.fixef_marker(1) == 3);
    expect_true(params.fixef_marker(2) == 4);

    expect_true(params.fixef_vary_marker(0) == 7);
    expect_true(params.fixef_vary_marker(1) == 9);
    expect_true(params.fixef_vary_marker(2) == 10);

    expect_true(params.vcov_marker() == 13);
    expect_true(params.vcov_vary() == 22);

    expect_true(params.vcov_surv() == 58);
    expect_true(params.n_params() == 58);
    expect_true(params.n_shared() == 6);
    expect_true(params.n_shared_surv() == 0);

    expect_true(params.va_mean() == 58);
    expect_true(params.va_vcov() == 64);

    expect_true(params.n_parms_w_va() == 100);

    // with the triangular matrices
    expect_true(params.fixef_marker<true>(0) == 0);
    expect_true(params.fixef_marker<true>(1) == 3);
    expect_true(params.fixef_marker<true>(2) == 4);

    expect_true(params.fixef_vary_marker<true>(0) == 7);
    expect_true(params.fixef_vary_marker<true>(1) == 9);
    expect_true(params.fixef_vary_marker<true>(2) == 10);

    expect_true(params.vcov_marker<true>() == 13);
    expect_true(params.vcov_vary<true>() == 19);

    expect_true(params.vcov_surv<true>() == 40);
    expect_true(params.n_params<true>() == 40);

    expect_true(params.va_mean<true>() == 40);
    expect_true(params.va_vcov<true>() == 46);

    expect_true(params.n_parms_w_va<true>() == 67);
  }

  test_that("Works only with survival outcomes") {
    subset_params params;

    params.add_surv({ 1L, 4L});
    params.add_surv({ 3L, 2L});

    expect_true(params.marker_info().size() == 0);
    expect_true(params.surv_info().size() == 2);

    expect_true(params.fixef_surv(0) == 0);
    expect_true(params.fixef_vary_surv(0) == 1);
    expect_true(params.association(0) == 5);

    expect_true(params.fixef_surv(1) == 5);
    expect_true(params.fixef_vary_surv(1) == 8);
    expect_true(params.association(1) == 10);

    expect_true(params.vcov_marker() == 10);
    expect_true(params.vcov_vary() == 10);
    expect_true(params.vcov_surv() == 10);
    expect_true(params.n_params() == 14);
    expect_true(params.n_shared() == 0);
    expect_true(params.n_shared_surv() == 2);

    expect_true(params.va_mean() == 14);
    expect_true(params.va_vcov() == 16);

    expect_true(params.n_parms_w_va() == 20);

    // with the triangular matrices
    expect_true(params.fixef_surv<true>(0) == 0);
    expect_true(params.fixef_vary_surv<true>(0) == 1);
    expect_true(params.association<true>(0) == 5);

    expect_true(params.fixef_surv<true>(1) == 5);
    expect_true(params.fixef_vary_surv<true>(1) == 8);
    expect_true(params.association<true>(1) == 10);

    expect_true(params.vcov_marker<true>() == 10);
    expect_true(params.vcov_vary<true>() == 10);
    expect_true(params.vcov_surv<true>() == 10);
    expect_true(params.n_params<true>() == 13);

    expect_true(params.va_mean<true>() == 13);
    expect_true(params.va_vcov<true>() == 15);

    expect_true(params.n_parms_w_va<true>() == 18);
  }
}
