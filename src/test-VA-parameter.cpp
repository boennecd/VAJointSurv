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

    expect_true(params.get_idx_error_term() == 13);
    expect_true(params.get_idx_shared_effect() == 22);

    expect_true(params.get_fixef_idx_survival(0) == 58);
    expect_true(params.get_varying_idx_survival(0) == 59);
    expect_true(params.get_idx_association_parameter(0) == 63);

    expect_true(params.get_fixef_idx_survival(1) == 66);
    expect_true(params.get_varying_idx_survival(1) == 69);
    expect_true(params.get_idx_association_parameter(1) == 71);

    expect_true(params.get_idx_shared_survival() == 74);
    expect_true(params.get_n_params() == 78);
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
  }

  test_that("Works only with survival outcomes") {
    subset_params params;

    params.add_survival({ 1L, 4L});
    params.add_survival({ 3L, 2L});

    expect_true(params.get_marker_info().size() == 0);
    expect_true(params.get_survival_info().size() == 2);

    expect_true(params.get_idx_error_term() == 0);
    expect_true(params.get_idx_shared_effect() == 0);

    expect_true(params.get_fixef_idx_survival(0) == 0);
    expect_true(params.get_varying_idx_survival(0) == 1);
    expect_true(params.get_idx_association_parameter(0) == 5);

    expect_true(params.get_fixef_idx_survival(1) == 5);
    expect_true(params.get_varying_idx_survival(1) == 8);
    expect_true(params.get_idx_association_parameter(1) == 10);

    expect_true(params.get_idx_shared_survival() == 10);
    expect_true(params.get_n_params() == 14);
  }
}
