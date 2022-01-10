#include "testthat-wrapper.h"
#include "VA-parameter.h"

context("subset_params works as expected") {
  test_that("Works with both markers and survival outcomes") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    params.add_surv({ 1L, 4L, {1, 1, 1}, true});
    params.add_surv({ 3L, 2L, {1, 1, 1}, true});

    expect_true(params.marker_info().size() == 3);
    expect_true(params.surv_info().size() == 2);
    for(auto &info : params.surv_info()){
      expect_true(info.n_associations[0] == 1);
      expect_true(info.n_associations[1] == 1);
      expect_true(info.n_associations[2] == 1);
    }

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

    expect_true(params.vcov_start() == 29);
    expect_true(params.vcov_marker() == 29);
    expect_true(params.vcov_vary() == 38);
    expect_true(params.vcov_surv() == 74);
    expect_true(params.vcov_end() == 78);

    expect_true(params.n_params() == 78);
    expect_true(params.n_shared() == 6);
    expect_true(params.n_shared_surv() == 2);

    expect_true(params.va_par_start() == params.va_mean());
    expect_true(params.va_mean() == 78);
    expect_true(params.va_mean_end() == 86);
    expect_true(params.va_vcov() == 86);
    expect_true(params.va_vcov_end() == 150);
    expect_true(params.va_par_end() == params.va_vcov_end());

    expect_true(params.n_params_w_va() == 150);
    expect_true(params.n_va_params() == 72);

    expect_true(params.frailty_offset(0) == 0);
    expect_true(params.frailty_offset(1) == 1);

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

    expect_true(params.vcov_start<true>() == 29);
    expect_true(params.vcov_marker<true>() == 29);
    expect_true(params.vcov_vary<true>() == 35);
    expect_true(params.vcov_surv<true>() == 56);
    expect_true(params.vcov_end<true>() == 59);

    expect_true(params.n_params<true>() == 59);

    expect_true(params.va_par_start<true>() == params.va_mean<true>());
    expect_true(params.va_mean<true>() == 59);
    expect_true(params.va_mean_end<true>() == 67);
    expect_true(params.va_vcov<true>() == 67);
    expect_true(params.va_vcov_end<true>() == 103);
    expect_true(params.va_par_end<true>() == params.va_vcov_end<true>());

    expect_true(params.n_params_w_va<true>() == 103);
    expect_true(params.n_va_params<true>() == 44);

    std::vector<std::string> const param_names{params.param_names(true)};
    expect_true(param_names.size() == 59);

    std::vector<std::string> const true_names
      {
        // the fixed effects
        "mark1_fixef1", "mark1_fixef2", "mark1_fixef3",
        "mark2_fixef1",
        "mark3_fixef1", "mark3_fixef2", "mark3_fixef3",
        // time-varying fixed effects
        "mark1_fixef_vary1", "mark1_fixef_vary2",
        "mark2_fixef_vary1",
        "mark3_fixef_vary1", "mark3_fixef_vary2", "mark3_fixef_vary3",
        // parameters for the first survival outcome
        "surv1_fixef1", "surv1_fixef_vary1", "surv1_fixef_vary2",  "surv1_fixef_vary3",
        "surv1_fixef_vary4", "surv1_assoc1_1", "surv1_assoc2_1", "surv1_assoc3_1",
        // parameters for the second survival outcome
        "surv2_fixef1", "surv2_fixef2", "surv2_fixef3", "surv2_fixef_vary1",
        "surv2_fixef_vary2", "surv2_assoc1_1",  "surv2_assoc2_1",  "surv2_assoc3_1",
        // vcov for the markers' error term
        "vcov_marker1", "vcov_marker2", "vcov_marker3", "vcov_marker4", "vcov_marker5", "vcov_marker6",
        // vcov for the shared effect
        "vcov_vary1", "vcov_vary2", "vcov_vary3", "vcov_vary4", "vcov_vary5", "vcov_vary6", "vcov_vary7", "vcov_vary8", "vcov_vary9", "vcov_vary10", "vcov_vary11", "vcov_vary12", "vcov_vary13", "vcov_vary14", "vcov_vary15", "vcov_vary16", "vcov_vary17", "vcov_vary18", "vcov_vary19", "vcov_vary20", "vcov_vary21",
        // vcov for the frailties
        "vcov_surv1", "vcov_surv2", "vcov_surv3",
      };

    for(size_t i = 0; i < true_names.size(); ++i)
      expect_true(param_names[i] == true_names[i]);

    std::vector<std::string> const va_param_names{params.va_param_names(true)};
    expect_true(va_param_names.size() == 44);

    std::vector<std::string> const va_true_names
    {
      "mark1_rng1", "mark1_rng2",
      "mark2_rng1", "mark2_rng2", "mark2_rng3",
      "mark3_rng1",
      "frailty1", "frailty2",
      "VA_vcov1", "VA_vcov2", "VA_vcov3", "VA_vcov4", "VA_vcov5", "VA_vcov6", "VA_vcov7", "VA_vcov8", "VA_vcov9", "VA_vcov10", "VA_vcov11", "VA_vcov12", "VA_vcov13", "VA_vcov14", "VA_vcov15", "VA_vcov16", "VA_vcov17", "VA_vcov18", "VA_vcov19", "VA_vcov20", "VA_vcov21", "VA_vcov22", "VA_vcov23", "VA_vcov24", "VA_vcov25", "VA_vcov26", "VA_vcov27", "VA_vcov28", "VA_vcov29", "VA_vcov30", "VA_vcov31", "VA_vcov32", "VA_vcov33", "VA_vcov34", "VA_vcov35", "VA_vcov36"
    };

    for(size_t i = 0; i < va_true_names.size(); ++i)
      expect_true(va_param_names[i] == va_true_names[i]);
  }

  test_that("Works with both markers and survival outcomes without one frailty") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    params.add_surv({ 1L, 4L, {1, 1, 1}, false});
    params.add_surv({ 3L, 2L, {1, 1, 1}, true});

    expect_true(params.marker_info().size() == 3);
    expect_true(params.surv_info().size() == 2);
    for(auto &info : params.surv_info()){
      expect_true(info.n_associations[0] == 1);
      expect_true(info.n_associations[1] == 1);
      expect_true(info.n_associations[2] == 1);
    }

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

    expect_true(params.vcov_start() == 29);
    expect_true(params.vcov_marker() == 29);
    expect_true(params.vcov_vary() == 38);
    expect_true(params.vcov_surv() == 74);
    expect_true(params.vcov_end() == 75);

    expect_true(params.n_params() == 75);
    expect_true(params.n_shared() == 6);
    expect_true(params.n_shared_surv() == 1);

    expect_true(params.va_par_start() == params.va_mean());
    expect_true(params.va_mean() == 75);
    expect_true(params.va_mean_end() == 82);
    expect_true(params.va_vcov() == 82);
    expect_true(params.va_vcov_end() == 131);
    expect_true(params.va_par_end() == params.va_vcov_end());

    expect_true(params.n_params_w_va() == 131);
    expect_true(params.n_va_params() == 56);

    expect_true(params.frailty_offset(0) == 0);
    expect_true(params.frailty_offset(1) == 0);

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

    expect_true(params.vcov_start<true>() == 29);
    expect_true(params.vcov_marker<true>() == 29);
    expect_true(params.vcov_vary<true>() == 35);
    expect_true(params.vcov_surv<true>() == 56);
    expect_true(params.vcov_end<true>() == 57);

    expect_true(params.n_params<true>() == 57);

    expect_true(params.va_par_start<true>() == params.va_mean<true>());
    expect_true(params.va_mean<true>() == 57);
    expect_true(params.va_mean_end<true>() == 64);
    expect_true(params.va_vcov<true>() == 64);
    expect_true(params.va_vcov_end<true>() == 92);
    expect_true(params.va_par_end<true>() == params.va_vcov_end<true>());

    expect_true(params.n_params_w_va<true>() == 92);
    expect_true(params.n_va_params<true>() == 35);

    std::vector<std::string> const param_names{params.param_names(true)};
    expect_true(param_names.size() == 57);

    std::vector<std::string> const true_names
      {
        // the fixed effects
        "mark1_fixef1", "mark1_fixef2", "mark1_fixef3",
        "mark2_fixef1",
        "mark3_fixef1", "mark3_fixef2", "mark3_fixef3",
        // time-varying fixed effects
        "mark1_fixef_vary1", "mark1_fixef_vary2",
        "mark2_fixef_vary1",
        "mark3_fixef_vary1", "mark3_fixef_vary2", "mark3_fixef_vary3",
        // parameters for the first survival outcome
        "surv1_fixef1", "surv1_fixef_vary1", "surv1_fixef_vary2",  "surv1_fixef_vary3",
        "surv1_fixef_vary4", "surv1_assoc1_1", "surv1_assoc2_1", "surv1_assoc3_1",
        // parameters for the second survival outcome
        "surv2_fixef1", "surv2_fixef2", "surv2_fixef3", "surv2_fixef_vary1",
        "surv2_fixef_vary2", "surv2_assoc1_1",  "surv2_assoc2_1",  "surv2_assoc3_1",
        // vcov for the markers' error term
        "vcov_marker1", "vcov_marker2", "vcov_marker3", "vcov_marker4", "vcov_marker5", "vcov_marker6",
        // vcov for the shared effect
        "vcov_vary1", "vcov_vary2", "vcov_vary3", "vcov_vary4", "vcov_vary5", "vcov_vary6", "vcov_vary7", "vcov_vary8", "vcov_vary9", "vcov_vary10", "vcov_vary11", "vcov_vary12", "vcov_vary13", "vcov_vary14", "vcov_vary15", "vcov_vary16", "vcov_vary17", "vcov_vary18", "vcov_vary19", "vcov_vary20", "vcov_vary21",
        // vcov for the frailties
        "vcov_surv1"
      };

    for(size_t i = 0; i < true_names.size(); ++i)
      expect_true(param_names[i] == true_names[i]);

    std::vector<std::string> const va_param_names{params.va_param_names(true)};
    expect_true(va_param_names.size() == 35);

    std::vector<std::string> const va_true_names
    {
      "mark1_rng1", "mark1_rng2",
      "mark2_rng1", "mark2_rng2", "mark2_rng3",
      "mark3_rng1",
      "frailty2",
      "VA_vcov1", "VA_vcov2", "VA_vcov3", "VA_vcov4", "VA_vcov5", "VA_vcov6", "VA_vcov7", "VA_vcov8", "VA_vcov9", "VA_vcov10", "VA_vcov11", "VA_vcov12", "VA_vcov13", "VA_vcov14", "VA_vcov15", "VA_vcov16", "VA_vcov17", "VA_vcov18", "VA_vcov19", "VA_vcov20", "VA_vcov21", "VA_vcov22", "VA_vcov23", "VA_vcov24", "VA_vcov25", "VA_vcov26", "VA_vcov27", "VA_vcov28"
    };

    for(size_t i = 0; i < va_true_names.size(); ++i)
      expect_true(va_param_names[i] == va_true_names[i]);
  }

  test_that("Works with both markers and survival outcomes and different number of association parameters") {
    subset_params params;
    params.add_marker({ 3L, 2L, 2L});
    params.add_marker({ 1L, 1L, 3L});
    params.add_marker({ 3L, 3L, 1L});

    params.add_surv({ 1L, 4L, {0, 0, 4}, true});
    params.add_surv({ 3L, 2L, {3, 2, 1}, true});

    expect_true(params.marker_info().size() == 3);
    expect_true(params.surv_info().size() == 2);
    expect_true(params.surv_info()[0].n_associations[0] == 0);
    expect_true(params.surv_info()[0].n_associations[1] == 0);
    expect_true(params.surv_info()[0].n_associations[2] == 4);
    expect_true(params.surv_info()[1].n_associations[0] == 3);
    expect_true(params.surv_info()[1].n_associations[1] == 2);
    expect_true(params.surv_info()[1].n_associations[2] == 1);

    expect_true(params.fixef_marker(0) == 0);
    expect_true(params.fixef_marker(1) == 3);
    expect_true(params.fixef_marker(2) == 4);

    expect_true(params.fixef_vary_marker(0) == 7);
    expect_true(params.fixef_vary_marker(1) == 9);
    expect_true(params.fixef_vary_marker(2) == 10);

    expect_true(params.fixef_surv(0) == 13);
    expect_true(params.fixef_vary_surv(0) == 14);
    expect_true(params.association(0) == 18);

    expect_true(params.fixef_surv(1) == 22);
    expect_true(params.fixef_vary_surv(1) == 25);
    expect_true(params.association(1) == 27);

    expect_true(params.vcov_start() == 33);
    expect_true(params.vcov_marker() == 33);
    expect_true(params.vcov_vary() == 42);
    expect_true(params.vcov_surv() == 78);
    expect_true(params.vcov_end() == 82);

    expect_true(params.n_params() == 82);
    expect_true(params.n_shared() == 6);
    expect_true(params.n_shared_surv() == 2);

    expect_true(params.va_par_start() == params.va_mean());
    expect_true(params.va_mean() == 82);
    expect_true(params.va_mean_end() == 90);
    expect_true(params.va_vcov() == 90);
    expect_true(params.va_vcov_end() == 154);
    expect_true(params.va_par_end() == params.va_vcov_end());

    expect_true(params.n_params_w_va() == 154);
    expect_true(params.n_va_params() == 72);

    expect_true(params.frailty_offset(0) == 0);
    expect_true(params.frailty_offset(1) == 1);

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

    expect_true(params.fixef_surv<true>(1) == 22);
    expect_true(params.fixef_vary_surv<true>(1) == 25);
    expect_true(params.association<true>(1) == 27);

    expect_true(params.vcov_start<true>() == 33);
    expect_true(params.vcov_marker<true>() == 33);
    expect_true(params.vcov_vary<true>() == 39);
    expect_true(params.vcov_surv<true>() == 60);
    expect_true(params.vcov_end<true>() == 63);

    expect_true(params.n_params<true>() == 63);

    expect_true(params.va_par_start<true>() == params.va_mean<true>());
    expect_true(params.va_mean<true>() == 63);
    expect_true(params.va_mean_end<true>() == 71);
    expect_true(params.va_vcov<true>() == 71);
    expect_true(params.va_vcov_end<true>() == 107);
    expect_true(params.va_par_end<true>() == params.va_vcov_end<true>());

    expect_true(params.n_params_w_va<true>() == 107);
    expect_true(params.n_va_params<true>() == 44);

    std::vector<std::string> const param_names{params.param_names(true)};
    expect_true(param_names.size() == 63);

    std::vector<std::string> const true_names
      {
        // the fixed effects
        "mark1_fixef1", "mark1_fixef2", "mark1_fixef3",
        "mark2_fixef1",
        "mark3_fixef1", "mark3_fixef2", "mark3_fixef3",
        // time-varying fixed effects
        "mark1_fixef_vary1", "mark1_fixef_vary2",
        "mark2_fixef_vary1",
        "mark3_fixef_vary1", "mark3_fixef_vary2", "mark3_fixef_vary3",
        // parameters for the first survival outcome
        "surv1_fixef1", "surv1_fixef_vary1", "surv1_fixef_vary2",  "surv1_fixef_vary3",
        "surv1_fixef_vary4", "surv1_assoc3_1", "surv1_assoc3_2", "surv1_assoc3_3", "surv1_assoc3_4",
        // parameters for the second survival outcome
        "surv2_fixef1", "surv2_fixef2", "surv2_fixef3", "surv2_fixef_vary1",
        "surv2_fixef_vary2", "surv2_assoc1_1", "surv2_assoc1_2", "surv2_assoc1_3",
        "surv2_assoc2_1",  "surv2_assoc2_2",  "surv2_assoc3_1",
        // vcov for the markers' error term
        "vcov_marker1", "vcov_marker2", "vcov_marker3", "vcov_marker4", "vcov_marker5", "vcov_marker6",
        // vcov for the shared effect
        "vcov_vary1", "vcov_vary2", "vcov_vary3", "vcov_vary4", "vcov_vary5", "vcov_vary6", "vcov_vary7", "vcov_vary8", "vcov_vary9", "vcov_vary10", "vcov_vary11", "vcov_vary12", "vcov_vary13", "vcov_vary14", "vcov_vary15", "vcov_vary16", "vcov_vary17", "vcov_vary18", "vcov_vary19", "vcov_vary20", "vcov_vary21",
        // vcov for the frailties
        "vcov_surv1", "vcov_surv2", "vcov_surv3",
      };

    for(size_t i = 0; i < true_names.size(); ++i)
      expect_true(param_names[i] == true_names[i]);

    std::vector<std::string> const va_param_names{params.va_param_names(true)};
    expect_true(va_param_names.size() == 44);

    std::vector<std::string> const va_true_names
    {
      "mark1_rng1", "mark1_rng2",
      "mark2_rng1", "mark2_rng2", "mark2_rng3",
      "mark3_rng1",
      "frailty1", "frailty2",
      "VA_vcov1", "VA_vcov2", "VA_vcov3", "VA_vcov4", "VA_vcov5", "VA_vcov6", "VA_vcov7", "VA_vcov8", "VA_vcov9", "VA_vcov10", "VA_vcov11", "VA_vcov12", "VA_vcov13", "VA_vcov14", "VA_vcov15", "VA_vcov16", "VA_vcov17", "VA_vcov18", "VA_vcov19", "VA_vcov20", "VA_vcov21", "VA_vcov22", "VA_vcov23", "VA_vcov24", "VA_vcov25", "VA_vcov26", "VA_vcov27", "VA_vcov28", "VA_vcov29", "VA_vcov30", "VA_vcov31", "VA_vcov32", "VA_vcov33", "VA_vcov34", "VA_vcov35", "VA_vcov36"
    };

    for(size_t i = 0; i < va_true_names.size(); ++i)
      expect_true(va_param_names[i] == va_true_names[i]);
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

    expect_true(params.vcov_start() == 13);
    expect_true(params.vcov_marker() == 13);
    expect_true(params.vcov_vary() == 22);
    expect_true(params.vcov_surv() == 58);
    expect_true(params.vcov_end() == 58);

    expect_true(params.n_params() == 58);
    expect_true(params.n_shared() == 6);
    expect_true(params.n_shared_surv() == 0);

    expect_true(params.va_par_start() == params.va_mean());
    expect_true(params.va_mean() == 58);
    expect_true(params.va_mean_end() == 64);
    expect_true(params.va_vcov() == 64);
    expect_true(params.va_vcov_end() == 100);
    expect_true(params.va_par_end() == params.va_vcov_end());

    expect_true(params.n_params_w_va() == 100);
    expect_true(params.n_va_params() == 42);

    // with the triangular matrices
    expect_true(params.fixef_marker<true>(0) == 0);
    expect_true(params.fixef_marker<true>(1) == 3);
    expect_true(params.fixef_marker<true>(2) == 4);

    expect_true(params.fixef_vary_marker<true>(0) == 7);
    expect_true(params.fixef_vary_marker<true>(1) == 9);
    expect_true(params.fixef_vary_marker<true>(2) == 10);

    expect_true(params.vcov_start<true>() == 13);
    expect_true(params.vcov_marker<true>() == 13);
    expect_true(params.vcov_vary<true>() == 19);
    expect_true(params.vcov_surv<true>() == 40);
    expect_true(params.vcov_end<true>() == 40);
    expect_true(params.n_params<true>() == 40);

    expect_true(params.va_par_start<true>() == params.va_mean<true>());
    expect_true(params.va_mean<true>() == 40);
    expect_true(params.va_mean_end<true>() == 46);
    expect_true(params.va_vcov<true>() == 46);
    expect_true(params.va_vcov_end<true>() == 67);
    expect_true(params.va_par_end<true>() == params.va_vcov_end<true>());

    expect_true(params.n_params_w_va<true>() == 67);
    expect_true(params.n_va_params<true>() == 27);

    std::vector<std::string> const param_names{params.param_names(true)};
    expect_true(param_names.size() == 40);

    std::vector<std::string> const true_names
    {
      // the fixed effects
      "mark1_fixef1", "mark1_fixef2", "mark1_fixef3",
      "mark2_fixef1",
      "mark3_fixef1", "mark3_fixef2", "mark3_fixef3",
      // time-varying fixed effects
      "mark1_fixef_vary1", "mark1_fixef_vary2",
      "mark2_fixef_vary1",
      "mark3_fixef_vary1", "mark3_fixef_vary2", "mark3_fixef_vary3",
      // vcov for the markers' error term
      "vcov_marker1", "vcov_marker2", "vcov_marker3", "vcov_marker4", "vcov_marker5", "vcov_marker6",
      // vcov for the shared effect
      "vcov_vary1", "vcov_vary2", "vcov_vary3", "vcov_vary4", "vcov_vary5", "vcov_vary6", "vcov_vary7", "vcov_vary8", "vcov_vary9", "vcov_vary10", "vcov_vary11", "vcov_vary12", "vcov_vary13", "vcov_vary14", "vcov_vary15", "vcov_vary16", "vcov_vary17", "vcov_vary18", "vcov_vary19", "vcov_vary20", "vcov_vary21",
    };

    for(size_t i = 0; i < true_names.size(); ++i)
      expect_true(param_names[i] == true_names[i]);

    std::vector<std::string> const va_param_names{params.va_param_names(true)};
    expect_true(va_param_names.size() == 27);

    std::vector<std::string> const va_true_names
    {
      "mark1_rng1", "mark1_rng2",
      "mark2_rng1", "mark2_rng2", "mark2_rng3",
      "mark3_rng1",
      "VA_vcov1", "VA_vcov2", "VA_vcov3", "VA_vcov4", "VA_vcov5", "VA_vcov6", "VA_vcov7", "VA_vcov8", "VA_vcov9", "VA_vcov10", "VA_vcov11", "VA_vcov12", "VA_vcov13", "VA_vcov14", "VA_vcov15", "VA_vcov16", "VA_vcov17", "VA_vcov18", "VA_vcov19", "VA_vcov20", "VA_vcov21",
    };

    for(size_t i = 0; i < va_true_names.size(); ++i)
      expect_true(va_param_names[i] == va_true_names[i]);
  }

  test_that("Works only with survival outcomes") {
    subset_params params;

    params.add_surv({ 1L, 4L, {}, true});
    params.add_surv({ 3L, 2L, {}, true});

    expect_true(params.marker_info().size() == 0);
    expect_true(params.surv_info().size() == 2);

    expect_true(params.fixef_surv(0) == 0);
    expect_true(params.fixef_vary_surv(0) == 1);
    expect_true(params.association(0) == 5);

    expect_true(params.fixef_surv(1) == 5);
    expect_true(params.fixef_vary_surv(1) == 8);
    expect_true(params.association(1) == 10);

    expect_true(params.vcov_start() == 10);
    expect_true(params.vcov_marker() == 10);
    expect_true(params.vcov_vary() == 10);
    expect_true(params.vcov_surv() == 10);
    expect_true(params.vcov_end() == 14);

    expect_true(params.n_params() == 14);
    expect_true(params.n_shared() == 0);
    expect_true(params.n_shared_surv() == 2);

    expect_true(params.va_par_start() == params.va_mean());
    expect_true(params.va_mean() == 14);
    expect_true(params.va_mean_end() == 16);
    expect_true(params.va_vcov() == 16);
    expect_true(params.va_vcov_end() == 20);
    expect_true(params.va_par_end() == params.va_vcov_end());

    expect_true(params.n_params_w_va() == 20);
    expect_true(params.n_va_params() == 6);

    expect_true(params.frailty_offset(0) == 0);
    expect_true(params.frailty_offset(1) == 1);

    // with the triangular matrices
    expect_true(params.fixef_surv<true>(0) == 0);
    expect_true(params.fixef_vary_surv<true>(0) == 1);
    expect_true(params.association<true>(0) == 5);

    expect_true(params.fixef_surv<true>(1) == 5);
    expect_true(params.fixef_vary_surv<true>(1) == 8);
    expect_true(params.association<true>(1) == 10);

    expect_true(params.vcov_start<true>() == 10);
    expect_true(params.vcov_marker<true>() == 10);
    expect_true(params.vcov_vary<true>() == 10);
    expect_true(params.vcov_surv<true>() == 10);
    expect_true(params.vcov_end<true>() == 13);

    expect_true(params.n_params<true>() == 13);

    expect_true(params.va_par_start<true>() == params.va_mean<true>());
    expect_true(params.va_mean<true>() == 13);
    expect_true(params.va_mean_end<true>() == 15);
    expect_true(params.va_vcov<true>() == 15);
    expect_true(params.va_vcov_end<true>() == 18);
    expect_true(params.va_par_end<true>() == params.va_vcov_end<true>());

    expect_true(params.n_params_w_va<true>() == 18);
    expect_true(params.n_va_params<true>() == 5);

    std::vector<std::string> const param_names{params.param_names(true)};
    expect_true(param_names.size() == 13);

    std::vector<std::string> const true_names
    {
      // parameters for the first survival outcome
      "surv1_fixef1", "surv1_fixef_vary1", "surv1_fixef_vary2",  "surv1_fixef_vary3",
      "surv1_fixef_vary4",
      // parameters for the second survival outcome
      "surv2_fixef1", "surv2_fixef2", "surv2_fixef3", "surv2_fixef_vary1",
      "surv2_fixef_vary2",
      // vcov for the frailties
      "vcov_surv1", "vcov_surv2", "vcov_surv3"
    };

    for(size_t i = 0; i < true_names.size(); ++i)
      expect_true(param_names[i] == true_names[i]);

    std::vector<std::string> const va_param_names{params.va_param_names(true)};
    expect_true(va_param_names.size() == 5);

    std::vector<std::string> const va_true_names
    {
      "frailty1", "frailty2",
      "VA_vcov1", "VA_vcov2", "VA_vcov3",
    };

    for(size_t i = 0; i < va_true_names.size(); ++i)
      expect_true(va_param_names[i] == va_true_names[i]);
  }

  test_that("Works only with survival outcomes without both frailties") {
    subset_params params;

    params.add_surv({ 1L, 4L, {}, false});
    params.add_surv({ 3L, 2L, {}, false});

    expect_true(params.marker_info().size() == 0);
    expect_true(params.surv_info().size() == 2);

    expect_true(params.fixef_surv(0) == 0);
    expect_true(params.fixef_vary_surv(0) == 1);
    expect_true(params.association(0) == 5);

    expect_true(params.fixef_surv(1) == 5);
    expect_true(params.fixef_vary_surv(1) == 8);
    expect_true(params.association(1) == 10);

    expect_true(params.vcov_start() == 10);
    expect_true(params.vcov_marker() == 10);
    expect_true(params.vcov_vary() == 10);
    expect_true(params.vcov_surv() == 10);
    expect_true(params.vcov_end() == 10);

    expect_true(params.n_params() == 10);
    expect_true(params.n_shared() == 0);
    expect_true(params.n_shared_surv() == 0);

    expect_true(params.va_par_start() == params.va_mean());
    expect_true(params.va_mean() == 10);
    expect_true(params.va_mean_end() == 10);
    expect_true(params.va_vcov() == 10);
    expect_true(params.va_vcov_end() == 10);
    expect_true(params.va_par_end() == params.va_vcov_end());

    expect_true(params.n_params_w_va() == 10);
    expect_true(params.n_va_params() == 0);

    expect_true(params.frailty_offset(0) == 0);
    expect_true(params.frailty_offset(1) == 0);

    // with the triangular matrices
    expect_true(params.fixef_surv<true>(0) == 0);
    expect_true(params.fixef_vary_surv<true>(0) == 1);
    expect_true(params.association<true>(0) == 5);

    expect_true(params.fixef_surv<true>(1) == 5);
    expect_true(params.fixef_vary_surv<true>(1) == 8);
    expect_true(params.association<true>(1) == 10);

    expect_true(params.vcov_start<true>() == 10);
    expect_true(params.vcov_marker<true>() == 10);
    expect_true(params.vcov_vary<true>() == 10);
    expect_true(params.vcov_surv<true>() == 10);
    expect_true(params.vcov_end<true>() == 10);

    expect_true(params.n_params<true>() == 10);

    expect_true(params.va_par_start<true>() == params.va_mean<true>());
    expect_true(params.va_mean<true>() == 10);
    expect_true(params.va_mean_end<true>() == 10);
    expect_true(params.va_vcov<true>() == 10);
    expect_true(params.va_vcov_end<true>() == 10);
    expect_true(params.va_par_end<true>() == params.va_vcov_end<true>());

    expect_true(params.n_params_w_va<true>() == 10);
    expect_true(params.n_va_params<true>() == 0);

    std::vector<std::string> const param_names{params.param_names(true)};
    expect_true(param_names.size() == 10);

    std::vector<std::string> const true_names
    {
      // parameters for the first survival outcome
      "surv1_fixef1", "surv1_fixef_vary1", "surv1_fixef_vary2",  "surv1_fixef_vary3",
      "surv1_fixef_vary4",
      // parameters for the second survival outcome
      "surv2_fixef1", "surv2_fixef2", "surv2_fixef3", "surv2_fixef_vary1",
      "surv2_fixef_vary2"
    };

    for(size_t i = 0; i < true_names.size(); ++i)
      expect_true(param_names[i] == true_names[i]);

    std::vector<std::string> const va_param_names{params.va_param_names(true)};
    expect_true(va_param_names.size() == 0);
  }
}
