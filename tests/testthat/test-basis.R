test_that("The C++ version of poly gives the right result", {
  in_x <- 2:5
  obj_truth <- poly(in_x, degree = 3)
  out_x <- 1:6
  truth <- predict(obj_truth, out_x)

  # without an intercept
  obj_cpp <- poly_term(in_x, degree = 3)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(truth), ignore_attr = TRUE)

  # with an intercept
  obj_cpp <- poly_term(in_x, degree = 3, intercept = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, rbind(1, t(truth)), ignore_attr = TRUE)

  # without an intercept and raw is TRUE
  obj_cpp <- poly_term(in_x, degree = 3, raw = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(outer(out_x, 1:3, `^`)))

  # with an intercept and raw is TRUE
  obj_cpp <- poly_term(in_x, degree = 3, raw = TRUE, intercept = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(outer(out_x, 0:3, `^`)))

  # without an intercept and degree == 0
  obj_cpp <- poly_term(in_x, degree = 0)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, matrix(0, 0, length(out_x)))

  # with an intercept and degree == 0
  obj_cpp <- poly_term(in_x, degree = 0, intercept = TRUE)
  expect_s3_class(obj_cpp, "poly_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, matrix(1, 1, length(out_x)))
})

test_that("The C++ version of bs gives the right result", {
  # without an intercept
  in_x <- 2:5
  obj_truth <- bs(in_x, df = 4)
  out_x <- c(2, 2.5, 3, 3.5, 4, 4.5, 5)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- bs_term(in_x, df = 4)
  expect_s3_class(obj_cpp, "bs_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)

  # with an intercept
  obj_truth <- bs(in_x, df = 4, intercept = TRUE)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- bs_term(in_x, df = 4, intercept = TRUE)
  expect_s3_class(obj_cpp, "bs_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)
})

test_that("The C++ version of ns gives the right result", {
  # without an intercept
  in_x <- 2:5
  obj_truth <- ns(in_x, df = 4)
  out_x <- c(2, 2.5, 3, 3.5, 4, 4.5, 5)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- ns_term(in_x, df = 4)
  expect_s3_class(obj_cpp, "ns_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)

  # with an intercept
  obj_truth <- ns(in_x, df = 4, intercept = TRUE)
  truth <- predict(obj_truth, out_x)

  obj_cpp <- ns_term(in_x, df = 4, intercept = TRUE)
  expect_s3_class(obj_cpp, "ns_term")
  expect_equal(obj_cpp$time, in_x)

  expansion <- VAJointSurv:::eval_expansion(obj_cpp, out_x)
  expect_equal(expansion, t(truth),  ignore_attr = TRUE)
})
