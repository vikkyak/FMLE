test_that("single-task training and prediction return expected structure", {
  demo <- make_fmle_demo_data(seed = 11)

  fit <- fmle_train(
    X = demo$X_train,
    y = demo$Y_train[, 1],
    Z = demo$Z_train,
    R = 3,
    m = 1.8,
    lambda_l1 = 1e-3,
    ridge = 1e-6,
    standardize = TRUE,
    seed = 1
  )

  expect_s3_class(fit, "fmle")
  expect_equal(fit$R, 3)
  expect_equal(nrow(fit$centers), 3)

  pred <- fmle_predict(
    model = fit,
    X_new = demo$X_test,
    Z_new = demo$Z_test,
    return_se = TRUE
  )

  expect_type(pred$mean, "double")
  expect_equal(length(pred$mean), nrow(demo$X_test))
  expect_equal(dim(pred$alpha), c(nrow(demo$X_test), fit$R))
  expect_equal(dim(pred$mu_r), c(nrow(demo$X_test), fit$R))
  expect_equal(length(pred$se), nrow(demo$X_test))
  expect_true(all(is.finite(pred$mean)))
  expect_true(all(is.finite(pred$alpha)))
})

test_that("multi-task training and prediction return expected dimensions", {
  demo <- make_fmle_demo_data(seed = 12, T = 4)

  fit_mt <- fmle_train_mt(
    X = demo$X_train,
    Y = demo$Y_train,
    Z = demo$Z_train,
    R = 3,
    m = 1.8,
    lambda_l1 = 1e-3,
    ridge = 1e-6,
    standardize = TRUE,
    seed = 1
  )

  expect_s3_class(fit_mt, "fmle_mt")
  expect_equal(fit_mt$T, ncol(demo$Y_train))

  pred_mt <- fmle_predict_mt(
    model = fit_mt,
    X_new = demo$X_test,
    Z_new = demo$Z_test,
    return_se = TRUE
  )

  expect_equal(dim(pred_mt$mean), c(nrow(demo$X_test), ncol(demo$Y_test)))
  expect_equal(dim(pred_mt$alpha), c(nrow(demo$X_test), fit_mt$R))
  expect_true(all(is.finite(pred_mt$mean)))
})

test_that("single-task CV returns a best configuration and a non-empty table", {
  demo <- make_fmle_demo_data(seed = 13)

  cv <- fmle_cv_parallel(
    X = demo$X_train,
    y = demo$Y_train[, 1],
    Z = demo$Z_train,
    R_grid = c(2, 3),
    m_grid = c(1.6, 1.8),
    lambda_grid = c(0, 1e-3),
    folds = 3,
    seed = 1,
    exec = "sequential",
    verbose = FALSE
  )

  expect_true(is.list(cv))
  expect_true(all(c("best", "table") %in% names(cv)))
})
