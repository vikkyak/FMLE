# FMLE

**FMLE** is an R package for regime-aware prediction of protein abundance from single-cell transcriptomic features using a fuzzy mixture of linear experts. It combines fuzzy c-means gating in a low-dimensional latent space with expert-specific linear predictors on high-dimensional gene expression features, enabling both interpretable and accurate protein prediction across heterogeneous cellular regimes.

The package supports:

- single-task prediction (`fmle_train()`, `fmle_predict()`)
- multi-task prediction (`fmle_train_mt()`, `fmle_predict_mt()`)
- cross-validation over the number of experts, fuzzifier, and L1 penalty (`fmle_cv_parallel()`, `fmle_cv_mt_parallel()`)
- fuzzy c-means gating (`fcm_fit()`)
- predictive uncertainty decomposition from the fitted experts

## Installation

```r
# install.packages("remotes")
remotes::install_local("FMLE")
# or
remotes::install_github("vikkyak/FMLE")
```



## Quickstart

The package ships with a tiny demo object for examples and the vignette.

## Single-task example
```r
library(FMLE)

demo <- readRDS(system.file("extdata", "fmle_demo.rds", package = "FMLE"))

X_train <- demo$X_train
X_test  <- demo$X_test
Y_train <- demo$Y_train
Y_test  <- demo$Y_test
Z_train <- demo$Z_train
Z_test  <- demo$Z_test

cv <- fmle_cv_parallel(
  X = X_train,
  y = Y_train[, 1],
  Z = Z_train,
  R_grid = c(2, 3),
  m_grid = c(1.6, 1.8),
  lambda_grid = c(0, 1e-3),
  folds = 3,
  seed = 1,
  exec = "sequential",
  verbose = FALSE
)

best <- cv$best

fit <- fmle_train(
  X = X_train,
  y = Y_train[, 1],
  Z = Z_train,
  R = best$R,
  m = best$m,
  lambda_l1 = best$lambda,
  ridge = 1e-6,
  standardize = TRUE,
  seed = 1
)

pred <- fmle_predict(
  model = fit,
  X_new = X_test,
  Z_new = Z_test,
  return_se = TRUE
)

cor(pred$mean, Y_test[, 1])
mean((pred$mean - Y_test[, 1])^2)

```

## Multi-task example

```r
fit_mt <- fmle_train_mt(
  X = X_train,
  Y = Y_train,
  Z = Z_train,
  R = 2,
  m = 1.8,
  lambda_l1 = 1e-3,
  ridge = 1e-6,
  standardize = TRUE,
  seed = 1
)

pred_mt <- fmle_predict_mt(
  model = fit_mt,
  X_new = X_test,
  Z_new = Z_test,
  return_se = TRUE
)

dim(pred_mt$mean)
colMeans((pred_mt$mean - Y_test)^2)
```

## Important preprocessing note

For **single-task FMLE**, `fmle_cv_parallel()` internally applies cap/log/scale preprocessing to the response before fold-wise fitting and evaluation. In contrast, `fmle_train()` fits the response exactly as supplied.

After selecting `(R, m, lambda)` by cross-validation, refit the full model using the response scale you intend to use for the final model and evaluation.