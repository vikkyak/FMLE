# FMLE

**FMLE** implements a fuzzy mixture of linear experts for predicting protein abundance from single-cell transcriptomic features, with shared fuzzy gates and regime-specific linear experts.

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
remotes::install_github("YOUR_GITHUB_USERNAME/FMLE")
