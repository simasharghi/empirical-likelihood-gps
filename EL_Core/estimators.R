# ============================================================
# estimators.R
# Estimator functions (EL-kernel and Normal) 
# ============================================================

# Required packages for these estimators:
# - mgcv (for gam)
# - survey (for svyglm)

# ---------
# Utilities
# ---------

.check_required_cols <- function(mydata, treat_col, outcome_col) {
  if (!is.data.frame(mydata)) stop("mydata must be a data.frame")
  if (!treat_col %in% names(mydata)) stop("Treatment column not found: ", treat_col)
  if (!outcome_col %in% names(mydata)) stop("Outcome column not found: ", outcome_col)
}

.make_treatment_df <- function(mydata, outcome_col) {
  # Treatment model should not include outcome
  mydata[, setdiff(names(mydata), outcome_col), drop = FALSE]
}

.safe_ci <- function(fit, coef_name = "Treat") {
  ci <- tryCatch(confint(fit), error = function(e) NULL)
  if (is.null(ci)) return(c(NA_real_, NA_real_))
  if (!is.matrix(ci)) ci <- as.matrix(ci)
  if (!coef_name %in% rownames(ci)) return(c(NA_real_, NA_real_))
  c(ci[coef_name, 1], ci[coef_name, 2])
}

.stabilize_ratio <- function(num, den, eps = 1e-12, trim_quantile = NULL) {
  num <- pmax(as.numeric(num), eps)
  den <- pmax(as.numeric(den), eps)
  w <- num / den
  
  # remove pathological values
  w[!is.finite(w)] <- NA_real_
  
  if (!is.null(trim_quantile)) {
    cap <- as.numeric(stats::quantile(w, probs = trim_quantile, na.rm = TRUE))
    w <- pmin(w, cap)
  }
  w
}

# Normal density at observed values (standardized)
.normal_density_standardized <- function(x, center, scale) {
  z <- (x - center) / scale
  stats::dnorm(z, mean = 0, sd = 1)
}

# -----------------------
# Estimator implementations
# -----------------------

# 2: EL No Boost, Linear 1 (All covariates)
my.fun.EL.No.Boost.1 <- function(mydata,
                                 form = "T~.",
                                 treat_col = "T",
                                 outcome_col = "Y",
                                 weight_floor = 1e-12,
                                 trim_quantile = NULL,
                                 kde_bw = "nrd0") {
  
  .check_required_cols(mydata, treat_col, outcome_col)
  
  Tvec <- mydata[[treat_col]]
  Yvec <- mydata[[outcome_col]]
  
  # Numerator: EL-weighted KDE of T
  num_obj <- el_weighted_kde_at_points(Tvec, constraints = c("mean", "var"), bw = kde_bw)
  ps_num <- num_obj$density_at_points
  
  # Denominator model: linear regression for m(X)
  treat_df <- .make_treatment_df(mydata, outcome_col)
  fit_mx <- stats::lm(stats::as.formula(form), data = treat_df)
  
  mu_hat <- stats::predict(fit_mx, newdata = treat_df)
  resid <- Tvec - mu_hat
  
  # Denominator: EL-weighted KDE of residuals
  den_obj <- el_weighted_kde_at_points(resid, constraints = c("mean", "var"), bw = kde_bw)
  ps_den <- den_obj$density_at_points
  
  w <- .stabilize_ratio(ps_num, ps_den, eps = weight_floor, trim_quantile = trim_quantile)
  
  dataset <- data.frame(Y = Yvec, Treat = Tvec, weight_model = w)
  design <- survey::svydesign(ids = ~1, weights = ~weight_model, data = dataset)
  
  fit <- survey::svyglm(Y ~ Treat, family = gaussian(), design = design)
  
  attr(fit, "diagnostics") <- list(
    max_weight = max(w, na.rm = TRUE),
    min_weight = min(w, na.rm = TRUE),
    n_na_weight = sum(is.na(w)),
    min_ps_den = min(ps_den, na.rm = TRUE),
    min_ps_num = min(ps_num, na.rm = TRUE),
    el_num_converged = num_obj$el$converged,
    el_den_converged = den_obj$el$converged
  )
  
  fit
}

# 3: EL No Boost, GAM 2 (Variable selection via select=TRUE)
my.fun.EL.No.Boost.2 <- function(mydata,
                                 form = "T~s(X1)+s(X2)+X3+X4+s(X5)+s(X6)+s(X7)+X8+X9+X10",
                                 treat_col = "T",
                                 outcome_col = "Y",
                                 weight_floor = 1e-12,
                                 trim_quantile = NULL,
                                 kde_bw = "nrd0") {
  
  .check_required_cols(mydata, treat_col, outcome_col)
  
  Tvec <- mydata[[treat_col]]
  Yvec <- mydata[[outcome_col]]
  
  # Numerator: EL-weighted KDE of T
  num_obj <- el_weighted_kde_at_points(Tvec, constraints = c("mean", "var"), bw = kde_bw)
  ps_num <- num_obj$density_at_points
  
  treat_df <- .make_treatment_df(mydata, outcome_col)
  fit_mx <- mgcv::gam(stats::as.formula(form), data = treat_df, select = TRUE)
  
  mu_hat <- stats::predict(fit_mx, newdata = treat_df, type = "response")
  resid <- Tvec - mu_hat
  
  # Denominator: EL-weighted KDE of residuals
  den_obj <- el_weighted_kde_at_points(resid, constraints = c("mean", "var"), bw = kde_bw)
  ps_den <- den_obj$density_at_points
  
  w <- .stabilize_ratio(ps_num, ps_den, eps = weight_floor, trim_quantile = trim_quantile)
  
  dataset <- data.frame(Y = Yvec, Treat = Tvec, weight_model = w)
  design <- survey::svydesign(ids = ~1, weights = ~weight_model, data = dataset)
  
  fit <- survey::svyglm(Y ~ Treat, family = gaussian(), design = design)
  
  attr(fit, "diagnostics") <- list(
    max_weight = max(w, na.rm = TRUE),
    min_weight = min(w, na.rm = TRUE),
    n_na_weight = sum(is.na(w)),
    min_ps_den = min(ps_den, na.rm = TRUE),
    min_ps_num = min(ps_num, na.rm = TRUE),
    el_num_converged = num_obj$el$converged,
    el_den_converged = den_obj$el$converged
  )
  
  fit
}

# 4: Normal No Boost, Linear 1 (All covariates)
my.fun.Norm.No.Boost.1 <- function(mydata,
                                   form = "T~.",
                                   treat_col = "T",
                                   outcome_col = "Y",
                                   weight_floor = 1e-12,
                                   trim_quantile = NULL) {
  
  .check_required_cols(mydata, treat_col, outcome_col)
  
  Tvec <- mydata[[treat_col]]
  Yvec <- mydata[[outcome_col]]
  
  # Numerator density for T: standardize by sample mean/sd (robust + clear)
  muT <- mean(Tvec, na.rm = TRUE)
  sdT <- sd(Tvec, na.rm = TRUE)
  ps_num <- .normal_density_standardized(Tvec, center = muT, scale = sdT)
  
  # Denominator model
  treat_df <- .make_treatment_df(mydata, outcome_col)
  fit_mx <- stats::lm(stats::as.formula(form), data = treat_df)
  
  mu_hat <- stats::predict(fit_mx, newdata = treat_df)
  resid <- Tvec - mu_hat
  
  # Residual scale: df-corrected linear-model sigma (defendable)
  N <- sum(is.finite(resid))
  k <- length(stats::coef(fit_mx))
  s <- sqrt(sum(resid^2, na.rm = TRUE) / max(1, (N - k)))
  
  ps_den <- stats::dnorm(resid / s, mean = 0, sd = 1)
  
  w <- .stabilize_ratio(ps_num, ps_den, eps = weight_floor, trim_quantile = trim_quantile)
  
  dataset <- data.frame(Y = Yvec, Treat = Tvec, weight_model = w)
  design <- survey::svydesign(ids = ~1, weights = ~weight_model, data = dataset)
  fit <- survey::svyglm(Y ~ Treat, family = gaussian(), design = design)
  
  attr(fit, "diagnostics") <- list(
    max_weight = max(w, na.rm = TRUE),
    min_weight = min(w, na.rm = TRUE),
    n_na_weight = sum(is.na(w)),
    min_ps_den = min(ps_den, na.rm = TRUE),
    min_ps_num = min(ps_num, na.rm = TRUE)
  )
  
  fit
}

# 5: Normal No Boost, GAM 2 (Variable selection)
my.fun.Norm.No.Boost.2 <- function(mydata,
                                   form = "T~s(X1)+s(X2)+X3+X4+s(X5)+s(X6)+s(X7)+X8+X9+X10",
                                   treat_col = "T",
                                   outcome_col = "Y",
                                   weight_floor = 1e-12,
                                   trim_quantile = NULL) {
  
  .check_required_cols(mydata, treat_col, outcome_col)
  
  Tvec <- mydata[[treat_col]]
  Yvec <- mydata[[outcome_col]]
  
  muT <- mean(Tvec, na.rm = TRUE)
  sdT <- sd(Tvec, na.rm = TRUE)
  ps_num <- .normal_density_standardized(Tvec, center = muT, scale = sdT)
  
  treat_df <- .make_treatment_df(mydata, outcome_col)
  fit_mx <- mgcv::gam(stats::as.formula(form), data = treat_df, select = TRUE)
  
  mu_hat <- stats::predict(fit_mx, newdata = treat_df, type = "response")
  resid <- Tvec - mu_hat
  
  # GAM residual scale: use sample sd (defendable and avoids EDF debates in a public repo)
  s <- sd(resid, na.rm = TRUE)
  ps_den <- stats::dnorm(resid / s, mean = 0, sd = 1)
  
  w <- .stabilize_ratio(ps_num, ps_den, eps = weight_floor, trim_quantile = trim_quantile)
  
  dataset <- data.frame(Y = Yvec, Treat = Tvec, weight_model = w)
  design <- survey::svydesign(ids = ~1, weights = ~weight_model, data = dataset)
  fit <- survey::svyglm(Y ~ Treat, family = gaussian(), design = design)
  
  attr(fit, "diagnostics") <- list(
    max_weight = max(w, na.rm = TRUE),
    min_weight = min(w, na.rm = TRUE),
    n_na_weight = sum(is.na(w)),
    min_ps_den = min(ps_den, na.rm = TRUE),
    min_ps_num = min(ps_num, na.rm = TRUE)
  )
  
  fit
}
