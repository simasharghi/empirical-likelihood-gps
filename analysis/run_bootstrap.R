# analysis/run_bootstrap.R
# Clean bootstrap driver for EL-GPS estimators
# - restartable (skips sims already processed)
# - robust to estimator failures
# - percentile bootstrap CI for treatment effect
# - saves one .rds per simulation file (recommended)

# ---------------------------
# User inputs (edit these)
# ---------------------------

# Folder that contains simulation .RData files (each file loads an object 'daData')
sim_data_dir <- "Saved_datasets"   # <-- change to your path
# Folder where results will be written
out_dir      <- "Saved_Analysis"   # <-- change to your path

# Number of bootstrap replicates
B <- 200
# Seed for reproducibility (set to NULL to skip)
seed <- 123

# Treatment effect coefficient name in the outcome model:
# In your estimators you fit svyglm(Y ~ Treat), so coefficient is typically "Treat"
treat_coef_name <- "Treat"

# Estimators to run (functions must be loaded in the session)
estimator_names <- c(
  "my.fun.EL.No.Boost.1",
  "my.fun.EL.No.Boost.2",
  "my.fun.Norm.No.Boost.1",
  "my.fun.Norm.No.Boost.2",
  "my.fun.EL.Boost",
  "my.fun.Normal.Boost"
)

# ---------------------------
# Load functions
# ---------------------------
# These should exist in your repo:
source("EL_Core/el_core.R")
source("EL_Core/estimators.R")

# ---------------------------
# Helpers
# ---------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_fit <- function(estimator_name, dat) {
  # Returns list(ok, fit, err)
  res <- tryCatch(
    {
      f <- get(estimator_name, mode = "function")
      fit <- f(mydata = dat)
      list(ok = TRUE, fit = fit, err = NULL)
    },
    error = function(e) list(ok = FALSE, fit = NULL, err = conditionMessage(e))
  )
  res
}

extract_treat <- function(fit, coef_name = "Treat") {
  # Works for svyglm objects; returns numeric scalar or NA
  if (is.null(fit)) return(NA_real_)
  cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
  if (is.null(cf)) return(NA_real_)
  if (!(coef_name %in% names(cf))) return(NA_real_)
  as.numeric(cf[[coef_name]])
}

extract_wald_ci <- function(fit, coef_name = "Treat") {
  # Uses confint() but indexes by name, not numeric position
  # Returns c(lo, hi) or c(NA, NA)
  if (is.null(fit)) return(c(NA_real_, NA_real_))
  ci <- tryCatch(stats::confint(fit), error = function(e) NULL)
  if (is.null(ci)) return(c(NA_real_, NA_real_))
  if (is.vector(ci) && length(ci) == 2 && !is.null(names(ci))) {
    # rare edge case; treat as not supported
    return(c(NA_real_, NA_real_))
  }
  if (!(coef_name %in% rownames(ci))) return(c(NA_real_, NA_real_))
  c(as.numeric(ci[coef_name, 1]), as.numeric(ci[coef_name, 2]))
}

bootstrap_theta <- function(estimator_name, dat, B = 200, seed = NULL,
                            coef_name = "Treat") {
  n <- nrow(dat)
  if (!is.null(seed)) set.seed(seed)
  
  # Point estimate
  base <- safe_fit(estimator_name, dat)
  theta_hat <- extract_treat(base$fit, coef_name)
  wald_ci   <- extract_wald_ci(base$fit, coef_name)
  
  # Bootstrap draws
  thetas <- rep(NA_real_, B)
  ok_vec <- rep(FALSE, B)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    fb <- safe_fit(estimator_name, dat[idx, , drop = FALSE])
    thetas[b] <- extract_treat(fb$fit, coef_name)
    ok_vec[b] <- is.finite(thetas[b])
  }
  
  thetas_ok <- thetas[is.finite(thetas)]
  n_ok <- length(thetas_ok)
  
  # Bootstrap SE (standard)
  se_boot_sd <- if (n_ok >= 2) stats::sd(thetas_ok) else NA_real_
  
  # Bootstrap RMSE around theta_hat (optional diagnostic)
  se_boot_rmse <- if (n_ok >= 1 && is.finite(theta_hat)) {
    sqrt(mean((thetas_ok - theta_hat)^2))
  } else NA_real_
  
  # Percentile CI (recommended default)
  ci_boot <- if (n_ok >= 10) {
    as.numeric(stats::quantile(thetas_ok, c(0.025, 0.975), names = FALSE))
  } else c(NA_real_, NA_real_)
  
  list(
    estimator = estimator_name,
    theta_hat = theta_hat,
    ci_wald_lo = wald_ci[1],
    ci_wald_hi = wald_ci[2],
    B = B,
    n_boot_ok = n_ok,
    se_boot_sd = se_boot_sd,
    se_boot_rmse = se_boot_rmse,
    ci_boot_lo = ci_boot[1],
    ci_boot_hi = ci_boot[2],
    boot_thetas = thetas,     # keep for internal use; you can drop for GitHub outputs
    base_ok = base$ok,
    base_err = base$err
  )
}

run_one_simfile <- function(sim_file, out_dir, B, seed, estimator_names,
                            coef_name = "Treat") {
  out_file <- file.path(out_dir, paste0(basename(sim_file), ".rds"))
  if (file.exists(out_file)) return(invisible(NULL))  # restartable
  
  env <- new.env(parent = emptyenv())
  load(sim_file, envir = env)
  
  if (!exists("daData", envir = env)) {
    warning("File does not contain 'daData': ", sim_file)
    return(invisible(NULL))
  }
  
  daData <- get("daData", envir = env)
  
  if (is.null(daData$mydata.A)) {
    warning("daData$mydata.A is missing: ", sim_file)
    return(invisible(NULL))
  }
  
  dat <- daData$mydata.A
  
  # Run all estimators
  res_list <- lapply(estimator_names, function(m) {
    bootstrap_theta(m, dat, B = B, seed = seed, coef_name = coef_name)
  })
  
  # Tidy summary table
  summary_df <- do.call(rbind, lapply(res_list, function(r) {
    data.frame(
      sim_file = basename(sim_file),
      estimator = r$estimator,
      theta_hat = r$theta_hat,
      ci_wald_lo = r$ci_wald_lo,
      ci_wald_hi = r$ci_wald_hi,
      B = r$B,
      n_boot_ok = r$n_boot_ok,
      se_boot_sd = r$se_boot_sd,
      se_boot_rmse = r$se_boot_rmse,
      ci_boot_lo = r$ci_boot_lo,
      ci_boot_hi = r$ci_boot_hi,
      base_ok = r$base_ok,
      base_err = ifelse(is.null(r$base_err), NA_character_, r$base_err),
      stringsAsFactors = FALSE
    )
  }))
  
  # Save: summary + (optional) full bootstrap draws
  out_obj <- list(
    summary = summary_df,
    details = res_list
  )
  saveRDS(out_obj, out_file)
  
  invisible(out_file)
}

# ---------------------------
# Main
# ---------------------------

ensure_dir(out_dir)

sim_files <- list.files(sim_data_dir, full.names = TRUE)
sim_files <- sample(sim_files, replace = FALSE)

cat("Found", length(sim_files), "simulation files.\n")
cat("Writing results to:", out_dir, "\n")

for (k in seq_along(sim_files)) {
  f <- sim_files[[k]]
  cat(sprintf("[%d/%d] %s\n", k, length(sim_files), basename(f)))
  run_one_simfile(
    sim_file = f,
    out_dir = out_dir,
    B = B,
    seed = if (!is.null(seed)) seed + k else NULL,
    estimator_names = estimator_names,
    coef_name = treat_coef_name
  )
}

cat("Done.\n")
