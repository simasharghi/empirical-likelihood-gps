# analysis/make_plots_highlight_case1.R
# Highlight Case 1 plots (m(x)=A, epsilon=Normal)
# Reads saved bootstrap outputs (per-simulation .RData files that contain `temp` data frame).

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ---- User settings (EDIT THESE) ----
# Point to where your saved analysis outputs live (NOT committed to GitHub).
# Recommended: set this to a local path on your machine.
base_dir <- "/Volumes/projects/TIDES-II/AirPollution/Causality/simulations_new/Sim_subset/HighlightSim1"
n_values <- c(500, 1000)

# True effect used in simulations (adjust if needed)
theta_true <- 0.4

# Where to save figures (relative to repo)
fig_dir <- file.path("analysis", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: map estimator names to display labels ----
label_estimator_mx <- function(estimator) {
  out <- rep(NA_character_, length(estimator))
  out[estimator == "my.fun.EL.Boost"]         <- "GBM"
  out[estimator == "my.fun.EL.No.Boost.1"]    <- "Linear.1"
  out[estimator == "my.fun.EL.No.Boost.2"]    <- "Linear.2"
  out[estimator == "my.fun.Norm.No.Boost.1"]  <- "Linear.1"
  out[estimator == "my.fun.Norm.No.Boost.2"]  <- "Linear.2"
  out[estimator == "my.fun.Normal.Boost"]     <- "GBM"
  out
}

label_modeled_epsilon <- function(estimator) {
  out <- rep(NA_character_, length(estimator))
  out[estimator %in% c("my.fun.EL.Boost","my.fun.EL.No.Boost.1","my.fun.EL.No.Boost.2")] <- "W.EL.Kernel"
  out[estimator %in% c("my.fun.Norm.No.Boost.1","my.fun.Norm.No.Boost.2","my.fun.Normal.Boost")] <- "Normal"
  out
}

# ---- Helper: load and stack per-simulation results into one data.frame ----
read_saved_analysis <- function(dir_path) {
  files <- list.files(dir_path, full.names = TRUE)
  files <- files[grepl("\\.RData$|\\.rda$|\\.RDA$", files, ignore.case = TRUE)]
  
  if (length(files) == 0) stop("No .RData files found in: ", dir_path)
  
  res <- lapply(files, function(f) {
    e <- new.env(parent = emptyenv())
    load(f, envir = e)
    
    # Expect `temp` to exist (your saved object)
    if (!exists("temp", envir = e)) return(NULL)
    e$temp
  })
  
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0) stop("Files loaded but no `temp` object found in: ", dir_path)
  
  do.call(rbind, res)
}

# ---- Summarize metrics for plotting ----
summarize_case <- function(df, n_label, theta_true) {
  # df contains stacked `temp` rows across simulations (and across estimators)
  # Expected columns (based on your code):
  # estimator, fit.A.Coef, C.I.Low.A, C.I.UP.A, C.I.Low.A.Boot, C.I.UP.A.Boot
  
  split_by_est <- split(df, df$estimator)
  
  out <- do.call(rbind, lapply(split_by_est, function(x) {
    data.frame(
      m.x = "A",
      epsilon = "Normal",
      n = n_label,
      estimator = unique(x$estimator),
      num.sims = nrow(x),
      BIAS = mean(x$fit.A.Coef - theta_true),
      SD = sd(x$fit.A.Coef),
      RMSE = sqrt(mean((x$fit.A.Coef - theta_true)^2)),
      MSE = mean((x$fit.A.Coef - theta_true)^2),
      SEE.Boot = mean((x$C.I.UP.A.Boot - x$fit.A.Coef)/1.96, na.rm = TRUE),
      SEE.LS   = mean((x$C.I.UP.A      - x$fit.A.Coef)/1.96, na.rm = TRUE),
      Coverage.Prob.Boot = mean(x$C.I.Low.A.Boot <= theta_true & x$C.I.UP.A.Boot >= theta_true, na.rm = TRUE),
      Avg.Length.Boot    = mean(x$C.I.UP.A.Boot - x$C.I.Low.A.Boot, na.rm = TRUE),
      Coverage.Prob.LS   = mean(x$C.I.Low.A <= theta_true & x$C.I.UP.A >= theta_true, na.rm = TRUE),
      Avg.Length.LS      = mean(x$C.I.UP.A - x$C.I.Low.A, na.rm = TRUE)
    )
  }))
  
  rownames(out) <- NULL
  
  out$Estimator.of.m.x <- label_estimator_mx(out$estimator)
  out$Modeled.Epsilon  <- label_modeled_epsilon(out$estimator)
  
  out
}

# ---- Plot theme (simple, readable) ----
my_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
  axis.text.x = element_text(face = "bold", size = 12, colour = "black"),
  axis.text.y = element_text(face = "bold", size = 12, colour = "black"),
  axis.title  = element_text(face = "bold", size = 13, colour = "black"),
  strip.text  = element_text(size = 12, colour = "black")
)

# ---- Main: build Case 1 object ----
Case.1 <- do.call(rbind, lapply(n_values, function(n) {
  dir_path <- file.path(base_dir, paste0("Saved_Analysis_n", n))
  df <- read_saved_analysis(dir_path)
  summarize_case(df, n_label = paste0("n=", n), theta_true = theta_true)
}))

Case.1$n <- factor(Case.1$n, levels = c("n=500","n=1000"))

# ---- Plots (Bias, RMSE, Coverage, Avg Length) ----
# NOTE: Using shape mapping. We avoid hard-coding shape values to keep it robust.
p_bias <- ggplot(Case.1, aes(x = Modeled.Epsilon, y = BIAS)) +
  geom_point(aes(shape = Estimator.of.m.x), size = 4) +
  facet_grid(. ~ n) +
  labs(title = "Bias", x = "", y = "BIAS") +
  my_theme

p_rmse <- ggplot(Case.1, aes(x = Modeled.Epsilon, y = RMSE)) +
  geom_point(aes(shape = Estimator.of.m.x), size = 4) +
  facet_grid(. ~ n) +
  labs(title = "RMSE", x = "", y = "RMSE") +
  my_theme

p_cov <- ggplot(Case.1, aes(x = Modeled.Epsilon, y = Coverage.Prob.Boot)) +
  geom_point(aes(shape = Estimator.of.m.x), size = 4) +
  facet_grid(. ~ n) +
  labs(title = "Bootstrap CI Coverage", x = "", y = "Coverage Probability") +
  my_theme

p_len <- ggplot(Case.1, aes(x = Modeled.Epsilon, y = Avg.Length.Boot)) +
  geom_point(aes(shape = Estimator.of.m.x), size = 4) +
  facet_grid(. ~ n) +
  labs(title = "Bootstrap CI Average Length", x = "", y = "Average Length") +
  my_theme

combinedCase1 <- (p_bias + p_rmse + p_cov + p_len) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(combinedCase1)

# ---- Save ----
png_file <- file.path(fig_dir, "highlight_case1_summary.png")
ggsave(filename = png_file, plot = combinedCase1, width = 11, height = 7, dpi = 300)

pdf_file <- file.path(fig_dir, "highlight_case1_summary.pdf")
ggsave(filename = pdf_file, plot = combinedCase1, width = 11, height = 7)

message("Saved: ", png_file)
message("Saved: ", pdf_file)
