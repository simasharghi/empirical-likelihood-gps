# analysis/simulate_data.R
# Simulate datasets for EL-GPS continuous treatment experiments
# Saves one .RData per replication containing daData = list(mydata.A, mydata.B, mydata.C)

# -----------------------
# USER SETTINGS
# -----------------------
out_dir <- "Saved_datasets_n1000"  # e.g., "Saved_datasets_n500" or "Saved_datasets_n1000"
sim_number <- 1000                # number of replications
N <- 1000                         # sample size per replication
seed <- 13555

# epsilon type: "normal", "mixture", or "bimodal"
epsilon_type <- "bimodal"

# For mixture epsilon
sigma1 <- 1
sigma2 <- 10
mix_p  <- 0.95

# For bimodal epsilon
bimodal_p <- 0.30    # P(component = -1.5)
bimodal_mu1 <- -1.5
bimodal_mu2 <-  1.5
bimodal_sd  <-  1

# True treatment effect in Y model
beta_T <- 0.4

# -----------------------
# Libraries (keep minimal)
# -----------------------
# (No heavy libraries needed just to simulate)
# -----------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

draw_epsilon <- function(N, type = c("normal", "mixture", "bimodal"),
                         sigma1 = 1, sigma2 = 10, mix_p = 0.95,
                         bimodal_p = 0.30, mu1 = -1.5, mu2 = 1.5, sd = 1) {
  type <- match.arg(type)
  if (type == "normal") {
    eps <- rnorm(N, mean = 0, sd = 1)
  } else if (type == "mixture") {
    eps <- ifelse(runif(N) < mix_p,
                  rnorm(N, mean = 0, sd = sigma1),
                  rnorm(N, mean = 0, sd = sigma2))
    eps <- eps - mean(eps)  # center mixture
  } else if (type == "bimodal") {
    eps <- ifelse(rbinom(N, 1, bimodal_p) == 1,
                  rnorm(N, mean = mu1, sd = sd),
                  rnorm(N, mean = mu2, sd = sd))
    eps <- eps - mean(eps)  # center bimodal for stability
  }
  eps
}

simulate_one <- function(N, epsilon_type, beta_T,
                         sigma1, sigma2, mix_p,
                         bimodal_p, bimodal_mu1, bimodal_mu2, bimodal_sd) {
  
  # Covariates
  X1  <- rnorm(N, 0, 1)
  X2  <- rnorm(N, 0, 1)
  X3  <- rbinom(N, 1, 0.50)
  X4  <- rbinom(N, 1, 0.30)
  X5  <- rnorm(N, 0, 1)
  X6  <- rnorm(N, 0, 1)
  X7  <- rnorm(N, 0, 1)
  X8  <- rbinom(N, 1, 0.50)
  X9  <- rbinom(N, 1, 0.50)
  X10 <- rbinom(N, 1, 0.50)
  
  x <- data.frame(X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
  
  # Epsilon
  eps <- draw_epsilon(
    N, type = epsilon_type,
    sigma1 = sigma1, sigma2 = sigma2, mix_p = mix_p,
    bimodal_p = bimodal_p, mu1 = bimodal_mu1, mu2 = bimodal_mu2, sd = bimodal_sd
  )
  
  # m(x) scenarios
  mx.A <- 6 + 0.30*X1 + 0.65*X2 - 0.35*X3 - 0.40*X4
  
  mx.B <- 6 + 0.30*(X1 > 0.5) + 0.65*(X2 < 0) - 0.35*(X3 == 1) - 0.40*(X4 == 0) +
    0.65*(X1 > 0) * (X1 > 1)
  
  mx.C <- 6 + 0.30*(X1 > 0.5) + 0.65*(X2 < 0) - 0.35*(X3 == 1) - 0.40*(X4 == 0) +
    0.65*(X1 > 0) * (X1 > 1) +
    0.30*(X1 > 0) * (X4 == 1) -
    0.65*(X2 > 0.30) * (X3 == 0)
  
  # Treatments
  T.A <- mx.A + eps
  T.B <- mx.B + eps
  T.C <- mx.C + eps
  
  # Outcomes (NOTE: fix exponent: X1^2 not X1^{2})
  Y.A <- 3.85 + beta_T*T.A + 0.30*X1 + 0.36*X2 + 0.73*X3 - 0.20*X4 + 0.25*(X1^2) + rnorm(N, 0, 1)
  Y.B <- 3.85 + beta_T*T.B + 0.30*X1 + 0.36*X2 + 0.73*X3 - 0.20*X4 + 0.25*(X1^2) + rnorm(N, 0, 1)
  Y.C <- 3.85 + beta_T*T.C + 0.30*X1 + 0.36*X2 + 0.73*X3 - 0.20*X4 + 0.25*(X1^2) + rnorm(N, 0, 1)
  
  mydata.A <- data.frame(T = T.A, x, Y = Y.A)
  mydata.B <- data.frame(T = T.B, x, Y = Y.B)
  mydata.C <- data.frame(T = T.C, x, Y = Y.C)
  
  list(mydata.A = mydata.A, mydata.B = mydata.B, mydata.C = mydata.C)
}

# -----------------------
# Run simulation
# -----------------------
set.seed(seed)
ensure_dir(out_dir)

for (j in seq_len(sim_number)) {
  daData <- simulate_one(
    N = N,
    epsilon_type = epsilon_type,
    beta_T = beta_T,
    sigma1 = sigma1, sigma2 = sigma2, mix_p = mix_p,
    bimodal_p = bimodal_p, bimodal_mu1 = bimodal_mu1, bimodal_mu2 = bimodal_mu2, bimodal_sd = bimodal_sd
  )
  
  save(daData, file = file.path(out_dir, sprintf("ZhuData_sim%04d.RData", j)))
  if (j %% 50 == 0) message("Saved replication: ", j)
}

message("Done. Files written to: ", out_dir)
