
# ============================================================
# el_core.R
# Core empirical likelihood (EL) utilities for moment constraints
# Used to compute EL weights and EL-weighted kernel density estimates
#
# This file is designed to be GitHub-safe and data-agnostic.
# ============================================================

#' EL ratio statistic (-2 log R)
#' g: n x d matrix of constraints evaluated at each observation
#' lambda: d-vector of Lagrange multipliers
AELRS <- function(g, lambda) {
  2 * sum(log(1 + g %*% lambda))
}

#' Equation for lambda: sum_i g_i / (1 + g_i^T lambda) = 0
EquationForLambda <- function(g, lambda) {
  apply(g / as.vector(1 + g %*% lambda), 2, sum)
}

#' Jacobian of EquationForLambda
D_EquationForLambda <- function(g, lambda) {
  -t(g) %*% (g / as.vector(1 + g %*% lambda)^2)
}

#' Newton-type solver for lambda with backtracking to enforce feasibility
#' Returns list(lambda=..., converged=TRUE/FALSE, error=..., n_iter=...)
SeekingForLambda <- function(g, max.iter = 200, tol = 1e-8, trace_ = FALSE) {
  g <- as.matrix(g)
  p <- ncol(g)
  
  lambda <- rep(0, p)
  gamma <- 1
  k <- 0
  out.iter <- 0
  
  converged <- FALSE
  last_err <- NA_real_
  
  repeat {
    resid <- EquationForLambda(g, lambda)
    last_err <- max(abs(resid))
    
    if (trace_) {
      cat("#", out.iter, "max|resid|=", last_err, "lambda=", paste(round(lambda, 6), collapse = ","), "\n")
    }
    
    if (last_err < tol) {
      converged <- TRUE
      break
    }
    
    if (out.iter >= max.iter) {
      converged <- FALSE
      break
    }
    
    # Newton step; guard against singular Jacobian
    J <- D_EquationForLambda(g, lambda)
    Delta <- tryCatch(
      -drop(solve(J, resid, tol = 1e-40)),
      error = function(e) rep(NA_real_, p)
    )
    
    if (any(!is.finite(Delta))) {
      converged <- FALSE
      break
    }
    
    delta <- gamma * Delta
    inner.iter <- 0
    
    # backtracking line search to keep (1 + g %*% (lambda + delta)) > 0
    repeat {
      if (inner.iter >= max.iter) break
      
      prob <- drop(1 + g %*% (lambda + delta))
      if (any(prob <= 0) || !is.finite(AELRS(g, lambda + delta)) || AELRS(g, lambda + delta) < AELRS(g, lambda)) {
        delta <- delta / 2
        inner.iter <- inner.iter + 1
      } else {
        break
      }
    }
    
    lambda <- lambda + delta
    k <- k + 1
    gamma <- 1 / sqrt(k)
    out.iter <- out.iter + 1
  }
  
  list(
    lambda = drop(lambda),
    converged = converged,
    max_abs_resid = last_err,
    n_iter = out.iter
  )
}

#' Compute EL (or adjusted EL) weights under mean constraints.
#' x: n x d matrix; mu: scalar or d-vector target (typically 0)
#'
#' Returns:
#'  - weights: EL weights (sum to 1)
#'  - lambda, converged, diagnostics
AELforMean <- function(x,
                       mu = 0,
                       an = NULL,
                       adjust = FALSE,
                       max.iter = 200,
                       tol = 1e-8,
                       trace_ = FALSE) {
  
  x <- as.matrix(x)
  n <- nrow(x)
  d <- ncol(x)
  
  # handle mu scalar or vector
  if (length(mu) == 1) mu <- rep(mu, d)
  if (length(mu) != d) stop("Length of mu must be 1 or equal to ncol(x).")
  
  # build constraint matrix g
  if (adjust) {
    if (is.null(an)) an <- log(n) / 2
    xbar <- colMeans(x, na.rm = TRUE)
    
    g <- matrix(NA_real_, n + 1, d)
    g[1:n, ] <- t(t(x) - mu)
    g[n + 1, ] <- -an * (xbar - mu)
  } else {
    g <- t(t(x) - mu)
  }
  
  sol <- SeekingForLambda(g, max.iter = max.iter, tol = tol, trace_ = trace_)
  lambda <- sol$lambda
  
  # If not converged, still attempt weights if feasible; otherwise return NA weights
  denom <- as.vector(1 + g %*% lambda)
  feasible <- all(is.finite(denom)) && all(denom > 0)
  
  if (!feasible) {
    w <- rep(NA_real_, nrow(g))
    wsum <- NA_real_
  } else {
    w <- (1 / denom) / nrow(g)
    wsum <- sum(w)
    # normalize explicitly to avoid numeric drift
    w <- w / wsum
  }
  
  Wstar <- if (feasible) AELRS(g, lambda) else NA_real_
  err_lambda <- if (feasible) EquationForLambda(g, lambda) else rep(NA_real_, d)
  
  list(
    AELRS = Wstar,
    method = if (!adjust) "EL" else paste0("AEL (an=", an, ")"),
    weights = w,
    lambda = lambda,
    converged = sol$converged && feasible,
    max_abs_resid = sol$max_abs_resid,
    n_iter = sol$n_iter,
    error_in_equation = err_lambda
  )
}

#' EL-weighted kernel density evaluation at observed points
#' z: numeric vector
#' constraints: character vector of constraints, currently supports c("mean","var")
#' returns list(density_at_z, weights, kde_object, el_object)
el_weighted_kde_at_points <- function(z,
                                      constraints = c("mean", "var"),
                                      bw = "nrd0") {
  z <- as.numeric(z)
  z <- z[is.finite(z)]
  if (length(z) < 5) stop("Not enough finite values in z to estimate density.")
  
  m <- mean(z, na.rm = TRUE)
  v <- var(z, na.rm = TRUE)
  
  # Build constraint matrix for EL
  G <- NULL
  if ("mean" %in% constraints) {
    G <- cbind(G, z - m)
  }
  if ("var" %in% constraints) {
    G <- cbind(G, (z - m)^2 - v)
  }
  if (is.null(G)) stop("No constraints specified.")
  
  el <- AELforMean(x = G, mu = 0, adjust = FALSE)
  w <- el$weights
  
  if (any(!is.finite(w)) || any(w < 0)) stop("EL weights not feasible (non-finite or negative).")
  
  den <- density(z, weights = w, bw = bw, na.rm = TRUE)
  den_at <- approx(den$x, den$y, xout = z, rule = 2)$y
  
  list(
    density_at_points = as.numeric(den_at),
    weights = w,
    kde = den,
    el = el,
    z_used = z
  )
}
