## Parallel RobustCDC runs
##
## This example mirrors the small runnable example in ?RobustCDC, but evaluates
## candidate K values in parallel by calling robustCDC_MCMC() directly.

## Optional: uncomment if you use a custom R library on a cluster.
## .libPaths("~/rlibs")

library(RCDC)
library(parallel)

seed <- 123
set.seed(seed)

n <- 60
K0 <- 3
M <- 4
burn <- 200
niter <- 500

rho <- c(0.10, 0.30)
zeta <- c(1.0, 1.0)

hyper <- list(
  sigma = 10,
  a_alpha = 1,
  b_alpha = 1,
  a_beta = 1,
  b_beta = 1,
  mu0 = rep(0, K0),
  iota0 = rep(1, K0),
  a_tau = rep(1, K0),
  b_tau = rep(0.01, K0),
  F_cat = c(K0, K0),
  gamma0 = c(1, 1),
  d = 1
)

B <- matrix(-1.5, K0, K0)
diag(B) <- 0

mu <- diag(K0) - matrix(1 / K0, K0, K0)
mu <- mu * (10 / sqrt(2))

tau2 <- diag(3^2, K0)
theta <- rnorm(n, mean = zeta[1], sd = zeta[2])
alpha <- matrix(rho[1], K0, K0)
beta <- matrix(rho[2], K0, K0)
cat_rate <- c(0.7, 0.7)

sim <- sim_data_new_thetainput(
  seed = seed,
  n = n,
  K = K0,
  M = M,
  p_cont = K0,
  p_cat = 2,
  F_s = K0,
  theta_true = theta,
  B = B,
  alpha = alpha,
  beta = beta,
  mu = mu,
  tau2 = tau2,
  correct_rate = cat_rate
)

candidate_K <- 2:5
n_cores <- length(candidate_K)

out_dir <- file.path(tempdir(), "RCDC_parallel_example")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

run_one_K <- function(Ki) {
  tryCatch({
    out_file <- file.path(out_dir, sprintf("K%d.rds", Ki))

    if (file.exists(out_file)) {
      message(sprintf("K %d already exists, skipping: %s", Ki, out_file))
      return(NULL)
    }

    fit <- robustCDC_MCMC(
      A_tilde = sim$A_tilde,
      X_cont = sim$X_cont,
      X_cat = sim$X_cat,
      K = Ki,
      n_iter = niter,
      burn_in = burn,
      hyper = hyper
    )

    log_lik_total <- cbind(
      fit$loglik_net,
      if (!is.null(fit$loglik_cont)) fit$loglik_cont,
      if (!is.null(fit$loglik_cat)) fit$loglik_cat
    )

    waic_val <- if (anyNA(log_lik_total)) {
      NA_real_
    } else {
      suppressWarnings(loo::waic(log_lik_total))$estimates["waic", "Estimate"]
    }

    z_post <- fit$z[(burn + 1):niter, ]
    z_hat <- apply(
      z_post,
      2,
      function(x) as.integer(names(which.max(table(x))))
    )

    tmp_file <- paste0(out_file, ".tmp")
    saveRDS(list(K = Ki, z = z_hat, waic = waic_val), file = tmp_file)
    file.rename(tmp_file, out_file)

    invisible(NULL)
  }, error = function(e) {
    message(sprintf("ERROR: K %d failed: %s", Ki, e$message))
    invisible(NULL)
  })
}

mclapply(
  candidate_K,
  run_one_K,
  mc.cores = n_cores,
  mc.set.seed = TRUE
)

res <- lapply(candidate_K, function(Ki) {
  readRDS(file.path(out_dir, sprintf("K%d.rds", Ki)))
})

waic_table <- data.frame(
  K = vapply(res, `[[`, integer(1), "K"),
  WAIC = vapply(res, `[[`, numeric(1), "waic")
)
waic_table$Delta_WAIC <- waic_table$WAIC - min(waic_table$WAIC, na.rm = TRUE)

best <- res[[which.min(waic_table$WAIC)]]
z_hat <- best$z
z_true <- sim$truth$z

print(waic_table)
cat("ARI:", aricode::ARI(z_hat, z_true), "\n")
cat("NMI:", aricode::NMI(z_hat, z_true), "\n")
cat("Selected K:", best$K, "\n")
