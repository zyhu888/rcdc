## Parallel RobustCDC runs
##
## RobustCDC() evaluates candidate values of K serially. For large analyses,
## run robustCDC_MCMC() separately for each candidate K, save one result per
## task, and select the final clustering by the smallest WAIC.

.libPaths("~/rlibs")

library(RCDC)
library(parallel)

hyper <- list(
  sigma = 10,
  a_alpha = 1,
  b_alpha = 1,
  a_beta = 1,
  b_beta = 1,
  mu0 = c(0, 0, 0),
  iota0 = c(1, 1, 1),
  a_tau = c(1, 1, 1),
  b_tau = c(0.01, 0.01, 0.01),
  F_cat = c(8),
  gamma0 = c(1),
  d = 1
)

candidate_K <- 5:30
n_iter <- 2000
burn_in <- 1000
n_cores <- length(candidate_K)

task_index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
id_run <- seq_along(network_list)
id <- id_run[task_index]

out_dir <- "/path/to/results"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(id)

run_one_K <- function(Ki) {
  tryCatch({
    out_file <- file.path(out_dir, sprintf("id%dK%d.rds", id, Ki))

    if (file.exists(out_file)) {
      message(sprintf("id %d K %d already exists, skipping: %s", id, Ki, out_file))
      return(NULL)
    }

    fit <- robustCDC_MCMC(
      A_tilde = network_list[[id]]$network_array,
      X_cont = as.matrix(network_list[[id]]$cov[, 2:4]),
      X_cat = NULL,
      K = Ki,
      n_iter = n_iter,
      burn_in = burn_in,
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

    z_post <- fit$z[(burn_in + 1):n_iter, ]
    z_hat <- apply(
      z_post,
      2,
      function(x) as.integer(names(which.max(table(x))))
    )

    tmp_file <- paste0(out_file, ".tmp")
    saveRDS(list(z = z_hat, waic = waic_val), file = tmp_file)
    file.rename(tmp_file, out_file)

    invisible(NULL)
  }, error = function(e) {
    message(sprintf("ERROR: id %d K %d failed: %s", id, Ki, e$message))
    invisible(NULL)
  })
}

mclapply(
  candidate_K,
  run_one_K,
  mc.cores = n_cores,
  mc.set.seed = TRUE
)
