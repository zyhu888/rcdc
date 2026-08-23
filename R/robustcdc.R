## ---- MCMC function ----
robustCDC_MCMC <- function(
    A_tilde,     
    X_cont = NULL,
    X_cat  = NULL,
    K,
    n_iter = 3000,
    burn_in = 500,
    hyper,
    init = NULL # list of hyperparameters
) {
  ## ------------------------------------------------
  ## 0. validate dimensions
  ## ------------------------------------------------
  if (length(dim(A_tilde)) != 3) stop("A_tilde must be a 3D array n   n   M")
  
  n <- dim(A_tilde)[1]
  M <- dim(A_tilde)[3]  
  
  if (dim(A_tilde)[2] != n) stop("A_tilde must be n   n   M")
  
  if (!is.null(X_cont)) {
    X_cont <- as.matrix(X_cont)
    storage.mode(X_cont) <- "numeric"
    if (nrow(X_cont) != n) {
      stop("X_cont must have the same number of rows as nodes in A_tilde")
    }
  }
  
  if (!is.null(X_cat)) {
    X_cat <- as.matrix(X_cat)
    if (nrow(X_cat) != n) {
      stop("X_cat must have the same number of rows as nodes in A_tilde")
    }
  }
  n_cont <- if (is.null(X_cont)) 0 else ncol(X_cont)
  n_cat  <- if (is.null(X_cat))  0 else ncol(X_cat)
  
  ## ------------------------------------------------
  ## 1. hyperparameters
  ## ------------------------------------------------
  sigma       <- hyper$sigma
  
  a_alpha <- hyper$a_alpha
  b_alpha <- hyper$b_alpha  
  a_beta  <- hyper$a_beta
  b_beta  <- hyper$b_beta   
  
  # for continuous covariates
  mu0     <- hyper$mu0      
  iota0   <- hyper$iota0    
  a_tau   <- hyper$a_tau    
  b_tau   <- hyper$b_tau    
  
  # for categorical covariates
  F_cat   <- hyper$F_cat    
  gamma0  <- hyper$gamma0   
  
  # for mixing proportions  
  d   <- hyper$d    
  
  
  ## ------------------------------------------------
  ## 2. Pre-compute constants (outside MCMC loop)
  ## ------------------------------------------------
  # Observation statistics (constant across iterations)
  Count_1 <- rowSums(A_tilde, dims = 2)  # sum_m I(A_tilde^(m)_ij = 1)
  Count_0 <- M - Count_1                  # sum_m I(A_tilde^(m)_ij = 0)
  
  # Upper triangle indices 
  upper_idx <- which(upper.tri(matrix(0, n, n), diag = FALSE)) # column order
  # include (i<j): column-major order
  edge_indices <- arrayInd(upper_idx, .dim = c(n, n)) # column order
  # include (k=k) and (k<l): column-major order
  comm_indices <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  
  # Design matrix Y for   and B sampling
  # Each row y_ij corresponds to edge (i,j) with i < j
  # Dimension: n(n-1)/2   (K(K-1)/2 + n)
  # Note: delta = (b,  ), B only has off-diagonal entries (B_kk = 0)
  n_edges <- length(upper_idx)
  n_b_params <- K * (K - 1) / 2  # Only different community pairs
  n_params <- n_b_params + n # length of all   + B
  Y <- matrix(0, n_edges, n_params)
  
  # Fill   part (columns n_b_params+1 : n_params)
  Y[cbind(1:n_edges, n_b_params + edge_indices[, 1])] <- 1
  Y[cbind(1:n_edges, n_b_params + edge_indices[, 2])] <- 1
  #The B-part (columns 1 : n_b_params) depends on z, will fill inside loop
  
  ## ------------------------------------------------
  ## 3. Initialization parameters
  ## ------------------------------------------------
  
  if (is.null(init)) {
    ## ----- Initialize A and z (le2018estimating) ----
    A       <- (Count_1 >= (M / 2)) * 1L
    if (any(rowSums(A)==0)) {iso <- which(rowSums(A)==0); n <- nrow(A); j <- ifelse(iso==n,1,iso+1); A[cbind(iso,j)] <- 1; A[cbind(j,iso)] <- 1}
    z <- as.numeric(
      suppressMessages(
        C4(A, K = K)
      )$cluster
    )
    
    z       <- canonical_relabel(z,K)$z
    ## ----- Initialize theta and B (Le-style parameter mapping) -----
    Psi_matrix     <- matrix(0, K, K)
    edge_comm_idx <- cbind(pmin(z[edge_indices[,1]], z[edge_indices[,2]]), 
                           pmax(z[edge_indices[,1]], z[edge_indices[,2]]))
    
    # Estimate Psi_kl based on empirical edge density in each block
    for (p in 1:nrow(comm_indices)) {
      k <- comm_indices[p, 1]
      l <- comm_indices[p, 2]
      
      # Logical mask for the current block (k, l)
      block_mask <- (edge_comm_idx[,1] == k & edge_comm_idx[,2] == l)
      
      # Calculate empirical edge probability (w_hat) from majority vote matrix A 
      w_hat <- mean(A[upper_idx][block_mask])
      if (is.nan(w_hat) || !any(block_mask)) {
        val <- 0
      } else {
        w_hat <- pmin(pmax(w_hat, 1e-6), 1 - 1e-6)
        val <- qlogis(w_hat)
      }
      # Map to DCSBM parameter Psi_kl using logit link
      Psi_matrix[k, l] <- Psi_matrix[l, k] <- val
    }
    
    # 1) community-level theta from diagonal (uses B_kk = 0 constraint)
    theta_block <- 0.5 * diag(Psi_matrix)     # length K, indexed by community label 1..K
    
    # 2) node-level theta aligned to node order (1..n) using z
    theta <- theta_block[z]                   # length n, indexed by node i
    
    # 3) community-level B reconstructed from psi = theta_k + theta_l + B_kl
    B <- Psi_matrix - outer(theta_block, theta_block, "+")
    
    ## ----- Initialize alpha and beta (le2018estimating) -----
    alpha <- matrix(0, K, K)
    beta  <- matrix(0, K, K)
    for (p in 1:nrow(comm_indices)) {
      k <- comm_indices[p, 1]
      l <- comm_indices[p, 2]
      
      block_mask <- (edge_comm_idx[,1] == k & edge_comm_idx[,2] == l)
      # Extract relevant statistics for the current block
      S_b <- Count_1[upper_idx][block_mask] # Observed counts of 1's
      A_b <- A[upper_idx][block_mask]       # Current latent adjacency A
      
      # Update alpha (False Positive rate): P(observed=1 | A=0) 
      # Sum of observed 1's where A=0 divided by total possible observations
      idx_A0 <- (A_b == 0)
      if (sum(idx_A0) > 0) {
        alpha[k,l] <- alpha[l,k] <- sum(S_b[idx_A0]) / (M * sum(idx_A0))
      } else {
        # Fallback to prior mean if the block contains no A=1 edges (e.g., 2 communities have no connection)
        alpha[k,l] <- alpha[l,k] <- a_alpha / (a_alpha + b_alpha)
      }
      
      # Update beta (False Negative rate): P(observed=0 | A=1)
      # Sum of observed 0's where A=1 divided by total possible observations
      idx_A1 <- (A_b == 1)
      if (sum(idx_A1) > 0) {
        beta[k,l]  <- beta[l,k]  <- sum(M - S_b[idx_A1]) / (M * sum(idx_A1))
      } else {
        # Fallback to prior mean
        beta[k,l]  <- beta[l,k]  <- a_beta / (a_beta + b_beta)
      }
    }
    
    ## ----- Initialize covariate -----
    if (n_cont > 0) {
      
      #  _{k } = average of X within each community
      mu <- rowsum(X_cont, z) / tabulate(z, nbins = K)
      
      #   _{k } = average squared deviation within each community
      tau2 <- rowsum( (X_cont - mu[z, ])^2 , z) / tabulate(z, nbins = K)
      
    } else {
      mu   <- NULL
      tau2 <- NULL
    }
    
    if (n_cat > 0) {
      
      xi <- vector("list", n_cat)
      
      for (s in 1:n_cat) {
        
        # count_{k,f} = number of nodes in community k with category f
        count_matrix <- matrix(
          tabulate(interaction(z, X_cat[, s]), nbins = K * F_cat[s]),
          nrow = K, ncol = F_cat[s], byrow = FALSE
        )
        
        #  _{ks } = normalized frequencies in each community
        xi[[s]] <- count_matrix / rowSums(count_matrix)
      }
      
    } else {
      xi <- NULL
    }
    
  } else {
    
    ## ----- user-supplied initialization -----
    A      <- init$A
    z      <- init$z
    theta  <- init$theta
    B      <- init$B
    alpha  <- init$alpha
    beta   <- init$beta
    
    mu     <- if (n_cont > 0) init$mu   else NULL
    tau2   <- if (n_cont > 0) init$tau2 else NULL
    xi     <- if (n_cat  > 0) init$xi   else NULL
  }
  
  
  ## ------------------------------------------------
  ## 4. storage container
  ## ------------------------------------------------
  n_save <- n_iter
  samples_z     <- matrix(NA, n_save, n)
  
  keep_idx <- (burn_in + 1L):n_iter
  n_keep <- length(keep_idx)
  
  samples_loglik_net <- matrix(NA_real_, nrow = n_keep, ncol = n_edges)
  
  
  samples_loglik_cont <- if (n_cont > 0) matrix(0, nrow = n_keep, ncol = n) else NULL
  samples_loglik_cat  <- if (n_cat  > 0) matrix(0, nrow = n_keep, ncol = n) else NULL
  
  eps <- 1e-12
  clip01 <- function(p) pmin(pmax(p, eps), 1 - eps)
  
  ## ------------------------------------------------
  ## 5. MCMC 
  ## ------------------------------------------------
  for (iter in 1:n_iter) {
    
    ## ===== 1. Update A =====
    # 1.1 Construct Linear Predictor Psi (Prior Log-Odds)
    #  _ij =  _i +  _j + B_{z_i, z_j}
    B_mat <- B[z, z]                    # Map B matrix to node pairs
    Psi <- outer(theta, theta, FUN = `+`) + B_mat # Do not need indicator as long as B has 0 diagnal
    diag(Psi) <- 0                      # No self-loops
    
    # 1.2 Error rates 
    # All map to node pairs
    Alpha_mat <- alpha[z, z]
    Beta_mat  <- beta[z, z]
    
    # 1.3 Calculate Log-Likelihoods
    # p*_{ij,1} = exp( )   (1- )^n1    ^n0
    p_star_1 <- exp(Psi) * (1 - Beta_mat)^Count_1 * Beta_mat^Count_0
    # p*_{ij,0} =  ^n1   (1- )^n0
    p_star_0 <- Alpha_mat^Count_1 * (1 - Alpha_mat)^Count_0
    
    # 1.4 Compute Posterior Probability P(A_ij=1 | Data)
    Prob_A1 <- p_star_1 / (p_star_1 + p_star_0)
    
    # 1.5 Sample A Efficiently (Upper Triangle Only for Undirected Graph)
    # Initialize new adjacency matrix
    A_new <- matrix(0, n, n)
    # Sample Bernoulli for upper triangle (vectorized)
    A_new[upper_idx] <- rbinom(n_edges, 1, Prob_A1[upper_idx])
    # Symmetrize: copy upper triangle to lower triangle
    A_new <- A_new + t(A_new)
    # Update A
    A <- A_new
    
    ## ===== 2. Update W (PG) =====
    # Sample Polya-Gamma augmentation variables W_{ij} ~ PG(1,  _{ij})
    # Construct symmetric W matrix
    W_new <- matrix(0, n, n)
    W_new[upper_idx] <- pgdraw::pgdraw(rep(1, n_edges), Psi[upper_idx])
    W_new <- W_new + t(W_new)
    # Update W
    W <- W_new
    
    ## ===== 3. Update   and B (joint sampling) =====
    # 3.1 Construct kappa vector:  _ij = A_ij - 1/2
    kappa <- (A - 0.5)[upper_idx]
    
    # 3.2 Update Y matrix B-part (depends on current z)
    # Reset B columns to 0
    Y[, 1:n_b_params] <- 0
    
    # Only edges between DIFFERENT communities involve B parameters
    # First get all node pair community assignment(k,l) and let k <=l
    k_vec <- pmin(z[edge_indices[, 1]], z[edge_indices[, 2]])
    l_vec <- pmax(z[edge_indices[, 1]], z[edge_indices[, 2]])
    # Only fill for k != l (inter-community edges)
    is_inter <- k_vec != l_vec
    
    if (any(is_inter)) {
      # Map (k, l) to linear index in upper triangle of B (k < l, excluding diagonal)
      # left part of Y follows standard R "upper.tri" order (Column-Major): 
      # (1,2), (1,3), (2,3), (1,4)...
      # Formula: (l-1)*(l-2)/2 + k
      
      b_idx_vec <- (l_vec[is_inter] - 1) * (l_vec[is_inter] - 2) / 
        2 + k_vec[is_inter]
      
      # Set 1s in Y
      Y[cbind(which(is_inter), b_idx_vec)] <- 1
    }
    # By zeyu: Check all here, correct, follows formula in manuscript:  =Y 
    # 3.3 Construct diagonal weight matrix  
    Omega <- W[upper_idx]
    
    # 3.4 Posterior covariance and mean for   = (b,  )
    # V = [Y'   Y + (1/  )I]^{-1}
    # V^-1 = Y'   Y + (1/  )I
    V_inv <- crossprod(Y * sqrt(Omega)) + (1/sigma^2) * diag(n_params)
    
    # Cholesky decomposition for stable inversion(better than solve)
    # Decompose the symmetric positive-definite matrix V^-1 into 
    # V_inv_chol' * V_inv_chol
    V_inv_chol <- chol(V_inv)
    # Use V_inv_chol to compute V, stable and efficient
    V <- chol2inv(V_inv_chol)
    
    # Posterior mean   = V * Y' *  
    eta <- V %*% crossprod(Y, kappa)
    
    # 3.5 Partition based on   = (b,  )
    eta_b <- eta[1:n_b_params]
    eta_theta <- eta[(n_b_params+1):n_params]
    
    V_bb <- V[1:n_b_params, 1:n_b_params, drop = FALSE]
    V_btheta <- V[1:n_b_params, (n_b_params+1):n_params, drop = FALSE]
    V_theta_theta <- V[(n_b_params+1):n_params, (n_b_params+1):n_params]
    
    # 3.6 Sample   (marginal)
    theta_chol <- chol(V_theta_theta)
    theta <- c(nimble::rmnorm_chol(1, mean = eta_theta,
                                   cholesky = theta_chol, 
                                   prec_param = FALSE))
    
    ## Sample b |  , A, w, z with constraint b_kl <= 0
    
    # 3.7 Conditional distribution
    # Compute V_{  }^{-1} (  -  _ ) and V_{  }^{-1}* V_{ b}using solve()
    # Conditional mean of b:  <-  =  _b + V_{b } * V_{  }^{-1} * (  -  _ )
    b_mean <- eta_b + V_btheta %*% solve(V_theta_theta, theta - eta_theta)
    # Conditional covariance of b: b_cov = V_{bb} - V_{b } * V_{  }^{-1} * V_{ b}
    b_cov <- V_bb - V_btheta %*% solve(V_theta_theta, t(V_btheta))
    # symmetry for numerical stability
    b_cov <- (b_cov + t(b_cov))/2
    
    # 3.8 Joint sampling of b from truncated MVN
    # lower bound = -Inf, upper bound = 0  (elementwise)
    b <- as.numeric(TruncatedNormal::rtmvnorm(
      n = 1,
      mu = b_mean,
      sigma = b_cov, # here sigma is covariance matrix, not percison matrix
      lb = rep(-Inf, n_b_params),
      ub = rep(0,     n_b_params)
    ))
    
    # 3.9 Reconstruct symmetric B matrix
    B <- matrix(0, K, K)
    B[upper.tri(B, diag = FALSE)] <- b
    B <- B + t(B)
    
    ## ===== 4. Update  ,   =====
    
    # extract upper triangle A
    A_upper <- A[upper_idx]
    
    # pair_mask[e, p] = TRUE if edge e belongs to community pair comm_indices[p,]
    pair_mask <- sapply(1:nrow(comm_indices), function(p) {
      (k_vec == comm_indices[p,1]) & (l_vec == comm_indices[p,2])
    }) # this is also column order
    
    ## 4.1 alpha update (A=0)
    A0_mask <- (A_upper == 0)
    
    # For each (k,l): sum observed 1's and 0's over edges where A=0
    # FP = observed 1 but true A=0
    false_pos <- colSums( (A0_mask * Count_1[upper_idx]) * pair_mask )
    # TN = observed 0 and true A=0
    true_neg  <- colSums( (A0_mask * Count_0[upper_idx]) * pair_mask )
    
    # vectorized posterior update (column ordered)
    alpha_vec <- rbeta(nrow(comm_indices), 
                       a_alpha + false_pos, b_alpha + true_neg)
    
    # place  _kl back into K K matrix (column order)
    alpha <- matrix(0, K, K)
    alpha[upper.tri(alpha, diag = TRUE)] <- alpha_vec
    alpha <- alpha + t(alpha) - diag(diag(alpha))
    
    ## 4.2 beta update (A=1) 
    A1_mask <- (A_upper == 1)
    
    # For each (k,l): sum observed 0's and 1's over edges where A=1
    # FN = observed 0 but A=1
    false_neg <- colSums( (A1_mask * Count_0[upper_idx]) * pair_mask )
    # TP = observed 1 and A=1
    true_pos  <- colSums( (A1_mask * Count_1[upper_idx]) * pair_mask )
    
    # vectorized posterior update (column ordered)
    beta_vec <- rbeta(nrow(comm_indices), a_beta + false_neg, b_beta + true_pos)
    
    # place  _kl back into K K matrix
    beta <- matrix(0, K, K)
    beta[upper.tri(beta, diag = TRUE)] <- beta_vec
    beta <- beta + t(beta) - diag(diag(beta))
    
    ## ===== 5. Update continuous covariate parameters =====
    # Calculate community sample sizes Nk
    Nk <- tabulate(z, nbins = K)
    
    if (n_cont > 0) {
      
      # X_cont: n   n_cont
      
      # 5.1 Calculate sum of covariates within each community
      # sum_x is a K   n_cont matrix:
      # rows (k = 1,...,K)     communities
      # columns (r = 1,...,n_cont(R))   continuous covariate dimensions
      # element (k,r):
      # sum_x[k, r] =  _{i=1}^n I{z(i) = k}   x^{cont}_{i r}
      sum_x <- rowsum(X_cont, z)
      
      # 5.2 update  _{k,r}
      # posterior variance:  2_kr / (1/ 0_r + Nk)
      post_var_mu <- tau2 / (outer(rep(1, K), 1 / iota0) + Nk)
      
      # posterior mean:
      #   ( 0/iota0 + sum x) / (1/ 0 + Nk)
      post_mean_mu <- (outer(rep(1, K), mu0 / iota0) + sum_x) / 
        (outer(rep(1, K), 1 / iota0) + Nk)
      
      # sample   (elementwise normal)   = post_mean_mu + post_std_mu * Z
      mu <- post_mean_mu + sqrt(post_var_mu) * matrix(rnorm(K * n_cont), 
                                                      K, n_cont)
      
      # 5.3 update  ^2_{k,r}
      #   (x_ir -  _kr)^2
      # sum_sq is a K   n_cont matrix:
      # rows (k = 1,...,K)          communities
      # columns (r = 1,...,n_cont)   continuous covariate dimensions
      # element (k,r):
      # sum_sq[k, r] =  _{i=1}^n I{ z(i) = k }   ( x^{cont}_{i r}    _{k r} )^2
      sq_dev <- (X_cont - mu[z, ])^2
      sum_sq <- rowsum(sq_dev, z)
      
      # posterior shape: a_ r + (1 + N_k) / 2
      post_shape_tau <- outer(rep(1, K), a_tau) + (1 + Nk) / 2
      
      # posterior scale: b_ r + ( _kr -  _0r) /(2 _0r) + (1/2) (x_ir -  _kr) 
      post_scale_tau <- outer(rep(1, K), b_tau) +
        (mu - outer(rep(1, K), mu0))^2 / 
        (2 * outer(rep(1, K), iota0)) +
        0.5 * sum_sq
      
      # sample  ^2 from IG via 1 / Gamma
      tau2 <- 1 / matrix(
        rgamma(K * n_cont,
               shape = post_shape_tau,
               rate  = post_scale_tau),
        K, n_cont
      )
    }
    
    ## ===== 6. Update categorical covariate parameters =====
    if (n_cat > 0) {
      
      # xi is a list of length n_cat(S):
      # each element xi[[s]] corresponds to the s-th categorical covariate
      # xi[[s]] is a K   F_s matrix:
      # rows (k = 1,...,K)        communities
      # columns (f = 1,...,F_s)   category levels
      xi <- lapply(seq_len(n_cat), function(s) {
        
        # 6.1 community   category counts
        # count_matrix[k, f] =  _{i=1}^n I{ z(i) = k } * I{ X^{cat}_{i s} = f }
        count_matrix <- matrix(
          tabulate(
            interaction(z, X_cat[, s], drop = FALSE),
            nbins = K * F_cat[s]
          ),
          nrow = K,
          ncol = F_cat[s],
          byrow = FALSE
        )
        # For s-th categorical variables, count_matrix is a K   F_s matrix,
        # rows = K(communities), columns = F_s(levels)
        
        # 6.2 Dirichlet posterior parameters
        # Posterior: _{k s} ~ Dirichlet(  _{0 s}+n_{ks1}, ...,  _{0 s}+n_{ksF_s} )
        # n_{ksf} = count_matrix[k, f]
        post_dirichlet <- count_matrix + gamma0[s]
        
        # 6.3 sample  _ks via normalized Gamma
        # Draw: g_{ksf} ~ Gamma(  _{0s} + n_{ksf}, 1 ),  independently
        G <- matrix(
          rgamma(K * F_cat[s], shape = post_dirichlet, rate = 1),
          K, F_cat[s]
        )
        
        #    _{ksf} = g_{ksf} /  _{f=1}^{F_s} g_{ksf}
        G / rowSums(G)
        # For s-th categorical variables, G is also a K   F_s matrix:
        # For each line k, it means in community k, the probabilities of the 
        # s-th categorical variable taking each value of F_s.
      })
    }
    
    ## ===== 7. Update pi =====
    # Dirichlet posterior parameter
    alpha_pi <- d / K + Nk
    
    # Sample pi from Dirichlet via normalized Gamma
    pi_vec <- rgamma(K, shape = alpha_pi, rate = 1)
    pi_vec <- pi_vec / sum(pi_vec)
    
    
    ## ===== 8. Update z =====
    for (i in seq_len(n)){
      zi_old <- z[i]
      idx    <- which(seq_len(n) != i)
      zj     <- z[idx]
      
      # Nk without node i
      Nk_minus_i <- tabulate(zj, nbins = K)
      
      ## 8.1 Observation likelihood 
      # A_i[j] = A_{ij}
      # true edge between node i and node j (j   i)
      A_i  <- A[i, idx] # 1   n
      # C1_i[j] =  {m=1}^M I{ (A~)^{(m)}_{ij} = 1 }
      # number of noisy observations where edge (i,j) is observed as 1
      C1_i <- Count_1[i, idx] # 1   n
      # C0_i[j] =  _{m=1}^M I{ \tilde A^{(m)}_{ij} = 0 }
      # number of noisy observations where edge (i,j) is observed as 0
      C0_i <- Count_0[i, idx] # 1   n
      
      #beta_mat[k, j]  =  _{k z(j)}
      beta_mat  <- beta[, zj, drop = FALSE]   # K   (n-1)
      # alpha_mat[k, j] =  _{k, z(j)}:
      alpha_mat <- alpha[, zj, drop = FALSE]  # K   (n-1)
      
      # Expand vectors to K   (n-1) matrices
      # Each row corresponds to a community k and aligns with 
      # alpha_mat and beta_mat for computation
      A_i_mat  <- outer(rep(1, K), A_i)   # Each row is A_i
      C1_i_mat <- outer(rep(1, K), C1_i)  # Each row is C1_i
      C0_i_mat <- outer(rep(1, K), C0_i)  # Each row is C0_i
      
      # A_i_mat[k,j]      corresponds to A_{ij}
      # C1_i_mat[k,j]     corresponds to C1_{ij}
      # C0_i_mat[k,j]     corresponds to C0_{ij}
      # beta_mat[k,j]     corresponds to  _{k,z(j)}
      # alpha_mat[k,j]    corresponds to  _{k,z(j)}
      # Taking log, the posterior of obs layer is:
      #    _{j i} {
      # A_{ij} [ C1_{ij} log(1  _{k,z(j)}) + C0_{ij} log( _{k,z(j)}) ] +
      # (1 A_{ij}) [ C1_{ij} log( _{k,z(j)}) + C0_{ij} log(1  _{k,z(j)}) ]
      # }.
      log_obs <- rowSums(
        A_i_mat * (C1_i_mat * log(1 - beta_mat) + C0_i_mat * log(beta_mat)) + (1 - A_i_mat) * (C1_i_mat * log(alpha_mat) + C0_i_mat * log(1 - alpha_mat)))
      
      ## 8.2 PG latent term 
      # Expand   to K   (n-1) matrices, theta_mat[k, j] =  _j,
      theta_mat <- outer(rep(1, K), theta[idx]) # Each row is  _(-i)
      #  _ij^(k) =  _i +  _j + B_{k,z(j)}
      psi_k <- theta[i] + theta_mat + B[, zj, drop = FALSE] # K   (n-1)
      
      # Expand W[i, idx]
      W_i_mat <- outer(rep(1, K), W[i, idx]) # K   (n-1)
      
      # I didn't use defined  _ij in Sec 5.3.1 because it only contains upper triangle, here we want all j i
      # Taking log, the posterior of PG layer is:
      #  _{j i} [(A_{ij}   1/2)    _{ij}^{(k)}   0.5w_{ij}   ( _{ij}^{(k)})^2]
      log_pg <- rowSums(
        (A_i_mat - 0.5) * psi_k - 0.5 * W_i_mat * psi_k^2
      )
      
      ## 8.3 Continuous covariates
      log_cont <- if (n_cont > 0) {
        # diff_sq[k,r] = ( x^{cont}_{ir}    _{kr} )^2
        diff_sq <- sweep(mu, 2, X_cont[i, ], "-")^2
        rowSums(-0.5 * log(tau2) - diff_sq / (2 * tau2))
      } else {
        rep(0, K)
      }
      
      ## 8.4 Categorical covariates
      # Taking log, the posterior of categorical layer is:
      #  {s=1}^{n_cat} log(  _{k s X^{cat}_{i s}} )
      log_cat <- if (n_cat > 0) {
        Reduce(`+`, lapply(seq_len(n_cat), function(s) {
          # xi[[s]] is K   F_s, 
          # X_cat[i, s] is the observed level for node i and covariate s.
          
          # xi[[s]][, X_cat[i, s]] extracts a length-K vector:
          # (  _{1 s X^{cat}_{i s}}, ...,  _{K s X^{cat}_{i s}} )
          log(xi[[s]][, X_cat[i, s]])
        }))# Reduce works iff n_cat>=2
      } else {
        rep(0, K)
      }
      
      ## 8.5 Combine and sample
      log_post <- log(pi_vec) + log_obs + log_pg + log_cont + log_cat
      
      # Numerical stability
      post <- exp(log_post - max(log_post))
      post <- post / sum(post)
      
      zi_new <- sample.int(K, 1, prob = post)
      
      ## 8.6 Check and immediately update
      # only check if different community
      if (zi_new != zi_old) {
        # Calculate the size of the community after the update.
        Nk_after <- Nk_minus_i
        Nk_after[zi_new] <- Nk_after[zi_new] + 1
        
        # reject if any Nk <= 1 (all Nk>= 2 hold)
        if (all(Nk_after >= 2)) {
          z[i] <- zi_new  # update z[i]
        }
        # otherwise keep z[i] = zi_old
      }
    }
    
    ## ===== 9. Canonical relabeling =====
    relab <- canonical_relabel(z, K)
    z <- relab$z
    ord   <- relab$ord   
    
    # permute all label-indexed parameters to match new z labels
    pi_vec <- pi_vec[ord]
    
    B     <- B[ord, ord, drop = FALSE]
    alpha <- alpha[ord, ord, drop = FALSE]
    beta  <- beta[ord, ord, drop = FALSE]
    
    if (n_cont > 0) {
      mu   <- mu[ord, , drop = FALSE]
      tau2 <- tau2[ord, , drop = FALSE]
    }
    
    if (n_cat > 0) {
      xi <- lapply(xi, function(mat) mat[ord, , drop = FALSE])
    }
    
    ## ===== 10. Save samples (all iterations, relabeled) =====
    save_idx <- iter
    
    samples_z[save_idx, ]     <- z
    ## ===== 10b. Save pointwise log-likelihood for WAIC (burn-in removed) =====
    if (iter > burn_in) {
      kk <- iter - burn_in  # row index in loglik storage
      
      # ---------- (A) Network layer loglik: p(A_tilde | A, alpha, beta) ----------
      # edge community labels (k,l) for i<j
      k_vec_ll <- z[edge_indices[, 1]]
      l_vec_ll <- z[edge_indices[, 2]]
      
      # alpha/beta for each edge, clipped for stability
      alpha_edges <- clip01(alpha[cbind(k_vec_ll, l_vec_ll)])
      beta_edges  <- clip01(beta[cbind(k_vec_ll, l_vec_ll)])
      
      # latent A on upper triangle
      A_edges  <- A[upper_idx]
      C1_edges <- Count_1[upper_idx]
      C0_edges <- Count_0[upper_idx]
      
      log_p_A1 <- C1_edges * log(1 - beta_edges) + C0_edges * log(beta_edges)
      log_p_A0 <- C1_edges * log(alpha_edges)    + C0_edges * log(1 - alpha_edges)
      
      samples_loglik_net[kk, ] <- A_edges * log_p_A1 + (1 - A_edges) * log_p_A0
      
      # ---------- (B) Continuous covariates loglik: p(X_cont | z, mu, tau2) ----------
      if (n_cont > 0) {
        mu_nodes  <- mu[z, , drop = FALSE]                       # n   n_cont
        sd_nodes  <- sqrt(pmax(tau2[z, , drop = FALSE], eps))     # n   n_cont
        
        samples_loglik_cont[kk, ] <- rowSums(
          dnorm(X_cont, mean = mu_nodes, sd = sd_nodes, log = TRUE)
        )
      }
      
      # ---------- (C) Categorical covariates loglik: p(X_cat | z, xi) ----------
      if (n_cat > 0) {
        ll_cat_vec <- rep(0, n)
        for (s in seq_len(n_cat)) {
          # xi[[s]]: K   F_s
          probs <- xi[[s]][cbind(z, X_cat[, s])]
          ll_cat_vec <- ll_cat_vec + log(pmax(probs, eps))
        }
        samples_loglik_cat[kk, ] <- ll_cat_vec
      }
    }
  }
  
  
  return(list(
    z     = samples_z,
    burn_in = burn_in,
    loglik_net  = samples_loglik_net,
    loglik_cont = samples_loglik_cont,
    loglik_cat  = samples_loglik_cat
  ))
}



## ---- Canonical relabel function (Peng & Carvalho, 2016) ----
canonical_relabel <- function(z, K) {
  n <- length(z)
  first_idx <- rep(Inf, K)
  for (i in 1:n) {
    lab <- z[i]
    if (is.infinite(first_idx[lab])) {
      first_idx[lab] <- i
    }
  }
  ord  <- order(first_idx)      # ord[k] = old label that becomes new label k
  z_new <- match(z, ord)        # new label = position of old label in ord
  list(z = z_new, ord = ord)
}



## ---- final ----
RobustCDC <- function(
    A_tilde,
    X_cont = NULL,
    X_cat  = NULL,
    K = 2:10,
    n_iter = 3000,
    burn_in = 500,
    hyper,
    init = NULL,
    verbose = TRUE
) {
  
  ## ---- Input check ----
  if (!is.numeric(K) || any(K < 2)) {
    stop("K must be >=2")
  }
  
  K <- sort(unique(as.integer(K)))
  n_K <- length(K)
  
  if (verbose) {
    
    cat("====================================\n")
    cat("RCDC: Selecting K via WAIC\n")
    
    if (n_K==1){
      cat(sprintf("Fitting fixed K=%d\n",K))
    }else{
      cat(sprintf("Testing K values: %s\n",
                  paste(K,collapse=", ")))
    }
    
    cat(sprintf("MCMC: %d iterations (%d burn-in)\n",
                n_iter,burn_in))
    
    cat("====================================\n\n")
  }
  
  
  ## table for reporting only (very small memory)
  waic_values <- numeric(n_K)
  
  
  
  ## ---- BEST HOLDER ----
  best_fit <- NULL
  best_K <- NULL
  best_waic <- Inf
  
  
  ## ---- LOOP ----
  for (i in seq_along(K)) {
    
    Ki <- K[i]
    
    if (verbose)
      cat("Running K =",Ki,"\n")
    
    
    ## 1) Run MCMC
    fit <- robustCDC_MCMC(
      
      A_tilde = A_tilde,
      X_cont  = X_cont,
      X_cat   = X_cat,
      K       = Ki,
      n_iter  = n_iter,
      burn_in = burn_in,
      hyper   = hyper,
      init    = init
    )
    
    
    ## 2) WAIC
    log_lik_total <- cbind(
      
      fit$loglik_net,
      
      if(!is.null(fit$loglik_cont))
        fit$loglik_cont,
      
      if(!is.null(fit$loglik_cat))
        fit$loglik_cat
    )
    
    if(anyNA(log_lik_total))
      stop("NA detected in loglik.")
    
    waic_res <- suppressWarnings(
      loo::waic(log_lik_total)
    )
    
    waic_val <- waic_res$estimates["waic","Estimate"]
    waic_values[i] <- waic_val
    
    ## 3) Compare with BEST
    if (waic_val < best_waic){
      
      if(verbose)
        cat(" -> New BEST model (lower WAIC)\n")
      
      best_waic <- waic_val
      best_K <- Ki
      
      z_post <- fit$z[(burn_in+1):n_iter, ]
      z_hat <- apply(
        z_post,
        2,
        function(x) as.integer(names(which.max(table(x))))
      )
      rm(z_post)
      
    }else{
      
      if(verbose)
        cat(" -> Discarded (worse WAIC)\n")
    }
    rm(log_lik_total)
    rm(waic_res)
    rm(fit)
    gc(FALSE)
  }
  
  
  ## ---- Summary ----
  comparison <- data.frame(
    
    K = K,
    WAIC = waic_values,
    Delta_WAIC = waic_values - min(waic_values)
    
  )
  
  if(verbose){
    
    cat("\n=============================\n")
    cat("Selected K =",best_K,"\n")
    
    print(comparison)
    
    cat("=============================\n")
    
  }
  
  
  ## ---- Output ----
  out <- list(
    K = best_K,
    WAIC = best_waic,
    z = z_hat,
    waic_table = comparison
    
  )
  
  class(out) <- "robustCDC"
  
  return(out)
}



## ---- main function (save all parameters) ----
robustCDC_MCMC_saveallparameters <- function(
    A_tilde,     
    X_cont = NULL,
    X_cat  = NULL,
    K,
    n_iter = 3000,
    burn_in = 500,
    hyper,
    init = NULL # list of hyperparameters
) {
  ## ------------------------------------------------
  ## 0. validate dimensions
  ## ------------------------------------------------
  if (length(dim(A_tilde)) != 3) stop("A_tilde must be a 3D array n   n   M")
  
  n <- dim(A_tilde)[1]
  M <- dim(A_tilde)[3]  
  
  if (dim(A_tilde)[2] != n) stop("A_tilde must be n   n   M")
  
  n_cont <- if (is.null(X_cont)) 0 else ncol(X_cont)
  n_cat  <- if (is.null(X_cat))  0 else ncol(X_cat)
  
  ## ------------------------------------------------
  ## 1. hyperparameters
  ## ------------------------------------------------
  sigma       <- hyper$sigma
  
  a_alpha <- hyper$a_alpha
  b_alpha <- hyper$b_alpha  
  a_beta  <- hyper$a_beta
  b_beta  <- hyper$b_beta   
  
  # for continuous covariates
  mu0     <- hyper$mu0      
  iota0   <- hyper$iota0    
  a_tau   <- hyper$a_tau    
  b_tau   <- hyper$b_tau    
  
  # for categorical covariates
  F_cat   <- hyper$F_cat    
  gamma0  <- hyper$gamma0   
  
  # for mixing proportions  
  d   <- hyper$d    
  
  
  ## ------------------------------------------------
  ## 2. Pre-compute constants (outside MCMC loop)
  ## ------------------------------------------------
  # Observation statistics (constant across iterations)
  Count_1 <- rowSums(A_tilde, dims = 2)  # sum_m I(A_tilde^(m)_ij = 1)
  Count_0 <- M - Count_1                  # sum_m I(A_tilde^(m)_ij = 0)
  
  # Upper triangle indices 
  upper_idx <- which(upper.tri(matrix(0, n, n), diag = FALSE)) # column order
  # include (i<j): column-major order
  edge_indices <- arrayInd(upper_idx, .dim = c(n, n)) # column order
  # include (k=k) and (k<l): column-major order
  comm_indices <- which(upper.tri(matrix(0, K, K), diag = TRUE), arr.ind = TRUE)
  
  # Design matrix Y for   and B sampling
  # Each row y_ij corresponds to edge (i,j) with i < j
  # Dimension: n(n-1)/2   (K(K-1)/2 + n)
  # Note: delta = (b,  ), B only has off-diagonal entries (B_kk = 0)
  n_edges <- length(upper_idx)
  n_b_params <- K * (K - 1) / 2  # Only different community pairs
  n_params <- n_b_params + n # length of all   + B
  Y <- matrix(0, n_edges, n_params)
  
  # Fill   part (columns n_b_params+1 : n_params)
  Y[cbind(1:n_edges, n_b_params + edge_indices[, 1])] <- 1
  Y[cbind(1:n_edges, n_b_params + edge_indices[, 2])] <- 1
  #The B-part (columns 1 : n_b_params) depends on z, will fill inside loop
  
  ## ------------------------------------------------
  ## 3. Initialization parameters
  ## ------------------------------------------------
  
  if (is.null(init)) {
    ## ----- Initialize A and z (le2018estimating) ----
    A       <- (Count_1 >= (M / 2)) * 1L
    if (any(rowSums(A)==0)) {iso <- which(rowSums(A)==0); n <- nrow(A); j <- ifelse(iso==n,1,iso+1); A[cbind(iso,j)] <- 1; A[cbind(j,iso)] <- 1}
    z <- as.numeric(
      suppressMessages(
        C4(A, K = K)
      )$cluster
    )
    
    z       <- canonical_relabel(z,K)$z
    ## ----- Initialize theta and B (Le-style parameter mapping) -----
    Psi_matrix     <- matrix(0, K, K)
    edge_comm_idx <- cbind(pmin(z[edge_indices[,1]], z[edge_indices[,2]]), 
                           pmax(z[edge_indices[,1]], z[edge_indices[,2]]))
    
    # Estimate Psi_kl based on empirical edge density in each block
    for (p in 1:nrow(comm_indices)) {
      k <- comm_indices[p, 1]
      l <- comm_indices[p, 2]
      
      # Logical mask for the current block (k, l)
      block_mask <- (edge_comm_idx[,1] == k & edge_comm_idx[,2] == l)
      
      # Calculate empirical edge probability (w_hat) from majority vote matrix A 
      w_hat <- mean(A[upper_idx][block_mask])
      if (is.nan(w_hat) || !any(block_mask)) {
        val <- 0
      } else {
        w_hat <- pmin(pmax(w_hat, 1e-6), 1 - 1e-6)
        val <- qlogis(w_hat)
      }
      # Map to DCSBM parameter Psi_kl using logit link
      val <- qlogis(w_hat)
      Psi_matrix[k, l] <- Psi_matrix[l, k] <- val
    }
    
    # 1) community-level theta from diagonal (uses B_kk = 0 constraint)
    theta_block <- 0.5 * diag(Psi_matrix)     # length K, indexed by community label 1..K
    
    # 2) node-level theta aligned to node order (1..n) using z
    theta <- theta_block[z]                   # length n, indexed by node i
    
    # 3) community-level B reconstructed from psi = theta_k + theta_l + B_kl
    B <- Psi_matrix - outer(theta_block, theta_block, "+")
    
    ## ----- Initialize alpha and beta (le2018estimating) -----
    alpha <- matrix(0, K, K)
    beta  <- matrix(0, K, K)
    for (p in 1:nrow(comm_indices)) {
      k <- comm_indices[p, 1]
      l <- comm_indices[p, 2]
      
      block_mask <- (edge_comm_idx[,1] == k & edge_comm_idx[,2] == l)
      # Extract relevant statistics for the current block
      S_b <- Count_1[upper_idx][block_mask] # Observed counts of 1's
      A_b <- A[upper_idx][block_mask]       # Current latent adjacency A
      
      # Update alpha (False Positive rate): P(observed=1 | A=0) 
      # Sum of observed 1's where A=0 divided by total possible observations
      idx_A0 <- (A_b == 0)
      if (sum(idx_A0) > 0) {
        alpha[k,l] <- alpha[l,k] <- sum(S_b[idx_A0]) / (M * sum(idx_A0))
      } else {
        # Fallback to prior mean if the block contains no A=1 edges (e.g., 2 communities have no connection)
        alpha[k,l] <- alpha[l,k] <- a_alpha / (a_alpha + b_alpha)
      }
      
      # Update beta (False Negative rate): P(observed=0 | A=1)
      # Sum of observed 0's where A=1 divided by total possible observations
      idx_A1 <- (A_b == 1)
      if (sum(idx_A1) > 0) {
        beta[k,l]  <- beta[l,k]  <- sum(M - S_b[idx_A1]) / (M * sum(idx_A1))
      } else {
        # Fallback to prior mean
        beta[k,l]  <- beta[l,k]  <- a_beta / (a_beta + b_beta)
      }
    }
    
    ## ----- Initialize covariate -----
    if (n_cont > 0) {
      
      #  _{k } = average of X within each community
      mu <- rowsum(X_cont, z) / tabulate(z, nbins = K)
      
      #   _{k } = average squared deviation within each community
      tau2 <- rowsum( (X_cont - mu[z, ])^2 , z) / tabulate(z, nbins = K)
      
    } else {
      mu   <- NULL
      tau2 <- NULL
    }
    
    if (n_cat > 0) {
      
      xi <- vector("list", n_cat)
      
      for (s in 1:n_cat) {
        
        # count_{k,f} = number of nodes in community k with category f
        count_matrix <- matrix(
          tabulate(interaction(z, X_cat[, s]), nbins = K * F_cat[s]),
          nrow = K, ncol = F_cat[s], byrow = FALSE
        )
        
        #  _{ks } = normalized frequencies in each community
        xi[[s]] <- count_matrix / rowSums(count_matrix)
      }
      
    } else {
      xi <- NULL
    }
    
  } else {
    
    ## ----- user-supplied initialization -----
    A      <- init$A
    z      <- init$z
    theta  <- init$theta
    B      <- init$B
    alpha  <- init$alpha
    beta   <- init$beta
    
    mu     <- if (n_cont > 0) init$mu   else NULL
    tau2   <- if (n_cont > 0) init$tau2 else NULL
    xi     <- if (n_cat  > 0) init$xi   else NULL
  }
  
  
  ## ------------------------------------------------
  ## 4. storage container
  ## ------------------------------------------------
  n_save <- n_iter
  
  samples_A     <- array(NA, dim = c(n_save, n, n))
  samples_z     <- matrix(NA, n_save, n)
  samples_pi    <- matrix(NA, n_save, K)
  samples_theta <- matrix(NA, n_save, n)
  samples_psi   <- array(NA, dim = c(n_save, n, n))
  samples_B     <- array(NA, dim = c(n_save, K, K))
  samples_alpha <- array(NA, dim = c(n_save, K, K))
  samples_beta  <- array(NA, dim = c(n_save, K, K))
  
  # continuous covariates
  if (n_cont > 0) {
    samples_mu   <- array(NA, dim = c(n_save, K, n_cont))
    samples_tau2 <- array(NA, dim = c(n_save, K, n_cont))
  } else {
    samples_mu   <- NULL
    samples_tau2 <- NULL
  }
  
  # categorical covariates
  if (n_cat > 0) {
    samples_xi <- vector("list", n_cat)
    for (j in 1:n_cat) {
      samples_xi[[j]] <- array(NA, dim = c(n_save, K, F_cat[j]))
    }
  } else {
    samples_xi <- NULL
  }
  
  keep_idx <- (burn_in + 1L):n_iter
  n_keep <- length(keep_idx)
  
  samples_loglik_net <- matrix(NA_real_, nrow = n_keep, ncol = n_edges)
  
  
  samples_loglik_cont <- if (n_cont > 0) matrix(0, nrow = n_keep, ncol = n) else NULL
  samples_loglik_cat  <- if (n_cat  > 0) matrix(0, nrow = n_keep, ncol = n) else NULL
  
  eps <- 1e-12
  clip01 <- function(p) pmin(pmax(p, eps), 1 - eps)
  
  ## ------------------------------------------------
  ## 5. MCMC 
  ## ------------------------------------------------
  for (iter in 1:n_iter) {
    
    ## ===== 1. Update A =====
    # 1.1 Construct Linear Predictor Psi (Prior Log-Odds)
    #  _ij =  _i +  _j + B_{z_i, z_j}
    B_mat <- B[z, z]                    # Map B matrix to node pairs
    Psi <- outer(theta, theta, FUN = `+`) + B_mat # Do not need indicator as long as B has 0 diagnal
    diag(Psi) <- 0                      # No self-loops
    
    # 1.2 Error rates 
    # All map to node pairs
    Alpha_mat <- alpha[z, z]
    Beta_mat  <- beta[z, z]
    
    # 1.3 Calculate Log-Likelihoods
    # p*_{ij,1} = exp( )   (1- )^n1    ^n0
    p_star_1 <- exp(Psi) * (1 - Beta_mat)^Count_1 * Beta_mat^Count_0
    # p*_{ij,0} =  ^n1   (1- )^n0
    p_star_0 <- Alpha_mat^Count_1 * (1 - Alpha_mat)^Count_0
    
    # 1.4 Compute Posterior Probability P(A_ij=1 | Data)
    Prob_A1 <- p_star_1 / (p_star_1 + p_star_0)
    
    # 1.5 Sample A Efficiently (Upper Triangle Only for Undirected Graph)
    # Initialize new adjacency matrix
    A_new <- matrix(0, n, n)
    # Sample Bernoulli for upper triangle (vectorized)
    A_new[upper_idx] <- rbinom(n_edges, 1, Prob_A1[upper_idx])
    # Symmetrize: copy upper triangle to lower triangle
    A_new <- A_new + t(A_new)
    # Update A
    A <- A_new
    
    ## ===== 2. Update W (PG) =====
    # Sample Polya-Gamma augmentation variables W_{ij} ~ PG(1,  _{ij})
    # Construct symmetric W matrix
    W_new <- matrix(0, n, n)
    W_new[upper_idx] <- pgdraw::pgdraw(rep(1, n_edges), Psi[upper_idx])
    W_new <- W_new + t(W_new)
    # Update W
    W <- W_new
    
    ## ===== 3. Update   and B (joint sampling) =====
    # 3.1 Construct kappa vector:  _ij = A_ij - 1/2
    kappa <- (A - 0.5)[upper_idx]
    
    # 3.2 Update Y matrix B-part (depends on current z)
    # Reset B columns to 0
    Y[, 1:n_b_params] <- 0
    
    # Only edges between DIFFERENT communities involve B parameters
    # First get all node pair community assignment(k,l) and let k <=l
    k_vec <- pmin(z[edge_indices[, 1]], z[edge_indices[, 2]])
    l_vec <- pmax(z[edge_indices[, 1]], z[edge_indices[, 2]])
    # Only fill for k != l (inter-community edges)
    is_inter <- k_vec != l_vec
    
    if (any(is_inter)) {
      # Map (k, l) to linear index in upper triangle of B (k < l, excluding diagonal)
      # left part of Y follows standard R "upper.tri" order (Column-Major): 
      # (1,2), (1,3), (2,3), (1,4)...
      # Formula: (l-1)*(l-2)/2 + k
      
      b_idx_vec <- (l_vec[is_inter] - 1) * (l_vec[is_inter] - 2) / 
        2 + k_vec[is_inter]
      
      # Set 1s in Y
      Y[cbind(which(is_inter), b_idx_vec)] <- 1
    }
    # By zeyu: Check all here, correct, follows formula in manuscript:  =Y 
    # 3.3 Construct diagonal weight matrix  
    Omega <- W[upper_idx]
    
    # 3.4 Posterior covariance and mean for   = (b,  )
    # V = [Y'   Y + (1/  )I]^{-1}
    # V^-1 = Y'   Y + (1/  )I
    V_inv <- crossprod(Y * sqrt(Omega)) + (1/sigma^2) * diag(n_params)
    
    # Cholesky decomposition for stable inversion(better than solve)
    # Decompose the symmetric positive-definite matrix V^-1 into 
    # V_inv_chol' * V_inv_chol
    V_inv_chol <- chol(V_inv)
    # Use V_inv_chol to compute V, stable and efficient
    V <- chol2inv(V_inv_chol)
    
    # Posterior mean   = V * Y' *  
    eta <- V %*% crossprod(Y, kappa)
    
    # 3.5 Partition based on   = (b,  )
    eta_b <- eta[1:n_b_params]
    eta_theta <- eta[(n_b_params+1):n_params]
    
    V_bb <- V[1:n_b_params, 1:n_b_params, drop = FALSE]
    V_btheta <- V[1:n_b_params, (n_b_params+1):n_params, drop = FALSE]
    V_theta_theta <- V[(n_b_params+1):n_params, (n_b_params+1):n_params]
    
    # 3.6 Sample   (marginal)
    theta_chol <- chol(V_theta_theta)
    theta <- c(nimble::rmnorm_chol(1, mean = eta_theta,
                                   cholesky = theta_chol, 
                                   prec_param = FALSE))
    
    ## Sample b |  , A, w, z with constraint b_kl <= 0
    
    # 3.7 Conditional distribution
    # Compute V_{  }^{-1} (  -  _ ) and V_{  }^{-1}* V_{ b}using solve()
    # Conditional mean of b:  <-  =  _b + V_{b } * V_{  }^{-1} * (  -  _ )
    b_mean <- eta_b + V_btheta %*% solve(V_theta_theta, theta - eta_theta)
    # Conditional covariance of b: b_cov = V_{bb} - V_{b } * V_{  }^{-1} * V_{ b}
    b_cov <- V_bb - V_btheta %*% solve(V_theta_theta, t(V_btheta))
    # symmetry for numerical stability
    b_cov <- (b_cov + t(b_cov))/2
    
    # 3.8 Joint sampling of b from truncated MVN
    # lower bound = -Inf, upper bound = 0  (elementwise)
    b <- as.numeric(TruncatedNormal::rtmvnorm(
      n = 1,
      mu = b_mean,
      sigma = b_cov, # here sigma is covariance matrix, not percison matrix
      lb = rep(-Inf, n_b_params),
      ub = rep(0,     n_b_params)
    ))
    
    # 3.9 Reconstruct symmetric B matrix
    B <- matrix(0, K, K)
    B[upper.tri(B, diag = FALSE)] <- b
    B <- B + t(B)
    
    ## ===== 4. Update  ,   =====
    
    # extract upper triangle A
    A_upper <- A[upper_idx]
    
    # pair_mask[e, p] = TRUE if edge e belongs to community pair comm_indices[p,]
    pair_mask <- sapply(1:nrow(comm_indices), function(p) {
      (k_vec == comm_indices[p,1]) & (l_vec == comm_indices[p,2])
    }) # this is also column order
    
    ## 4.1 alpha update (A=0)
    A0_mask <- (A_upper == 0)
    
    # For each (k,l): sum observed 1's and 0's over edges where A=0
    # FP = observed 1 but true A=0
    false_pos <- colSums( (A0_mask * Count_1[upper_idx]) * pair_mask )
    # TN = observed 0 and true A=0
    true_neg  <- colSums( (A0_mask * Count_0[upper_idx]) * pair_mask )
    
    # vectorized posterior update (column ordered)
    alpha_vec <- rbeta(nrow(comm_indices), 
                       a_alpha + false_pos, b_alpha + true_neg)
    
    # place  _kl back into K K matrix (column order)
    alpha <- matrix(0, K, K)
    alpha[upper.tri(alpha, diag = TRUE)] <- alpha_vec
    alpha <- alpha + t(alpha) - diag(diag(alpha))
    
    ## 4.2 beta update (A=1) 
    A1_mask <- (A_upper == 1)
    
    # For each (k,l): sum observed 0's and 1's over edges where A=1
    # FN = observed 0 but A=1
    false_neg <- colSums( (A1_mask * Count_0[upper_idx]) * pair_mask )
    # TP = observed 1 and A=1
    true_pos  <- colSums( (A1_mask * Count_1[upper_idx]) * pair_mask )
    
    # vectorized posterior update (column ordered)
    beta_vec <- rbeta(nrow(comm_indices), a_beta + false_neg, b_beta + true_pos)
    
    # place  _kl back into K K matrix
    beta <- matrix(0, K, K)
    beta[upper.tri(beta, diag = TRUE)] <- beta_vec
    beta <- beta + t(beta) - diag(diag(beta))
    
    ## ===== 5. Update continuous covariate parameters =====
    # Calculate community sample sizes Nk
    Nk <- tabulate(z, nbins = K)
    
    if (n_cont > 0) {
      
      # X_cont: n   n_cont
      
      # 5.1 Calculate sum of covariates within each community
      # sum_x is a K   n_cont matrix:
      # rows (k = 1,...,K)     communities
      # columns (r = 1,...,n_cont(R))   continuous covariate dimensions
      # element (k,r):
      # sum_x[k, r] =  _{i=1}^n I{z(i) = k}   x^{cont}_{i r}
      sum_x <- rowsum(X_cont, z)
      
      # 5.2 update  _{k,r}
      # posterior variance:  2_kr / (1/ 0_r + Nk)
      post_var_mu <- tau2 / (outer(rep(1, K), 1 / iota0) + Nk)
      
      # posterior mean:
      #   ( 0/iota0 + sum x) / (1/ 0 + Nk)
      post_mean_mu <- (outer(rep(1, K), mu0 / iota0) + sum_x) / 
        (outer(rep(1, K), 1 / iota0) + Nk)
      
      # sample   (elementwise normal)   = post_mean_mu + post_std_mu * Z
      mu <- post_mean_mu + sqrt(post_var_mu) * matrix(rnorm(K * n_cont), 
                                                      K, n_cont)
      
      # 5.3 update  ^2_{k,r}
      #   (x_ir -  _kr)^2
      # sum_sq is a K   n_cont matrix:
      # rows (k = 1,...,K)          communities
      # columns (r = 1,...,n_cont)   continuous covariate dimensions
      # element (k,r):
      # sum_sq[k, r] =  _{i=1}^n I{ z(i) = k }   ( x^{cont}_{i r}    _{k r} )^2
      sq_dev <- (X_cont - mu[z, ])^2
      sum_sq <- rowsum(sq_dev, z)
      
      # posterior shape: a_ r + (1 + N_k) / 2
      post_shape_tau <- outer(rep(1, K), a_tau) + (1 + Nk) / 2
      
      # posterior scale: b_ r + ( _kr -  _0r) /(2 _0r) + (1/2) (x_ir -  _kr) 
      post_scale_tau <- outer(rep(1, K), b_tau) +
        (mu - outer(rep(1, K), mu0))^2 / 
        (2 * outer(rep(1, K), iota0)) +
        0.5 * sum_sq
      
      # sample  ^2 from IG via 1 / Gamma
      tau2 <- 1 / matrix(
        rgamma(K * n_cont,
               shape = post_shape_tau,
               rate  = post_scale_tau),
        K, n_cont
      )
    }
    
    ## ===== 6. Update categorical covariate parameters =====
    if (n_cat > 0) {
      
      # xi is a list of length n_cat(S):
      # each element xi[[s]] corresponds to the s-th categorical covariate
      # xi[[s]] is a K   F_s matrix:
      # rows (k = 1,...,K)        communities
      # columns (f = 1,...,F_s)   category levels
      xi <- lapply(seq_len(n_cat), function(s) {
        
        # 6.1 community   category counts
        # count_matrix[k, f] =  _{i=1}^n I{ z(i) = k } * I{ X^{cat}_{i s} = f }
        count_matrix <- matrix(
          tabulate(
            interaction(z, X_cat[, s], drop = FALSE),
            nbins = K * F_cat[s]
          ),
          nrow = K,
          ncol = F_cat[s],
          byrow = FALSE
        )
        # For s-th categorical variables, count_matrix is a K   F_s matrix,
        # rows = K(communities), columns = F_s(levels)
        
        # 6.2 Dirichlet posterior parameters
        # Posterior: _{k s} ~ Dirichlet(  _{0 s}+n_{ks1}, ...,  _{0 s}+n_{ksF_s} )
        # n_{ksf} = count_matrix[k, f]
        post_dirichlet <- count_matrix + gamma0[s]
        
        # 6.3 sample  _ks via normalized Gamma
        # Draw: g_{ksf} ~ Gamma(  _{0s} + n_{ksf}, 1 ),  independently
        G <- matrix(
          rgamma(K * F_cat[s], shape = post_dirichlet, rate = 1),
          K, F_cat[s]
        )
        
        #    _{ksf} = g_{ksf} /  _{f=1}^{F_s} g_{ksf}
        G / rowSums(G)
        # For s-th categorical variables, G is also a K   F_s matrix:
        # For each line k, it means in community k, the probabilities of the 
        # s-th categorical variable taking each value of F_s.
      })
    }
    
    ## ===== 7. Update pi =====
    # Dirichlet posterior parameter
    alpha_pi <- d / K + Nk
    
    # Sample pi from Dirichlet via normalized Gamma
    pi_vec <- rgamma(K, shape = alpha_pi, rate = 1)
    pi_vec <- pi_vec / sum(pi_vec)
    
    
    ## ===== 8. Update z =====
    for (i in seq_len(n)){
      zi_old <- z[i]
      idx    <- which(seq_len(n) != i)
      zj     <- z[idx]
      
      # Nk without node i
      Nk_minus_i <- tabulate(zj, nbins = K)
      
      ## 8.1 Observation likelihood 
      # A_i[j] = A_{ij}
      # true edge between node i and node j (j   i)
      A_i  <- A[i, idx] # 1   n
      # C1_i[j] =  {m=1}^M I{ (A~)^{(m)}_{ij} = 1 }
      # number of noisy observations where edge (i,j) is observed as 1
      C1_i <- Count_1[i, idx] # 1   n
      # C0_i[j] =  _{m=1}^M I{ \tilde A^{(m)}_{ij} = 0 }
      # number of noisy observations where edge (i,j) is observed as 0
      C0_i <- Count_0[i, idx] # 1   n
      
      #beta_mat[k, j]  =  _{k z(j)}
      beta_mat  <- beta[, zj, drop = FALSE]   # K   (n-1)
      # alpha_mat[k, j] =  _{k, z(j)}:
      alpha_mat <- alpha[, zj, drop = FALSE]  # K   (n-1)
      
      # Expand vectors to K   (n-1) matrices
      # Each row corresponds to a community k and aligns with 
      # alpha_mat and beta_mat for computation
      A_i_mat  <- outer(rep(1, K), A_i)   # Each row is A_i
      C1_i_mat <- outer(rep(1, K), C1_i)  # Each row is C1_i
      C0_i_mat <- outer(rep(1, K), C0_i)  # Each row is C0_i
      
      # A_i_mat[k,j]      corresponds to A_{ij}
      # C1_i_mat[k,j]     corresponds to C1_{ij}
      # C0_i_mat[k,j]     corresponds to C0_{ij}
      # beta_mat[k,j]     corresponds to  _{k,z(j)}
      # alpha_mat[k,j]    corresponds to  _{k,z(j)}
      # Taking log, the posterior of obs layer is:
      #    _{j i} {
      # A_{ij} [ C1_{ij} log(1  _{k,z(j)}) + C0_{ij} log( _{k,z(j)}) ] +
      # (1 A_{ij}) [ C1_{ij} log( _{k,z(j)}) + C0_{ij} log(1  _{k,z(j)}) ]
      # }.
      log_obs <- rowSums(
        A_i_mat * (C1_i_mat * log(1 - beta_mat) + C0_i_mat * log(beta_mat)) + (1 - A_i_mat) * (C1_i_mat * log(alpha_mat) + C0_i_mat * log(1 - alpha_mat)))
      
      ## 8.2 PG latent term 
      # Expand   to K   (n-1) matrices, theta_mat[k, j] =  _j,
      theta_mat <- outer(rep(1, K), theta[idx]) # Each row is  _(-i)
      #  _ij^(k) =  _i +  _j + B_{k,z(j)}
      psi_k <- theta[i] + theta_mat + B[, zj, drop = FALSE] # K   (n-1)
      
      # Expand W[i, idx]
      W_i_mat <- outer(rep(1, K), W[i, idx]) # K   (n-1)
      
      # I didn't use defined  _ij in Sec 5.3.1 because it only contains upper triangle, here we want all j i
      # Taking log, the posterior of PG layer is:
      #  _{j i} [(A_{ij}   1/2)    _{ij}^{(k)}   0.5w_{ij}   ( _{ij}^{(k)})^2]
      log_pg <- rowSums(
        (A_i_mat - 0.5) * psi_k - 0.5 * W_i_mat * psi_k^2
      )
      
      ## 8.3 Continuous covariates
      log_cont <- if (n_cont > 0) {
        # diff_sq[k,r] = ( x^{cont}_{ir}    _{kr} )^2
        diff_sq <- sweep(mu, 2, X_cont[i, ], "-")^2
        rowSums(-0.5 * log(tau2) - diff_sq / (2 * tau2))
      } else {
        rep(0, K)
      }
      
      ## 8.4 Categorical covariates
      # Taking log, the posterior of categorical layer is:
      #  {s=1}^{n_cat} log(  _{k s X^{cat}_{i s}} )
      log_cat <- if (n_cat > 0) {
        Reduce(`+`, lapply(seq_len(n_cat), function(s) {
          # xi[[s]] is K   F_s, 
          # X_cat[i, s] is the observed level for node i and covariate s.
          
          # xi[[s]][, X_cat[i, s]] extracts a length-K vector:
          # (  _{1 s X^{cat}_{i s}}, ...,  _{K s X^{cat}_{i s}} )
          log(xi[[s]][, X_cat[i, s]])
        }))# Reduce works iff n_cat>=2
      } else {
        rep(0, K)
      }
      
      ## 8.5 Combine and sample
      log_post <- log(pi_vec) + log_obs + log_pg + log_cont + log_cat
      
      # Numerical stability
      post <- exp(log_post - max(log_post))
      post <- post / sum(post)
      
      zi_new <- sample.int(K, 1, prob = post)
      
      ## 8.6 Check and immediately update
      # only check if different community
      if (zi_new != zi_old) {
        # Calculate the size of the community after the update.
        Nk_after <- Nk_minus_i
        Nk_after[zi_new] <- Nk_after[zi_new] + 1
        
        # reject if any Nk <= 1 (all Nk>= 2 hold)
        if (all(Nk_after >= 2)) {
          z[i] <- zi_new  # update z[i]
        }
        # otherwise keep z[i] = zi_old
      }
    }
    
    ## ===== 9. Canonical relabeling =====
    relab <- canonical_relabel(z, K)
    z <- relab$z
    ord   <- relab$ord   
    
    # permute all label-indexed parameters to match new z labels
    pi_vec <- pi_vec[ord]
    
    B     <- B[ord, ord, drop = FALSE]
    alpha <- alpha[ord, ord, drop = FALSE]
    beta  <- beta[ord, ord, drop = FALSE]
    
    if (n_cont > 0) {
      mu   <- mu[ord, , drop = FALSE]
      tau2 <- tau2[ord, , drop = FALSE]
    }
    
    if (n_cat > 0) {
      xi <- lapply(xi, function(mat) mat[ord, , drop = FALSE])
    }
    
    ## ===== 10. Save samples (all iterations, relabeled) =====
    save_idx <- iter
    
    samples_A[save_idx, , ]   <- A
    samples_z[save_idx, ]     <- z
    samples_pi[save_idx, ]    <- pi_vec
    samples_theta[save_idx, ] <- theta
    samples_B[save_idx, , ]   <- B
    samples_alpha[save_idx, , ] <- alpha
    samples_beta[save_idx, , ]  <- beta
    
    # save   using FINAL ( , B, z) 
    psi_save <- outer(theta, theta, "+") + B[z, z]
    diag(psi_save) <- 0
    samples_psi[save_idx, , ] <- psi_save
    
    # continuous covariates
    if (n_cont > 0) {
      samples_mu[save_idx, , ]   <- mu
      samples_tau2[save_idx, , ] <- tau2
    }
    
    # categorical covariates
    if (n_cat > 0) {
      for (j in 1:n_cat) {
        samples_xi[[j]][save_idx, , ] <- xi[[j]]
      }
    }
    ## ===== 10b. Save pointwise log-likelihood for WAIC (burn-in removed) =====
    if (iter > burn_in) {
      kk <- iter - burn_in  # row index in loglik storage
      
      # ---------- (A) Network layer loglik: p(A_tilde | A, alpha, beta) ----------
      # edge community labels (k,l) for i<j
      k_vec_ll <- z[edge_indices[, 1]]
      l_vec_ll <- z[edge_indices[, 2]]
      
      # alpha/beta for each edge, clipped for stability
      alpha_edges <- clip01(alpha[cbind(k_vec_ll, l_vec_ll)])
      beta_edges  <- clip01(beta[cbind(k_vec_ll, l_vec_ll)])
      
      # latent A on upper triangle
      A_edges  <- A[upper_idx]
      C1_edges <- Count_1[upper_idx]
      C0_edges <- Count_0[upper_idx]
      
      log_p_A1 <- C1_edges * log(1 - beta_edges) + C0_edges * log(beta_edges)
      log_p_A0 <- C1_edges * log(alpha_edges)    + C0_edges * log(1 - alpha_edges)
      
      samples_loglik_net[kk, ] <- A_edges * log_p_A1 + (1 - A_edges) * log_p_A0
      
      # ---------- (B) Continuous covariates loglik: p(X_cont | z, mu, tau2) ----------
      if (n_cont > 0) {
        mu_nodes  <- mu[z, , drop = FALSE]                       # n   n_cont
        sd_nodes  <- sqrt(pmax(tau2[z, , drop = FALSE], eps))     # n   n_cont
        
        samples_loglik_cont[kk, ] <- rowSums(
          dnorm(X_cont, mean = mu_nodes, sd = sd_nodes, log = TRUE)
        )
      }
      
      # ---------- (C) Categorical covariates loglik: p(X_cat | z, xi) ----------
      if (n_cat > 0) {
        ll_cat_vec <- rep(0, n)
        for (s in seq_len(n_cat)) {
          # xi[[s]]: K   F_s
          probs <- xi[[s]][cbind(z, X_cat[, s])]
          ll_cat_vec <- ll_cat_vec + log(pmax(probs, eps))
        }
        samples_loglik_cat[kk, ] <- ll_cat_vec
      }
    }
  }
  
  
  return(list(
    A     = samples_A,
    z     = samples_z,
    pi    = samples_pi,
    theta = samples_theta,
    psi = samples_psi,
    B     = samples_B,
    alpha = samples_alpha,
    beta  = samples_beta,
    mu    = if (n_cont > 0) samples_mu   else NULL,
    tau2  = if (n_cont > 0) samples_tau2 else NULL,
    xi    = if (n_cat > 0)  samples_xi   else NULL,
    burn_in = burn_in,
    loglik_net  = samples_loglik_net,
    loglik_cont = samples_loglik_cont,
    loglik_cat  = samples_loglik_cat
  ))
}


## ---- date generation, this can be used for different theta ----
dg <- function(
    seed,
    n,
    K,
    M,
    p_cont,
    p_cat,
    F_s,
    theta_true,
    B,          # K   K
    alpha,          # K   K
    beta,          # K   K
    mu,          # K   p_cont
    tau2,          # diag(p_cont)
    correct_rate   # length = p_cat
) {
  
  set.seed(seed)
  
  ## -------------------------
  ## 1. Community memberships
  ## -------------------------
  z_true  <- rep(1:K, each = n / K)
  pi_true <- rep(1 / K, K)
  
  ## -------------------------
  ## 2. Degree parameters  
  ## -------------------------
  #theta_true <- rnorm(n, mean = theta_mean, sd = theta_sd)
  theta_true <- theta_true
  ## -------------------------
  ## 3. Block matrix B
  ## -------------------------
  stopifnot(all(dim(B) == c(K, K)))
  B_true <- B
  
  ## -------------------------
  ## 4. Linear predictor  
  ## -------------------------
  B_mat_true <- B_true[z_true, z_true]
  psi_true   <- outer(theta_true, theta_true, "+") + B_mat_true
  diag(psi_true) <- 0
  
  ## -------------------------
  ## 5. True network A
  ## -------------------------
  A_true <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      prob <- plogis(psi_true[i, j])
      A_true[i, j] <- rbinom(1, 1, prob)
      A_true[j, i] <- A_true[i, j]
    }
  }
  
  ## -------------------------
  ## 6. Noise parameters
  ## -------------------------
  stopifnot(all(dim(alpha) == c(K, K)))
  alpha_true <- alpha
  
  
  stopifnot(all(dim(beta) == c(K, K)))
  beta_true <- beta
  
  A_tilde <- array(0, dim = c(n, n, M))
  
  for (m in 1:M) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- z_true[i]
        l <- z_true[j]
        
        obs <- if (A_true[i, j] == 1) {
          rbinom(1, 1, 1 - beta_true[k, l])
        } else {
          rbinom(1, 1, alpha_true[k, l])
        }
        
        A_tilde[i, j, m] <- obs
        A_tilde[j, i, m] <- obs
      }
    }
  }
  
  ## -------------------------
  ## 7. Continuous covariates
  ## -------------------------
  stopifnot(
    all(dim(mu)   == c(K, p_cont)),
    all(dim(tau2) == c(p_cont, p_cont))
  )
  
  X_cont <- matrix(0, n, p_cont)
  
  for (k in 1:K) {
    idx <- which(z_true == k)
    X_cont[idx, ] <- MASS::mvrnorm(
      n     = length(idx),
      mu    = mu[k, ],
      Sigma = tau2
    )
  }
  ## -------------------------
  ## 8. Categorical covariates
  ## -------------------------
  stopifnot(
    length(correct_rate) == p_cat,
    all(correct_rate > 0),
    all(correct_rate < 1)
  )
  
  X_cat   <- matrix(0, n, p_cat)
  xi_true <- array(0, dim = c(K, p_cat, F_s))
  
  ## construct xi_true
  for (j in 1:p_cat) {
    for (k in 1:K) {
      for (f in 1:F_s) {
        xi_true[k, j, f] <-
          ifelse(
            f == k,
            correct_rate[j],
            (1 - correct_rate[j]) / (F_s - 1)
          )
      }
    }
  }
  
  ## sample categorical covariates
  for (j in 1:p_cat) {
    for (k in 1:K) {
      idx <- which(z_true == k)
      X_cat[idx, j] <- sample(
        1:F_s,
        size = length(idx),
        replace = TRUE,
        prob = xi_true[k, j, ]
      )
    }
  }
  storage.mode(X_cat) <- "integer"
  
  ## -------------------------
  ## 9. Return
  ## -------------------------
  list(
    A_tilde = A_tilde,
    X_cont  = X_cont,
    X_cat   = X_cat,
    truth = list(
      z     = z_true,
      pi    = pi_true,
      theta = theta_true,
      B     = B_true,
      psi   = psi_true,
      A     = A_true,
      alpha = alpha_true,
      beta  = beta_true,
      mu    = mu,
      tau2  = diag(tau2),
      xi    = xi_true
    )
  )
}
