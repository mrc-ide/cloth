# TODO:
## 1. Pull out all pre-computable elements that wouldn't change with use case.
## 2. Factor in proper fitting regime from weave

test_fit <- function(obs_data, coordinates, hyperparameters, n, nt, period = 52, draw = FALSE){

  ## --- Kernels (separable; no full K ever built) ---------------------------
  time_mat  <- time_kernel(
    times = 1:nt,
    periodic_scale = hyperparameters[2],
    long_term_scale = hyperparameters[3],
    period = period
  )
  space_mat <- space_kernel(
    coordinates = coordinates, length_scale = hyperparameters[1]
  )

  obs_idx <- which(!is.na(obs_data$y_obs))
  N <- n * nt

  ## --- Working response on log scale + heteroscedastic nugget -------------
  # f_infer = log(y + 1) - mu_infer  (your design)
  y_work <- obs_data$f_infer[obs_idx]
  lam_hat <- exp(obs_data$mu_infer[obs_idx])
  # Poisson+log delta-method variance on the observed log scale
  noise_var <- lam_hat / (lam_hat + 1)^2

  # diag(K) for Jacobi preconditioner (diag(space) ⊗ diag(time))
  kdiag_full <- as.vector(kronecker(diag(space_mat), diag(time_mat)))

  ## --- Posterior mean via PCG on observed system ---------------------------
  # (S K S^T + D) alpha = y_work
  alpha <- weave:::pcg(y_work, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-6)
  # f_hat = K S^T alpha
  f_hat <- weave:::kron_mv(weave:::fill_vector(alpha, obs_idx, N), space_mat, time_mat)
  z_hat <- f_hat + obs_data$mu_infer
  out <- exp(z_hat)

  if(draw){
    m <- length(obs_idx)
    # Prior GP draw: eta ~ N(0, K) using separable sampling
    eta <- weave:::quick_mvnorm(space_mat, time_mat)                                   # uses chol(space), chol(time) only
    # Pseudo-noise on observed log scale
    eps <- stats::rnorm(m, sd = sqrt(noise_var))
    # Solve: (S K S^T + D) alpha_tilde = (S eta + eps) - y_work
    rhs <- (eta[obs_idx] + eps) - y_work
    alpha_tilde <- weave:::pcg(rhs, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = 1e-6)
    # Posterior GP draw for f (on centred log scale)
    f_draw <- eta - weave:::kron_mv(weave:::fill_vector(alpha_tilde, obs_idx, N), space_mat, time_mat)
    # Convert to z = f + mu, then lambda = exp(z)
    z_draw <- f_draw + obs_data$mu_infer
    out <- exp(z_draw)
  }

  return(out)
}

# Hacky fit example
fit_data <- obs_data
fit_data$y_obs[sample(1:nrow(fit_data), 0.2 * nrow(fit_data))] <- NA

opt <- function(par){
  pred <- test_fit(fit_data, coordinates, par, n, nt)
  sum(dpois(obs_data$y_obs, pred, log = TRUE), na.rm = TRUE)
}

fit_opt <- optim(
  par = c(1, 5, 100),
  fn = opt,
  method = "L-BFGS-B",
  lower = c(0.0001, 0.8, 52 * 1.5),
  upper = c(2, 10, 500),
  control = list(trace = 1, maxit = 100, fnscale = -1)
)
fit_opt$par
hyperparameters <- fit_opt$par
#hyperparameters <- c(length_scale, periodic_scale, long_term_scale)


# Outputs
obs_data$posterior_mean <- test_fit(obs_data, coordinates, hyperparameters, n, nt)

draws <- lapply(1:50, function(x){
  print(x)
  data.frame(draw = x, t = obs_data$t, id = obs_data$id, r = test_fit(obs_data, coordinates, hyperparameters, n, nt, draw = TRUE))
}) |>
  bind_rows()


# approximation
#obs_data$l <- test_fit(obs_data, coordinates, hyperparameters, n, nt, q = 0.025)
#obs_data$h<- test_fit(obs_data, coordinates, hyperparameters, n, nt, q = 0.975)

bounds <- draws |>
  summarise(
    l = quantile(r, 0.025),
    l2 = quantile(r, 0.25),
    u = quantile(r, 0.975),
    u2 = quantile(r, 0.75),
    .by = c("id", "t")
  )

ggplot() +
  #geom_line(data = draws, aes(x = t, y = r, group = draw), col = "orange") +
  geom_ribbon(data = bounds, aes(x = t, ymin = l, ymax = u, fill = id), alpha = 0.5) +
  geom_ribbon(data = bounds, aes(x = t, ymin = l2, ymax = u2, fill = id), alpha = 0.75) +
  geom_line(data = obs_data, aes(x = t, y = posterior_mean), col = "deeppink") +
  #geom_point(data = true_data, aes(x = t, y = y), size = 0.2, colour = "red") +
  geom_point(data = obs_data, aes(x = t, y = y_obs), size = 0.4, colour = "black") +
  geom_line(data = obs_data, aes(x = t, y = l), col = "chartreuse") +
  geom_line(data = obs_data, aes(x = t, y = h), col = "chartreuse") +
  facet_wrap( ~ id, scales = "free_y") +
  theme_bw()# + scale_y_log10()


# 1) Low-rank factors for space and time (use your favourite: pivoted chol or eig trunc)
lr_space <- function(space_mat, r_s) { # returns n x r_s
  eig <- eigen(space_mat, symmetric = TRUE)
  idx <- seq_len(r_s)
  eig$vectors[, idx, drop=FALSE] %*% diag(sqrt(pmax(eig$values[idx], 0)), r_s)
}
lr_time <- function(time_mat, r_t) { # returns nt x r_t
  eig <- eigen(time_mat, symmetric = TRUE)
  idx <- seq_len(r_t)
  eig$vectors[, idx, drop=FALSE] %*% diag(sqrt(pmax(eig$values[idx], 0)), r_t)
}

# 2) Build Kronecker low rank L = Ls ⊗ Lt  (N x r)
build_L <- function(Ls, Lt) {
  # Khatri–Rao-friendly: L = kron(Ls, Lt)
  kronecker(Ls, Lt)  # N x (r_s * r_t)
}

gp_var_diag_woodbury <- function(obs_idx, N, space_mat, time_mat, noise_var,
                                 r_s = 16L, r_t = 16L) {
  n  <- nrow(space_mat); nt <- nrow(time_mat)
  stopifnot(n * nt == N)

  Ls <- lr_space(space_mat, r_s)    # n x r_s
  Lt <- lr_time(time_mat,  r_t)     # nt x r_t
  L  <- build_L(Ls, Lt)             # N x r (r = r_s * r_t)

  # Observed rows and weights
  L_S   <- L[obs_idx, , drop = FALSE]                    # m x r
  Dinv  <- 1 / noise_var                                 # length m
  # G = L_S^T D^{-1} L_S  (r x r) via weighted crossprod
  G     <- crossprod(L_S * Dinv, L_S)                    # r x r
  # Solve (I + G) once; compute M = G (I+G)^{-1}
  IR    <- diag(ncol(G))
  C     <- chol(IR + G)                                  # r x r
  M     <- backsolve(C, forwardsolve(t(C), G))           # G %*% (I+G)^{-1}

  # diag(K) under low-rank: rowSums(L^2)
  kdiag <- rowSums(L * L)                                # length N

  # c_i are the rows of L (as column vectors); compute c_i^T M c_i
  # Efficiently: diag(L %*% M %*% t(L)) = rowSums((L %*% M) * L)
  LM    <- L %*% M                                       # N x r
  qform <- rowSums(LM * L)                               # length N

  pmax(kdiag - qform, 0)
}
