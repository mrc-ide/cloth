# helper: Kronecker diag without allocating a big dense kronecker()
.kdiag_from_factors <- function(space_diag, time_diag, n, nt) {
  # time varies fastest (sites × times)
  rep(space_diag, each = nt) * rep(time_diag, times = n)
}

# GP prior draw matching kron_mv() layout:
#   - X is (sites × times)
#   - return vec(t( Rs %*% X %*% t(Rt) ))
.quick_mvnorm_from_chol <- function(Rs, Rt, n, nt) {
  X <- matrix(stats::rnorm(n * nt), nrow = n, ncol = nt)  # white noise (sites × times)
  Y <- t(Rs) %*% X %*% Rt                                 # use t(Rs) on the left, Rt on the right
  as.vector(t(Y))                                          # vec(t(Y)): time-fast vector
}

# 1) Build a reusable sampler state -----------------------------------------
gp_build_state <- function(obs_data, coordinates, hyperparameters, n, nt, period = 52) {
  time_mat  <- time_kernel(times = 1:nt,
                           periodic_scale = hyperparameters[2],
                           long_term_scale = hyperparameters[3],
                           period = period)
  space_mat <- space_kernel(coordinates = coordinates,
                            length_scale = hyperparameters[1])

  obs_idx   <- which(!is.na(obs_data$y_obs))
  N         <- n * nt

  # Working response (log scale), heteroscedastic nugget
  y_work    <- obs_data$f_infer[obs_idx]
  lam_hat   <- exp(obs_data$mu_infer[obs_idx])
  noise_var <- lam_hat / (lam_hat + 1)^2

  # Preconditioner diag(K) = diag(space) ⊗ diag(time)
  kdiag_full <- .kdiag_from_factors(diag(space_mat), diag(time_mat), n, nt)

  # Optional: precompute Choleskys for fast prior draws
  Rt <- chol(time_mat)   # nt × nt
  Rs <- chol(space_mat)  #  n ×  n

  # A-solve closure for (S K S^T + D) y = rhs
  A_solve <- function(rhs_obs, tol = 1e-6) {
    weave:::pcg(rhs_obs, obs_idx, N, space_mat, time_mat, noise_var, kdiag_full, tol = tol)
  }

  list(
    time_mat = time_mat,
    space_mat = space_mat,
    Rt = Rt, Rs = Rs,
    obs_idx = obs_idx, N = N,
    y_work = y_work, noise_var = noise_var,
    kdiag_full = kdiag_full,
    A_solve = A_solve,
    mu_infer = obs_data$mu_infer,     # full-length vector for z = f + mu
    n = n, nt = nt
  )
}

# 2) Posterior mean once -----------------------------------------------------
gp_posterior_mean <- function(state, tol = 1e-6) {
  alpha <- state$A_solve(state$y_work, tol = tol)             # A^{-1} y
  f_hat <- weave:::kron_mv(
    weave:::fill_vector(alpha, state$obs_idx, state$N),
    state$space_mat, state$time_mat
  )
  z_hat <- f_hat + state$mu_infer
  exp(z_hat)  # lambda_hat
}

# 3) One posterior draw, reusing everything ---------------------------------
gp_draw <- function(state, tol = 1e-6) {
  m   <- length(state$obs_idx)

  # Prior draw using precomputed Choleskys (no fresh chol each time)
  #eta <- .quick_mvnorm_from_chol(state$Rs, state$Rt, state$n, state$nt)       # N-vector
  eta <- weave:::quick_mvnorm(state$space_mat, state$time_mat)

  # Pseudo-noise on observed log scale (heteroscedastic)
  eps <- stats::rnorm(m, sd = sqrt(state$noise_var))

  # Solve A alpha_tilde = (S eta + eps) - y_work
  rhs <- (eta[state$obs_idx] + eps) - state$y_work
  alpha_tilde <- state$A_solve(rhs, tol = tol)

  # Posterior draw for f (centred log scale): eta - K S^T alpha_tilde
  f_draw <- eta - weave:::kron_mv(
    weave:::fill_vector(alpha_tilde, state$obs_idx, state$N),
    state$space_mat, state$time_mat
  )

  z_draw <- f_draw + state$mu_infer
  exp(z_draw)
}

state <- gp_build_state(obs_data, coordinates, hyperparameters, n, nt)

t1 <- bench::mark(
  .quick_mvnorm_from_chol(state$Rs, state$Rt, state$n, state$nt),
  weave:::quick_mvnorm(state$space_mat, state$time_mat), check = FALSE
)
plot(t1)


# mean once
obs_data$posterior_mean <- gp_posterior_mean(state)

draws <- lapply(1:30, function(x){
  print(x)
  data.frame(draw = x, t = obs_data$t, id = obs_data$id, r =  gp_draw(state))
}) |>
  bind_rows() |>
  mutate(rp = rpois(n(), r))

bounds <- draws |>
  summarise(
    l = quantile(rp, 0.025),
    l2 = quantile(rp, 0.25),
    u = quantile(rp, 0.975),
    u2 = quantile(rp, 0.75),
    .by = c("id", "t")
  )


B <- 30L  # bump draws a bit for stable quantiles
lam_draws <- replicate(B, gp_draw(state), simplify = "array")

# Add observation noise:
count_draws <- lapply(1:100, function(x){
  count_draws <- lam_draws
  count_draws[] <- rpois(length(lam_draws), lambda = lam_draws)
  return(count_draws)
})
count_draws <-  do.call("cbind", count_draws)

bounds <- data.frame(
  id = obs_data$id,
  t = obs_data$t,
  l  = apply(count_draws, 1, quantile, 0.025),
  l2  = apply(count_draws, 1, quantile, 0.25),
  u  = apply(count_draws, 1, quantile, 0.975),
  u2  = apply(count_draws, 1, quantile, 0.75)
)

ggplot() +
  #geom_line(data = draws, aes(x = t, y = r, group = draw), col = "orange") +
  geom_ribbon(data = bounds, aes(x = t, ymin = l, ymax = u, fill = id), alpha = 0.5) +
  geom_ribbon(data = bounds, aes(x = t, ymin = l2, ymax = u2, fill = id), alpha = 0.75) +
  geom_line(data = obs_data, aes(x = t, y = posterior_mean), col = "deeppink") +
  geom_point(data = true_data, aes(x = t, y = y), size = 0.1, colour = "red") +
  geom_point(data = obs_data, aes(x = t, y = y_obs), size = 0.4, colour = "black") +
  facet_wrap( ~ id, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none")
