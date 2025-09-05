source("R/simulate.R")
library(weave)
library(progress)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

#' Case counts:
#'  y_st ~ Pois(λ_st)
#' Cases mean:
#'  λ_st = exp(z_st)
#' Cases mean (log scale):
#'  z_st = µ_s + f_st
#' Where we have a site specific intercept:
#'  µ_s
#' and a latent GP:
#'  f_st ~ MVN(0,Σ)
#' Spatiotemporal covariance:
#'  Σ = dist_k ⊗ time_k


# True parameters --------------------------------------------------------------
set.seed(321234)
# Number of sites
n = 25
# Number of timesteps
nt = 52 * 3
# Site mean case count
site_means = round(runif(n, 2, 100))
# length_scale: determines how quickly correlation decays with distance
#   - higher length_scale: correlation persists over longer distances (smoother spatial variation)
#   - lower length_scale: correlation decays rapidly, indicating localised variation
length_scale <- 2
# periodic_scale: strength of repeating (seasonal or cyclical) patterns
#   - higher periodic_scale: stronger seasonal patterns
#   - lower periodic_scale: weaker seasonal patterns
periodic_scale = 1
# long_term_scale: scale controlling decay rate of long-term temporal correlation
#   - higher long_term_scale: smoother long-term trends (correlation persists longer)
#   - lower long_term_scale: rapid loss of correlation, short-term variation dominates
long_term_scale = 100
# period: duration of the repeating cycle (e.g., 52 weeks for annual seasonality)
period = 52
# p_one: probability that a new cluster begins with a missing value.
#   - higher p_one: more missing data overall.
#   - lower p_one: less missing data overall.
p_one = 0.25
# p_switch: probability of switching between missing and observed states.
#   - lower p_switch: longer sequences (clusters) of missingness or non-missingness.
#   - higher p_switch: shorter clusters, more frequent switching between states.
p_switch = 0.2
# ------------------------------------------------------------------------------

# Simulated data ---------------------------------------------------------------
coordinates <- data.frame(
  id = factor(1:n),
  lat = 1:n,
  lon = 1:n,
  mu = log(site_means)
)

space_k <- space_kernel(
  coordinates = coordinates,
  length_scale = length_scale
)

time_k <- time_kernel(
  times = 1:nt,
  periodic_scale = periodic_scale,
  long_term_scale = long_term_scale,
  period = period
)

true_data <- simulate_data(
  n = n,
  nt = nt,
  coordinates = coordinates,
  space_k = space_k,
  time_k = time_k
)

obs_data <- observed_data(
  data = true_data,
  p_one = p_one,
  p_switch = p_switch
)

space_pd <- data.frame(
  distance = 1:100
) |>
  mutate(
    k = rbf_kernel(distance, theta = length_scale),
    group = "True"
  )

space_plot <- ggplot(space_pd, aes(x = distance, y = k, colour = group)) +
  geom_line() +
  theme_bw()

time_pd <- data.frame(
  week = 1:(nt)
) |>
  mutate(
    k = periodic_kernel(x = week, alpha = periodic_scale, period = period) *
      rbf_kernel(x = week, theta = long_term_scale),
    group = "True"
  )

time_plot <- ggplot(time_pd, aes(x = week, y = k, colour = group)) +
  geom_line() +
  theme_bw()

hf_labeller <- function(value) {
  paste("HF:", value)
}

sim_plot <- ggplot() +
  geom_point(data = true_data, aes(x = t, y = y), size = 1, colour = "red") +
  geom_point(data = obs_data, aes(x = t, y = y_obs), size = 1, colour = "black") +
  geom_line(data = true_data, aes(x = t, y = lambda)) +
  facet_wrap(~ id, scales = "free_y", labeller = labeller(id = hf_labeller)) +
  ylab("Cases") +
  xlab("Week") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "grey50"),
    strip.text = element_text(size = 8, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  ggtitle("Simulated data")
# ------------------------------------------------------------------------------

# Infer kernel hyper-parameters ------------------------------------------------
infer_space <- infer_space_kernel_params(obs_data, nt = nt, n = n, TRUE)
infer_time <- infer_time_kernel_params(obs_data, 52, nt = nt, n = n, plot = TRUE)

hyperparameters <- c(infer_space$length_scale, infer_time$periodic_scale, infer_time$long_term_scale)
#hyperparameters <- c(3, 1, 200)
## Options
# Fix these conservatively - simulated sensitvity analyses?
# Fit Bayesian with priors
# ------------------------------------------------------------------------------

# Fit --------------------------------------------------------------------------
system.time({
  fit_data <- fit(obs_data, coordinates, hyperparameters, n, nt)
})


fit_plot <- sim_plot +
  geom_ribbon(
    data = fit_data,
    aes(x = t, ymin = pred_Q2.5, ymax = pred_Q97.5, fill = id), alpha = 0.25
  ) +
  geom_ribbon(
    data = fit_data,
    aes(x = t, ymin = pred_Q25, ymax = pred_Q75, fill = id, alpha = 0.5)
  ) +
  geom_line(
    data = fit_data,
    aes(x = t, y = data_Q50, col = id), linewidth = 1
  ) +
  ggtitle("Filling missingness")
fit_plot

compare_pd <- data.frame(
  Modelled = fit_data$data_Q50,
  True = true_data$y,
  type = factor(ifelse(is.na(fit_data$y_obs), "Missing", "Observed"), levels = c("Observed", "Missing"))
)

comp_plot <- ggplot() +
  geom_point(data = compare_pd, aes(x = True, y = Modelled, colour = type, alpha = type), shape = 19) +
  geom_abline(intercept = 0, slope = 1) +
  scale_colour_manual(values = c("chartreuse3", "darkmagenta")) +
  scale_alpha_manual(values = c(0.25, 1)) +
  theme_bw()


