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
nt = 52 * 5
# Site mean case count
site_means = round(runif(n, 10, 100))
# length_scale: determines how quickly correlation decays with distance
#   - higher length_scale: correlation persists over longer distances (smoother spatial variation)
#   - lower length_scale: correlation decays rapidly, indicating localised variation
length_scale <- 10
# periodic_scale: strength of repeating (seasonal or cyclical) patterns
#   - higher periodic_scale: stronger seasonal patterns
#   - lower periodic_scale: weaker seasonal patterns
periodic_scale = 1
# long_term_scale: scale controlling decay rate of long-term temporal correlation
#   - higher long_term_scale: smoother long-term trends (correlation persists longer)
#   - lower long_term_scale: rapid loss of correlation, short-term variation dominates
long_term_scale = 300
# period: duration of the repeating cycle (e.g., 52 weeks for annual seasonality)
period = 52
# m_one: prob of outlier
m_one = 0.02
# m_switch: probability of switching between missing and observed states.
m_switch = 1
# Outlier factor
m_factor = 2
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

#make_outlier <- sample(1:nrow(true_data), 50)
#true_data$y[make_outlier] <- true_data$y[make_outlier] * 2

obs_data <- true_data |>
  observed_data(
    m_one = m_one,
    m_switch = m_switch,
    m_factor = m_factor
  )
#obs_data$outlier[make_outlier] <- 1

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
  geom_point(data = obs_data, aes(x = t, y = y_obs), size = 1, colour = "black") +
  geom_point(data = filter(obs_data, outlier == 1), aes(x = t, y = y_obs), size = 1, colour = "red") +
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

fit_data$surprisal[fit_data$surprisal < 0] <- 0
roc_obj <- pROC::roc(fit_data$outlier, fit_data$surprisal)
auc_value <- pROC::auc(roc_obj)
roc_pd <- as.data.frame(roc_obj[c("thresholds", "sensitivities", "specificities")])
roc_plot <- ggplot(data = roc_pd, aes(x = specificities, y = sensitivities)) +
  geom_abline(slope = 1, intercept = 1, linetype = 2, col = "grey50") +
  geom_line(col = "chartreuse3", linewidth = 1) +
  #geom_label(aes(x = 0.1, y = 0.1), label = auc_value) +
  xlim(1, 0) +
  xlab("Specificity") +
  ylab("Sensitivity") +
  coord_fixed() +
  theme_bw()

best_thr <- pROC::coords(roc_obj, x = "best", best.method = "youden", transpose = TRUE)[1]

obs_data$surprisal <- fit_data$surprisal

ggplot(data= obs_data, aes(x= factor(outlier), y = surprisal)) +
  geom_jitter(alpha = 0.3) +
  geom_boxplot(fill = NA)

sim_plot <- ggplot() +
  geom_line(data = true_data, aes(x = t, y = lambda)) +
  geom_point(data = obs_data, aes(x = t, y = y_obs, colour = surprisal, size = surprisal)) +
  geom_point(
    data = dplyr::filter(obs_data, outlier == 1),
    aes(x = t, y = y_obs),
    shape = 19, size = 1, colour = "chartreuse3"
  ) +
  scale_colour_viridis_c(option =  "C", direction = -1) +
  scale_size_continuous(range = c(0.1, 3)) +
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



