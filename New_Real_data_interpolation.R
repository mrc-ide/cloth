source("R/simulate.R")
library(weave)
library(progress)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(sf)
library(lubridate)

period <- 12

sle_raw <- readRDS("data/SLE_routine_HF_malaria_outputs_with_chiefdom_centroids.rds")

set.seed(12345)
n_hf <- 12

sle <- sle_raw |>
  filter(variable == "Malaria treated with ACT <24 hours 15+y_X") |>
  mutate(
    # Convert "month year" to a proper Date (first day of month)
    ym = my(paste(month, year)),
    # Find the first month across dataset
    m = interval(min(ym), ym) %/% months(1) + 1,
    lon = sapply(centroid, function(x){ x[2] }),
    lat = sapply(centroid, function(x){ x[1] })
  ) |>
  mutate(t = m, n = cases) |>
  select(name_2, name_3, HF, t, n, lat, lon) |>
  sf::st_drop_geometry() |>
  dplyr::mutate(
    lat = lat + runif(1, 0, 0.5),
    lon = lon + runif(1, 0, 0.5),
    .by = c("name_2", "name_3", "HF")
  ) |>
  dplyr::distinct(name_2, name_3, HF, t, .keep_all = TRUE)

obs_data <- weave:::data_process(sle, name_2, name_3, HF) |>
  rename(y_obs = n) |>
  dplyr::mutate(
    z_infer = log1p(.data$y_obs),
    mu_infer = mean(.data$z_infer, na.rm = TRUE),
    f_infer = .data$z_infer - .data$mu_infer,
    .by = id
  ) |>
  #filter(id %in% sample(unique(id), n_hf, replace = FALSE)) |>
  filter(name_2 %in% c("BOMBALI", "BO")) |>
  mutate(name = paste0(name_2, "|", name_3, "|",  HF))

mean(is.na(obs_data$y_obs))
length(unique(obs_data$id))

hf_labeller <- function(value) {
  paste("HF:", value)
}

sim_plot <- ggplot() +
  geom_point(data = obs_data, aes(x = t, y = y_obs), size = 1, colour = "black") +
  facet_wrap(~ name, scales = "free_y", labeller = labeller(id = hf_labeller)) +
  ylab("Cases") +
  xlab("Month") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", colour = "grey50"),
    strip.text = element_text(size = 8, face = "bold"),
    panel.spacing = unit(0.5, "lines")
  ) +
  ggtitle("HF data")

n <- length(levels(obs_data$id))
nt <- max(obs_data$t)
coordinates <- unique(obs_data[, c("id", "lat", "lon")])
estimate <- fit(obs_data, nt, period, n_sites = 10, mask_prop = 0.2,
                par0 = c(space = 0.00015, t_per = 0.5, t_long = 22),
                lower = c(space =  0.0001, t_per =  0.1, t_long = 12),
                upper = c(space = 0.001, t_per = 2, t_long = 3 * period))
hyperparameters <- estimate$par
#hyperparameters <- c(0.00001, 2, 18)


state <- gp_build_state(obs_data, coordinates, hyperparameters, n, nt, period)

# Posterior mean:
obs_data$posterior_mean <- gp_posterior_mean(state)

interval <- bounds(
  state,
  n_lambda = 25,
  n_draw = 100,
  quantiles = c(0.025, 0.25, 0.75, 0.975)
)

obs_data <- obs_data |>
  left_join(interval)

obs_data_sub <- obs_data |>
  filter(as.numeric(id) < 40)

fit_plot <- ggplot() +
  geom_vline(xintercept = 12 * 1:7, linetype = 2, col = "grey40") +
  geom_vline(xintercept =  77, col = "chartreuse3") +
  geom_vline(xintercept =  65, col = "dodgerblue") +
  geom_ribbon(data = obs_data_sub, aes(x = t, ymin = q0.025, ymax = q0.975, fill = id), alpha = 0.5) +
  geom_ribbon(data = obs_data_sub, aes(x = t, ymin = q0.25, ymax = q0.75, fill = id), alpha = 0.75) +
  geom_line(data = obs_data_sub, aes(x = t, y = posterior_mean), col = "deeppink") +
  geom_point(data = obs_data_sub, aes(x = t, y = y_obs), size = 0.4, colour = "black") +
  facet_wrap( ~ id, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none")
