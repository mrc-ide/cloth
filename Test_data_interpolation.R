source("R/simulate.R")
library(weave)
library(progress)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(sf)
library(lubridate)


sle_raw <- readRDS("data/SLE_routine_HF_malaria_outputs_with_chiefdom_centroids.rds")

set.seed(12345)
n_hf <- 36 * 2

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
  filter(id %in% sample(unique(id), n_hf, replace = FALSE)) |>
  #filter(name_2 %in% c("BO")) |>
  mutate(name = paste0(name_2, "|", HF))

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
# Infer kernel hyper-parameters ------------------------------------------------
n <- length(levels(obs_data$id))
nt <- max(obs_data$t)
infer_space <- infer_space_kernel_params(obs_data, nt = nt, n = n, TRUE)
infer_time <- infer_time_kernel_params(obs_data, 12, nt = nt, n = n, plot = TRUE)
hyperparameters <- c(0.0000001, infer_time$periodic_scale / 2, infer_time$long_term_scale)
# ------------------------------------------------------------------------------

# Fit --------------------------------------------------------------------------
coordinates <- unique(obs_data[, c("id", "lat", "lon")])
system.time({
  fit_data <- fit(obs_data, coordinates, hyperparameters, n, nt)
})


fit_plot <- sim_plot +
  geom_vline(xintercept = 12 * 1:7, linetype = 2, col = "grey40") +
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
# ------------------------------------------------------------------------------
