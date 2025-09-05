#' Generate clustered binary sequence
#'
#' Creates a binary vector of length `n` where consecutive values tend to
#' cluster. Values are freshly drawn using `p_one` when a new cluster starts and
#' otherwise repeat the previous value.
#'
#' @param n Number of draws to generate.
#' @param p_one Probability that a new cluster begins with a one. Also used for
#'   the first draw.
#' @param p_switch Probability of switching to a new cluster at each step.
#'
#' @return Numeric vector of 0s and 1s.
#'
#' @details Randomness is generated using [stats::rbinom()] and
#'   [stats::runif()]. Set a seed via [set.seed()] for reproducible results.
#'
#' @examples
#' set.seed(1)
#' generate_clustered_binary(5, 0.2, 0.1)
#'
#' @export
generate_clustered_binary <- function(n, p_one, p_switch) {
  result <- numeric(n)
  result[1] <- rbinom(1, 1, p_one)
  for (i in 2:n) {
    if (runif(1) < p_switch) {
      result[i] <- rbinom(1, 1, p_one)
    } else {
      result[i] <- result[i - 1]
    }
  }
  return(result)
}

#' Simulate spatiotemporal case counts
#'
#' This function generates synthetic case count data under a spatiotemporal
#' Poisson model with site-specific effects and a latent Gaussian process.
#'
#' The model is defined as:
#' \deqn{y_{st} \sim \text{Pois}(\lambda_{st})}
#' \deqn{\lambda_{st} = \exp(z_{st})}
#' \deqn{z_{st} = \mu_s + f_{st}}
#'
#' where:
#' - \eqn{\mu_s} is a site-specific intercept,
#' - \eqn{f_{st}} is a latent Gaussian process with mean zero and covariance
#'   \eqn{\Sigma},
#' - \eqn{\Sigma} has separable spatiotemporal structure given by
#'   \eqn{\Sigma = \text{dist}_k \otimes \text{time}_k}, the Kronecker product
#'   of a spatial kernel and a temporal kernel.
#'
#' This formulation allows counts at each site and time point to vary according
#' to both site-level heterogeneity and correlated spatiotemporal fluctuations.
#'
#' @param n Number of spatial locations.
#' @param nt Number of time points.
#' @param coordinates Data frame with columns `id`, `lat` and `lon` describing
#'   site locations.
#' @param space_k Spatial covariance matrix.
#' @param time_k Temporal covariance matrix.
#'
#' @return Tibble containing the latent effect `f`, mean `z`, rate `lambda` and
#'   simulated counts `y` for each site-time combination.
#'
#' @details Randomness comes from multivariate normal and Poisson draws. Use
#'   [set.seed()] before calling for reproducible simulations.
#'
#' @export
simulate_data <- function(
    n, nt,
    coordinates,
    space_k, time_k) {

  output_df <- tidyr::expand_grid(
    id = factor(1:n),
    t = 1:nt
  ) |>
    dplyr::left_join(
      coordinates, by = "id"
    ) |>
    dplyr::mutate(
      f = quick_mvnorm(space_k, time_k),
      z = .data$mu + .data$f,
      lambda = exp(.data$z),
      y = rpois(n * nt, .data$lambda)
    )

  return(output_df)
}

#' Observed data generation
#'
#' @param data Output of [simulate_data()].
#' @param p_one Probability that a new missingness cluster starts with a missing
#'   value.
#' @param p_switch Probability of switching between missing and observed
#'   clusters.
#' @param end Add NAs to end of sequence, for testing forecasting
#'
#' @return Tibble with observed counts `y_obs` and inferred `mu_infer`,
#'   `z_infer` and `f_infer`.
#'
#' @details Missingness is generated stochastically via
#'   [generate_clustered_binary()]. Set a seed with [set.seed()] to reproduce
#'   results.
#'
#' @export
observed_data <- function(data, p_one = 0, p_switch = 0, end = 0, m_one = 0, m_switch = 0, m_factor = 1) {
  data |>
    dplyr::mutate(
      y_obs = .data$y,
      # Add missingness
      missing = generate_clustered_binary(dplyr::n(), p_one, p_switch),
      y_obs = ifelse(missing == 1, NA, .data$y_obs),
      # Add NAs at end (for forecasting)
      y_obs = ifelse(t > (max(t) - end), NA, .data$y_obs),
      # Add outliers
      outlier = generate_clustered_binary(dplyr::n(), m_one, m_switch),
      y_obs = round(ifelse(outlier == 1, y_obs * sample(c(m_factor, 1 / m_factor), sum(outlier), replace = TRUE), y_obs))
    ) |>
    dplyr::mutate(
      z_infer = log1p(.data$y_obs),
      mu_infer = mean(.data$z_infer, na.rm = TRUE),
      f_infer = .data$z_infer - .data$mu_infer,
      .by = id
    ) |>
    dplyr::select(dplyr::any_of(c("id", "t", "lat", "lon", "y_obs", "mu_infer", "z_infer", "f_infer", "outlier")))
}
