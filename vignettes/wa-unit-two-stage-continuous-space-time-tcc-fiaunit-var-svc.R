library(stLMM)
library(RhpcBLASctl)
library(sf)
library(tidyverse)
library(knitr)
library(kableExtra)

source("article-utils.R")

set.seed(1)
theme_set(theme_bw(base_size = 12))

article_id <- "wa-unit-two-stage-continuous-space-time-tcc-fiaunit-var-svc"
data_dir <- "wa_data"

wa_counties <- read_rds(file.path(data_dir, "wa_counties.rds"))
unit_plots <- read_csv(file.path(data_dir, "wa_unit_plots.csv"), show_col_types = FALSE)
direct_estimates <- read_csv(file.path(data_dir, "wa_direct_estimates.csv"), show_col_types = FALSE)
prediction_grid <- read_csv(
  file.path(data_dir, "wa_unit_prediction_grid_5km.csv"),
  show_col_types = FALSE
)

wa_counties <- wa_counties |>
  mutate(
    county_fips = as.character(county_fips),
    fia_unit = factor(fia_unit)
  )

unit_plots <- unit_plots |>
  mutate(
    county_fips = as.character(county_fips),
    county_fips = factor(county_fips),
    fia_unit = factor(fia_unit),
    year = as.integer(year),
    time = as.numeric(time),
    x_km = x / 1000,
    y_km = y / 1000,
    biomass_positive = as.integer(agb_live_mg_ha > 0)
  ) |>
  arrange(county_fips, year, plot_id)

direct_estimates <- direct_estimates |>
  mutate(
    county_fips = as.character(county_fips),
    year = as.integer(year)
  )

plot_years <- sort(unique(unit_plots$year))
tcc_center <- mean(unit_plots$tcc_mean, na.rm = TRUE)
tcc_scale <- sd(unit_plots$tcc_mean, na.rm = TRUE)

unit_plots <- unit_plots |>
  mutate(tcc_mean_scaled = (tcc_mean - tcc_center) / tcc_scale)

county_levels <- levels(unit_plots$county_fips)
fia_unit_levels <- levels(unit_plots$fia_unit)
fia_unit_names <- unit_plots |>
  distinct(fia_unit, fia_unit_name) |>
  arrange(fia_unit)

prediction_grid <- prediction_grid |>
  mutate(
    county_fips = as.character(county_fips),
    county_fips = factor(county_fips, levels = county_levels),
    fia_unit = factor(fia_unit, levels = fia_unit_levels),
    year = as.integer(year),
    time = as.numeric(year),
    tcc_mean_scaled = (tcc_mean - tcc_center) / tcc_scale
  ) |>
  filter(year %in% plot_years) |>
  arrange(county_fips, year, x, y)

positive_plots <- unit_plots |>
  filter(biomass_positive == 1) |>
  mutate(biomass_sqrt = sqrt(agb_live_mg_ha))

data_summary <- tibble(
  quantity = c(
    "plot rows",
    "positive-biomass plot rows",
    "FIA units",
    "5 km prediction-grid rows",
    "measurement years"
  ),
  value = c(
    fmt_int(nrow(unit_plots)),
    fmt_int(nrow(positive_plots)),
    fmt_int(n_distinct(unit_plots$fia_unit)),
    fmt_int(nrow(prediction_grid)),
    paste(range(plot_years), collapse = "-")
  )
)

show_table(data_summary, caption = "Unit-response data for the continuous space-time model with spatial SVC.")

ggplot() +
  geom_sf(data = wa_counties, fill = "grey96", color = "grey70", linewidth = 0.2) +
  geom_point(
    data = unit_plots,
    aes(x, y, color = time),
    size = 0.35,
    alpha = 0.55
  ) +
  coord_sf(expand = FALSE) +
  scale_color_gradientn(colors = stlmm_palette(), name = "measurement time") +
  labs(x = NULL, y = NULL)

cor_at_range <- 0.05

decay_from_range <- function(range_05) {
  -log(cor_at_range) / range_05
}

## Broad process: approximately 25-100 km and 1-500 years.
phi_1_bounds <- decay_from_range(rev(c(25, 100)))
lambda_1_bounds <- decay_from_range(rev(c(1, 500)))

## Local process: approximately 1-25 km and 1-500 years.
phi_2_bounds <- decay_from_range(rev(c(1, 25)))
lambda_2_bounds <- decay_from_range(rev(c(1, 500)))

## Spatial SVC process: approximately 10-150 km.
phi_svc_bounds <- decay_from_range(rev(c(10, 150)))

n_samples <- 1000
burnin <- 1
chains <- 3
n_keep <- 25
posterior_thin <- max(1L, floor((n_samples - burnin) / n_keep))
posterior_sub_sample <- list(start = burnin + 1, thin = posterior_thin)
n_report <- 500

## Let stLMM use OpenMP threads while preventing BLAS from nesting threads
## inside the sampler and prediction calculations.
blas_set_num_threads(1)
n_omp_threads <- as.integer(Sys.getenv("STLMM_OMP_THREADS", unset = "25"))

fit_st_scale <- 2
pred_st_scale <- fit_st_scale

sqrt_response_sd <- sd(positive_plots$biomass_sqrt, na.rm = TRUE)



jitter_start <- function(value, amount = 0.01, bounds = NULL) {
  out <- rep(value, chains) * exp(rnorm(chains, 0, amount))

  if (!is.null(bounds)) {
    out <- pmin(pmax(out, bounds[1] + 1e-8), bounds[2] - 1e-8)
  }

  out
}


occurrence_start <- list(
  beta = c(
    `(Intercept)` = -0.18,
    tcc_mean_scaled = 10
  ),
  nngp_1 = list(
    sigma_sq = 80,
    alpha = 0.52,
    phi_1 = 0.031,
    lambda_1 = 0.039,
    phi_2 = 1.7,
    lambda_2 = 1.6
  )
)

biomass_start <- list(
  resid = list(
    tau_sq = c(
      `5` = 28,
      `6` = 25,
      `7` = 26,
      `8` = 12,
      `9` = 11
    )
  ),
  nngp_1 = list(
    sigma_sq = 2.3,
    alpha = 0.74,
    phi_1 = 0.037,
    lambda_1 = 0.059,
    phi_2 = 1.6,
    lambda_2 = 1.5
  ),
  nngp_2 = list(
    sigma_sq = 2.4,
    phi = 0.027
  )
)

occurrence_starting <- list(
  beta = occurrence_start$beta,
  nngp_1 = list(
    sigma_sq = jitter_start(occurrence_start$nngp_1$sigma_sq),
    alpha = jitter_start(occurrence_start$nngp_1$alpha, bounds = c(0, 1)),
    phi_1 = jitter_start(occurrence_start$nngp_1$phi_1, bounds = phi_1_bounds),
    lambda_1 = jitter_start(occurrence_start$nngp_1$lambda_1, bounds = lambda_1_bounds),
    phi_2 = jitter_start(occurrence_start$nngp_1$phi_2, bounds = phi_2_bounds),
    lambda_2 = jitter_start(occurrence_start$nngp_1$lambda_2, bounds = lambda_2_bounds)
  )
)

biomass_starting <- list(
  resid = biomass_start$resid,
  nngp_1 = list(
    sigma_sq = jitter_start(biomass_start$nngp_1$sigma_sq),
    alpha = jitter_start(biomass_start$nngp_1$alpha, bounds = c(0, 1)),
    phi_1 = jitter_start(biomass_start$nngp_1$phi_1, bounds = phi_1_bounds),
    lambda_1 = jitter_start(biomass_start$nngp_1$lambda_1, bounds = lambda_1_bounds),
    phi_2 = jitter_start(biomass_start$nngp_1$phi_2, bounds = phi_2_bounds),
    lambda_2 = jitter_start(biomass_start$nngp_1$lambda_2, bounds = lambda_2_bounds)
  ),
  nngp_2 = list(
    sigma_sq = jitter_start(biomass_start$nngp_2$sigma_sq),
    phi = jitter_start(biomass_start$nngp_2$phi, bounds = phi_svc_bounds)
  )
)

## occurrence_fit <- stLMM(
##   biomass_positive ~
##     tcc_mean_scaled +
##     nngp(x_km, y_km, time,
##          m = 15,
##          cov_model = "multi_res_sep_exp",
##          ordering = "maxmin",
##          st_scale = fit_st_scale),
##   data = unit_plots,
##   family = "binomial",
##   starting = occurrence_starting,
##   priors = list(
##     beta = normal(mean = 0, sd = 2.5),
##     nngp_1 = list(
##       sigma_sq = log_normal(meanlog = log(100), sdlog = 0.15),
##       alpha = uniform(0, 1),
##       phi_1 = uniform(phi_1_bounds[1], phi_1_bounds[2]),
##       lambda_1 = uniform(lambda_1_bounds[1], lambda_1_bounds[2]),
##       phi_2 = uniform(phi_2_bounds[1], phi_2_bounds[2]),
##       lambda_2 = uniform(lambda_2_bounds[1], lambda_2_bounds[2])
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   save_process = posterior_sub_sample,
##   n_omp_threads = n_omp_threads,
##   verbose = TRUE,
##   n_report = n_report
## )

## biomass_fit <- stLMM(
##   biomass_sqrt ~
##     tcc_mean_scaled +
##     nngp(x_km, y_km, time,
##          m = 15,
##          cov_model = "multi_res_sep_exp",
##          ordering = "maxmin",
##          st_scale = fit_st_scale) +
##     tcc_mean_scaled:nngp(x_km, y_km,
##          m = 15,
##          cov_model = "exp",
##          ordering = "maxmin") +
##     resid(model = "group", group = fia_unit),
##   data = positive_plots,
##   starting = biomass_starting,
##   priors = list(
##     beta = normal(mean = 0, sd = 2 * sqrt_response_sd),
##     resid = list(tau_sq = half_t(df = 3, scale = sqrt_response_sd)),
##     nngp_1 = list(
##       sigma_sq = half_t(df = 3, scale = sqrt_response_sd),
##       alpha = uniform(0, 1),
##       phi_1 = uniform(phi_1_bounds[1], phi_1_bounds[2]),
##       lambda_1 = uniform(lambda_1_bounds[1], lambda_1_bounds[2]),
##       phi_2 = uniform(phi_2_bounds[1], phi_2_bounds[2]),
##       lambda_2 = uniform(lambda_2_bounds[1], lambda_2_bounds[2])
##     ),
##     nngp_2 = list(
##       sigma_sq = half_t(df = 3, scale = 0.5 * sqrt_response_sd),
##       phi = uniform(phi_svc_bounds[1], phi_svc_bounds[2])
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   n_omp_threads = n_omp_threads,
##   verbose = TRUE,
##   n_report = n_report
## )



## occurrence_summary$parameters
## biomass_summary$parameters

## posterior_medians <- function(fit) {
##   draws <- as_samples(fit, burn = burnin, metadata = FALSE)
##   apply(draws, 2, stats::median)
## }
##
## param_value <- function(x, name) {
##   unname(x[name])
## }
##
## occurrence_medians <- posterior_medians(occurrence_fit)
## biomass_medians <- posterior_medians(biomass_fit)
##
## practical_ranges <- function(p, model, process) {
##   tibble(
##     model = model,
##     component = c("broad", "local"),
##     spatial_range_km = -log(cor_at_range) / c(
##       param_value(p, paste0(process, "_phi_1")),
##       param_value(p, paste0(process, "_phi_2"))
##     ),
##     temporal_range_years = -log(cor_at_range) / c(
##       param_value(p, paste0(process, "_lambda_1")),
##       param_value(p, paste0(process, "_lambda_2"))
##     ),
##     variance_weight = c(
##       param_value(p, paste0(process, "_alpha")),
##       1 - param_value(p, paste0(process, "_alpha"))
##     )
##   )
## }
##
## range_summary <- bind_rows(
##   practical_ranges(occurrence_medians, "positive-biomass occurrence", "nngp_1"),
##   practical_ranges(biomass_medians, "sqrt positive biomass intercept", "nngp_1"),
##   tibble(
##     model = "sqrt positive biomass TCC SVC",
##     component = "spatial",
##     spatial_range_km = -log(cor_at_range) / param_value(biomass_medians, "nngp_2_phi"),
##     temporal_range_years = NA_real_,
##     variance_weight = NA_real_
##   )
## ) |>
##   mutate(across(where(is.numeric), ~ round(.x, 2)))
##
## show_table(
##   range_summary,
##   caption = "Posterior median practical ranges for the intercept and spatial SVC covariance components."
## )

## biomass_draws <- as_samples(biomass_fit, burn = burnin, metadata = FALSE)
## resid_var_names <- paste0("tau_sq_", fia_unit_levels)
## if (!all(resid_var_names %in% colnames(biomass_draws))) {
##   resid_var_names <- fia_unit_levels
## }
## resid_sd_draws <- sqrt(biomass_draws[, resid_var_names, drop = FALSE])
##
## fia_resid_sd <- tibble(
##   fia_unit = factor(fia_unit_levels, levels = fia_unit_levels),
##   sd_mean = colMeans(resid_sd_draws),
##   sd_lower = apply(resid_sd_draws, 2, quantile, probs = 0.025),
##   sd_upper = apply(resid_sd_draws, 2, quantile, probs = 0.975)
## ) |>
##   left_join(fia_unit_names, by = "fia_unit")
##
## show_table(
##   fia_resid_sd |>
##     transmute(
##       `FIA unit` = fia_unit_name,
##       `posterior mean SD` = fmt_num(sd_mean, 2),
##       `lower` = fmt_num(sd_lower, 2),
##       `upper` = fmt_num(sd_upper, 2)
##     ),
##   caption = "FIA-unit residual standard deviations for square-root positive biomass."
## )

prediction_sample_matrix <- function(pred, sample, include_w = FALSE) {
  as.matrix(as_samples(
    pred,
    sample = sample,
    include_w = include_w,
    metadata = FALSE
  ))
}

prediction_w_matrix <- function(pred, term) {
  draws <- prediction_sample_matrix(pred, "mu", include_w = TRUE)
  w_cols <- grep(paste0("^w_", term, "_"), colnames(draws), value = TRUE)

  if (!length(w_cols)) {
    stop("No predicted process samples found for ", term, call. = FALSE)
  }

  draws[, w_cols, drop = FALSE]
}

aggregate_grid_draws <- function(draws, grid, prefix = "") {
  county_year <- grid |>
    distinct(county_fips, county, year) |>
    arrange(county_fips, year)

  county_samples <- matrix(
    NA_real_,
    nrow = nrow(draws),
    ncol = nrow(county_year)
  )

  for (j in seq_len(nrow(county_year))) {
    ii <- which(
      grid$county_fips == county_year$county_fips[j] &
        grid$year == county_year$year[j]
    )
    county_samples[, j] <- rowMeans(draws[, ii, drop = FALSE])
  }

  county_year |>
    bind_cols(
      summarize_draw_matrix(county_samples, prefix = prefix) |>
        select(-prediction_row)
    )
}

## occurrence_rec <- recover(
##   occurrence_fit,
##   sub_sample = posterior_sub_sample
## )
##
## biomass_rec <- recover(
##   biomass_fit,
##   sub_sample = posterior_sub_sample
## )



## occurrence_pred <- predict(
##   occurrence_rec,
##   newdata = prediction_grid,
##   y_samples = TRUE,
##   st_scale = pred_st_scale,
##   return_w_samples = FALSE
## )
##
## biomass_pred <- predict(
##   biomass_rec,
##   newdata = prediction_grid,
##   y_samples = TRUE,
##   return_w_samples = FALSE
## )
##
## svc_map_year <- 2024
## svc_prediction_grid <- prediction_grid |>
##   filter(year == svc_map_year)
##
## biomass_svc_pred <- predict(
##   biomass_rec,
##   newdata = svc_prediction_grid,
##   y_samples = FALSE,
##   return_w_samples = TRUE
## )



## direct_for_compare <- direct_estimates |>
##   filter(year %in% plot_years) |>
##   select(
##     county_fips, county, year, n,
##     direct_biomass_model, direct_biomass_se_model
##   )
##
## county_summary <- county_summary |>
##   mutate(county_fips = as.character(county_fips)) |>
##   left_join(
##     direct_for_compare,
##     by = c("county_fips", "county", "year")
##   )

## ggplot() +
##   geom_raster(
##     data = prediction_map,
##     aes(x, y, fill = theta_median)
##   ) +
##   geom_sf(data = wa_counties, fill = NA, color = "grey35", linewidth = 0.15) +
##   facet_wrap(~ year, nrow = 1) +
##   coord_sf(expand = FALSE) +
##   scale_fill_gradientn(
##     colors = stlmm_palette(),
##     name = "Mg/ha",
##     na.value = "grey92"
##   ) +
##   labs(
##     title = "posterior predictive median biomass",
##     x = NULL,
##     y = NULL
##   ) +
##   theme(
##     plot.title = element_text(size = 10, face = "bold", margin = margin(0, 0, 2, 0)),
##     plot.margin = margin(0, 0, 0, 0)
##   )

## ggplot() +
##   geom_raster(
##     data = tcc_coef_map,
##     aes(x, y, fill = tcc_coef_mean)
##   ) +
##   geom_sf(data = wa_counties, fill = NA, color = "grey35", linewidth = 0.15) +
##   coord_sf(expand = FALSE) +
##   scale_fill_gradientn(
##     colors = stlmm_continuous_palette(palette = "navia"),
##     name = "sqrt(Mg/ha)\nper SD TCC",
##     na.value = "grey92"
##   ) +
##   labs(
##     title = "posterior mean positive-stage TCC coefficient",
##     subtitle = paste("prediction grid,", svc_map_year),
##     x = NULL,
##     y = NULL
##   ) +
##   theme(
##     plot.title = element_text(size = 10, face = "bold", margin = margin(0, 0, 2, 0)),
##     plot.subtitle = element_text(size = 9, margin = margin(0, 0, 2, 0)),
##     plot.margin = margin(0, 0, 0, 0)
##   )

## show_table(
##   county_summary |>
##     filter(year == 2024) |>
##     arrange(desc(theta_mean)) |>
##     select(county, n, direct_biomass_model, theta_mean, theta_lower, theta_upper) |>
##     slice_head(n = 10) |>
##     mutate(
##       across(
##         c(direct_biomass_model, theta_mean, theta_lower, theta_upper),
##         ~ fmt_num(.x, 1)
##       )
##     ),
##   caption = "Highest 2024 county-year posterior predictive means from the continuous space-time model with spatial SVC."
## )

## profile_counties <- c("Clallam", "King", "Okanogan", "Yakima")
##
## profile_dat <- county_summary |>
##   filter(county %in% profile_counties) |>
##   mutate(county = factor(county, levels = profile_counties))
##
## ggplot(profile_dat, aes(year, theta_mean)) +
##   geom_ribbon(
##     aes(ymin = theta_lower, ymax = theta_upper),
##     fill = stlmm_color("primary"),
##     alpha = 0.15,
##     color = NA
##   ) +
##   geom_line(color = stlmm_color("primary"), linewidth = 0.75) +
##   geom_point(
##     aes(y = direct_biomass_model),
##     color = "grey20",
##     size = 1.4,
##     alpha = 0.75,
##     na.rm = TRUE
##   ) +
##   facet_wrap(~ county, scales = "free_y", ncol = 2) +
##   labs(
##     x = "year",
##     y = "mean agb_live_mg_ha (Mg/ha)"
##   )
