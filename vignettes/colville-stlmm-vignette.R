library(tidyverse)
library(sf)
library(Matrix)
library(stLMM)
library(RhpcBLASctl)

source("article-utils.R")

article_id <- "colville-stlmm-vignette"
colville_data_dir <- Sys.getenv(
  "STLMM_COLVILLE_DATA_DIR",
  unset = file.path(article_helper_dir, "colville_data", "derived")
)

set.seed(1)
theme_set(theme_bw(base_size = 12))

fia_plots <- st_read(file.path(colville_data_dir, "fia_plots_model_data.gpkg"),
                     layer = "fia_plots_model_data",
                     quiet = TRUE)

prediction_grid_eco <- st_read(file.path(colville_data_dir, "prediction_grid_eco_5km_1999_2026.gpkg"),
                               layer = "prediction_grid_eco_5km_1999_2026",
                               quiet = TRUE)

prediction_grid_nf <- st_read(file.path(colville_data_dir, "prediction_grid_nf_2km_1999_2026.gpkg"),
                              layer = "prediction_grid_nf_2km_1999_2026",
                              quiet = TRUE)

nf <- st_read(file.path(colville_data_dir, "colville_national_forest_northern_cascades.gpkg"),
              layer = "colville_national_forest_northern_cascades",
              quiet = TRUE)

eco <- st_read(file.path(colville_data_dir, "northern_cascades_ecoregion.gpkg"),
               layer = "northern_cascades_ecoregion",
               quiet = TRUE)

fia_dat <- fia_plots %>%
  st_drop_geometry() %>%
  mutate(
    forested = as.integer(forested)
  )

pred_eco_dat <- prediction_grid_eco %>%
  st_drop_geometry()

pred_nf_dat <- prediction_grid_nf %>%
  st_drop_geometry()

glimpse(fia_dat)
glimpse(pred_eco_dat)
glimpse(pred_nf_dat)

ggplot() +
  geom_sf(data = eco, fill = "grey94", color = "grey55", linewidth = 0.35) +
  geom_sf(data = nf, fill = NA, color = stlmm_color("secondary"), linewidth = 0.8) +
  geom_sf(data = fia_plots, aes(color = time), size = 0.7, alpha = 0.7) +
  scale_color_gradientn(colors = stlmm_palette(), name = "Measurement year") +
  coord_sf(datum = NA) +
  labs(x = NULL, y = NULL)

data_support_start <- min(fia_dat$year, na.rm = TRUE)
data_support_end <- max(fia_dat$year, na.rm = TRUE)
data_support_dense_end <- 2021
prediction_years <- sort(unique(pred_nf_dat$year))

plot_counts_by_year <- tibble(year = prediction_years) %>%
  left_join(
    fia_dat %>%
      group_by(year) %>%
      summarize(
        plots = n(),
        forested_plots = sum(forested == 1),
        .groups = "drop"
      ),
    by = "year"
  ) %>%
  mutate(across(c(plots, forested_plots), ~replace_na(.x, 0L))) %>%
  mutate(
    support = case_when(
      year < data_support_start ~ "no same-year FIA measurements",
      year == data_support_start ~ "first measurement year",
      year <= data_support_dense_end ~ "dense FIA support",
      year <= data_support_end ~ "sparse recent support",
      TRUE ~ "no same-year FIA measurements"
    )
  )

knitr::kable(
  plot_counts_by_year,
  caption = "FIA plot measurements by year in the Colville National Forest portion of the Northern Cascades ecoregion."
)

cor_at_range <- 0.05

decay_from_range <- function(range_05) {
  -log(cor_at_range) / range_05
}

## Broad process: approximately 25-100 km and 1-200 years.
phi_1_bounds <- decay_from_range(rev(c(25, 100)))
lambda_1_bounds <- decay_from_range(rev(c(1, 500)))

## Local process: approximately 1-25 km and 1-200 years.
phi_2_bounds <- decay_from_range(rev(c(1, 25)))
lambda_2_bounds <- decay_from_range(rev(c(1, 500)))

## Some chain thinning choices and metropolis controls for warmup.
n_chains <- 3
n_samples <- 2000
posterior_burn <- 1000
posterior_thin <- 10
prediction_draws_per_chain <- 25
prediction_thin <- ceiling((n_samples - posterior_burn) / prediction_draws_per_chain)
prediction_sub_sample <- list(start = posterior_burn + 1, thin = prediction_thin)
map_years <- c(1999, 2006, 2013, 2020, 2026)
n_report <- 1000

warmup_control <- list(
  batch_length = 25,
  min_batches = 10,
  max_batches = 20
)

## Request OpenMP threads for sampler and prediction work.
blas_set_num_threads(1)
n_omp_threads <- as.integer(Sys.getenv("STLMM_OMP_THREADS", unset = "4"))

## Space-time neighbor scaling for NNGP fitting and prediction. This treats one
## year as about 0.5 km when ordering nodes and finding neighbors.
fit_st_scale <- 0.5
pred_st_scale <- fit_st_scale



jitter_start <- function(value, amount = 0.02, bounds = NULL) {
  out <- rep(value, n_chains) * exp(rnorm(n_chains, 0, amount))

  if (!is.null(bounds)) {
    out <- pmin(pmax(out, bounds[1] + 1e-8), bounds[2] - 1e-8)
  }

  out
}

## forest_fit <- stLMM(
##   forested ~ nngp(x_km, y_km, time,
##                   m = 15,
##                   cov_model = "multi_res_sep_exp",
##                   ordering = "maxmin",
##                   st_scale = fit_st_scale),
##   data = fia_dat,
##   family = "binomial",
##   starting = list(
##     nngp_1 = list(
##       sigma_sq = jitter_start(100),
##       alpha = jitter_start(0.5, bounds = c(0, 1)),
##       phi_1 = jitter_start(0.075, bounds = phi_1_bounds),
##       lambda_1 = jitter_start(0.008, bounds = lambda_1_bounds),
##       phi_2 = jitter_start(1.5, bounds = phi_2_bounds),
##       lambda_2 = jitter_start(0.075, bounds = lambda_2_bounds)
##     )
##   ),
##   priors = list(
##     beta = normal(mean = 0, sd = 2),
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
##   chains = n_chains,
##   warmup = warmup_control,
##   save_process = prediction_sub_sample,
##   n_omp_threads = n_omp_threads,
##   n_report = n_report,
##   verbose = TRUE
## )
## 

biomass_dat <- fia_dat %>%
  filter(forested == 1) %>%
  mutate(biomass_sqrt = sqrt(biomass_Mg_ha))

## biomass_fit <- stLMM(
##   biomass_sqrt ~ nngp(x_km, y_km, time,
##                           m = 15,
##                           cov_model = "multi_res_sep_exp",
##                           ordering = "maxmin",
##                           st_scale = fit_st_scale),
##   data = biomass_dat,
##   family = "gaussian",
##   starting = list(
##     resid = list(
##       tau_sq = jitter_start(1)
##     ),
##     nngp_1 = list(
##       sigma_sq = jitter_start(45),
##       alpha = jitter_start(0.5, bounds = c(0, 1)),
##       phi_1 = jitter_start(0.05, bounds = phi_1_bounds),
##       lambda_1 = jitter_start(0.01, bounds = lambda_1_bounds),
##       phi_2 = jitter_start(2, bounds = phi_2_bounds),
##       lambda_2 = jitter_start(0.01, bounds = lambda_2_bounds)
##     )
##   ),
##   priors = list(
##     beta = normal(mean = mean(biomass_dat$biomass_sqrt), sd = 5),
##     resid = list(tau_sq = half_t(df = 3, scale = 1)),
##     nngp_1 = list(
##       sigma_sq = half_t(df = 3, scale = 2),
##       alpha = uniform(0, 1),
##       phi_1 = uniform(phi_1_bounds[1], phi_1_bounds[2]),
##       lambda_1 = uniform(lambda_1_bounds[1], lambda_1_bounds[2]),
##       phi_2 = uniform(phi_2_bounds[1], phi_2_bounds[2]),
##       lambda_2 = uniform(lambda_2_bounds[1], lambda_2_bounds[2])
##     )
##   ),
##   n_samples = n_samples,
##   chains = n_chains,
##   warmup = warmup_control,
##   n_omp_threads = n_omp_threads,
##   n_report = n_report,
##   verbose = TRUE
## )
## 



forest_summary
biomass_summary

posterior_medians <- function(fit) {
  samples <- as.matrix(as_mcmc(
    fit,
    burn = posterior_burn,
    thin = posterior_thin
  ))

  apply(samples, 2, stats::median)
}

param_value <- function(x, name) {
  unname(x[name])
}

forest_medians <- posterior_medians(forest_fit)
biomass_medians <- posterior_medians(biomass_fit)

practical_ranges <- function(p, model) {
  tibble(
    model = model,
    component = c("broad", "local"),
    spatial_range_km = -log(cor_at_range) / c(param_value(p, "nngp_1_phi_1"),
                                               param_value(p, "nngp_1_phi_2")),
    temporal_range_years = -log(cor_at_range) / c(param_value(p, "nngp_1_lambda_1"),
                                                  param_value(p, "nngp_1_lambda_2")),
    variance_weight = c(param_value(p, "nngp_1_alpha"),
                        1 - param_value(p, "nngp_1_alpha"))
  )
}

range_summary <- bind_rows(
  practical_ranges(forest_medians, "Forested probability"),
  practical_ranges(biomass_medians, "Forested-plot biomass")
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

knitr::kable(
  range_summary,
  caption = "Posterior median practical ranges for the broad and local space-time covariance components."
)

sigma_sq <- param_value(biomass_medians, "nngp_1_sigma_sq")
alpha <- param_value(biomass_medians, "nngp_1_alpha")
tau_sq <- param_value(biomass_medians, "tau_sq")

signal_noise_summary <- tibble(
  component = c("total process", "broad component", "local component"),
  variance_ratio = c(sigma_sq / tau_sq,
                     alpha * sigma_sq / tau_sq,
                     (1 - alpha) * sigma_sq / tau_sq)
) %>%
  mutate(variance_ratio = round(variance_ratio, 2))

knitr::kable(
  signal_noise_summary,
  caption = "Posterior median biomass process-to-residual variance ratios."
)

predict_support_joint <- function(pred_dat) {
  list(
    forest = predict(
      forest_rec,
      newdata = pred_dat,
      y_samples = TRUE,
      joint = TRUE,
      joint_method = "vecchia",
      pred_m = 15,
      pred_ordering = "maxmin",
      st_scale = pred_st_scale,
      return_w_samples = FALSE
    ),
    biomass = predict(
      biomass_rec,
      newdata = pred_dat,
      y_samples = TRUE,
      joint = TRUE,
      joint_method = "vecchia",
      pred_m = 15,
      pred_ordering = "maxmin",
      st_scale = pred_st_scale,
      return_w_samples = FALSE
    )
  )
}

predict_support_independent <- function(pred_dat) {
  list(
    forest = predict(
      forest_rec,
      newdata = pred_dat,
      y_samples = TRUE,
      joint = FALSE,
      pred_m = 15,
      st_scale = pred_st_scale,
      return_w_samples = FALSE
    ),
    biomass = predict(
      biomass_rec,
      newdata = pred_dat,
      y_samples = TRUE,
      joint = FALSE,
      pred_m = 15,
      st_scale = pred_st_scale,
      return_w_samples = FALSE
    )
  )
}


## forest_rec <- recover(forest_fit, sub_sample = prediction_sub_sample)
## biomass_rec <- recover(biomass_fit, sub_sample = prediction_sub_sample)
## 
## eco_pred <- predict_support_joint(pred_eco_dat)
## nf_pred <- predict_support_joint(pred_nf_dat)
## nf_pred_independent <- predict_support_independent(pred_nf_dat)

prediction_sample_matrix <- function(pred, sample) {
  as.matrix(as_samples(pred, sample = sample, metadata = FALSE))
}

prediction_product <- function(prediction_grid, pred) {
  forest_prob_draws <- prediction_sample_matrix(pred$forest, "mu")
  forest_draws <- prediction_sample_matrix(pred$forest, "y")
  biomass_sqrt_draws <- prediction_sample_matrix(pred$biomass, "y")

  biomass_positive_draws <- biomass_sqrt_draws^2
  biomass_ppd_draws <- forest_draws * biomass_positive_draws

  list(
    map = prediction_grid %>%
      mutate(
        forest_prob_median = apply(forest_prob_draws, 2, stats::median),
        biomass_positive_median = apply(biomass_positive_draws, 2, stats::median),
        biomass_ppd_median = apply(biomass_ppd_draws, 2, stats::median)
      ),
    biomass_ppd_draws = biomass_ppd_draws
  )
}



## eco_product <- prediction_product(prediction_grid_eco, eco_pred)
## nf_product <- prediction_product(prediction_grid_nf, nf_pred)
## nf_product_independent <- prediction_product(prediction_grid_nf, nf_pred_independent)
## 
## prediction_map_eco <- eco_product$map
## prediction_map_nf <- nf_product$map

prediction_map_eco_years <- prediction_map_eco %>%
  filter(year %in% map_years)

prediction_map_nf_years <- prediction_map_nf %>%
  filter(year %in% map_years)

biomass_map_limits <- range(
  c(prediction_map_eco$biomass_ppd_median,
    prediction_map_nf$biomass_positive_median,
    prediction_map_nf$biomass_ppd_median)
)

map_theme <- theme_void(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    plot.title.position = "plot"
  )

ggplot() +
  geom_sf(data = eco, fill = "grey95", color = NA) +
  geom_raster(data = prediction_map_eco_years,
              aes(x = x, y = y, fill = biomass_ppd_median)) +
  geom_sf(data = eco, fill = NA, color = "grey35", linewidth = 0.3) +
  geom_sf(data = nf, fill = NA, color = stlmm_color("secondary"), linewidth = 0.35) +
  facet_wrap(vars(year), nrow = 1) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "Mg/ha",
    limits = biomass_map_limits
  ) +
  coord_sf(datum = NA) +
  labs(title = "Ecoregion posterior predictive median biomass product") +
  map_theme

ggplot() +
  geom_sf(data = nf, fill = "grey95", color = NA) +
  geom_raster(data = prediction_map_nf_years,
              aes(x = x, y = y, fill = forest_prob_median)) +
  geom_sf(data = nf, fill = NA, color = "grey20", linewidth = 0.35) +
  facet_wrap(vars(year), nrow = 1) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "Median probability",
    limits = c(0, 1)
  ) +
  coord_sf(datum = NA) +
  labs(title = "Posterior predictive median FIA-forested probability") +
  map_theme

ggplot() +
  geom_sf(data = nf, fill = "grey95", color = NA) +
  geom_raster(data = prediction_map_nf_years,
              aes(x = x, y = y, fill = biomass_positive_median)) +
  geom_sf(data = nf, fill = NA, color = "grey20", linewidth = 0.35) +
  facet_wrap(vars(year), nrow = 1) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "Mg/ha",
    limits = biomass_map_limits
  ) +
  coord_sf(datum = NA) +
  labs(title = "Posterior predictive median positive biomass") +
  map_theme

ggplot() +
  geom_sf(data = nf, fill = "grey95", color = NA) +
  geom_raster(data = prediction_map_nf_years,
              aes(x = x, y = y, fill = biomass_ppd_median)) +
  geom_sf(data = nf, fill = NA, color = "grey20", linewidth = 0.35) +
  facet_wrap(vars(year), nrow = 1) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "Mg/ha",
    limits = biomass_map_limits
  ) +
  coord_sf(datum = NA) +
  labs(title = "Posterior predictive median biomass product") +
  map_theme

annual_summary <- function(biomass_ppd_draws, pred_dat, support, prediction_method) {
  grid_year <- pred_dat$year

  annual_draws <- sapply(sort(unique(grid_year)), function(yy) {
    rowMeans(biomass_ppd_draws[, grid_year == yy, drop = FALSE])
  })

  as_tibble(t(apply(
    annual_draws,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975)
  ))) %>%
    set_names(c("q025", "median", "q975")) %>%
    mutate(year = sort(unique(grid_year)),
           support = support,
           prediction_method = prediction_method) %>%
    relocate(support, prediction_method, year)
}

annual_summary_eco <- annual_summary(eco_product$biomass_ppd_draws,
                                     pred_eco_dat,
                                     "Northern Cascades ecoregion",
                                     "Joint Vecchia")

annual_summary_nf <- annual_summary(nf_product$biomass_ppd_draws,
                                    pred_nf_dat,
                                    "Colville National Forest",
                                    "Joint Vecchia")

annual_summary_nf_independent <- annual_summary(
  nf_product_independent$biomass_ppd_draws,
  pred_nf_dat,
  "Colville National Forest",
  "Independent NNGP"
)

annual_summaries <- bind_rows(annual_summary_eco,
                              annual_summary_nf,
                              annual_summary_nf_independent)

annual_summaries

support_levels <- c("Northern Cascades ecoregion", "Colville National Forest")

data_support_lines <- tibble(
  year = c(data_support_start, data_support_dense_end),
  label = c("first FIA measurements", "dense FIA support ends")
)

data_support_labels <- data_support_lines %>%
  mutate(support = factor("Northern Cascades ecoregion",
                          levels = support_levels))

annual_summaries_plot <- annual_summaries %>%
  mutate(
    support = factor(support, levels = support_levels),
    prediction_method = factor(
      prediction_method,
      levels = c("Joint Vecchia", "Independent NNGP")
    )
  )

prediction_method_colors <- c(
  "Joint Vecchia" = stlmm_color("primary"),
  "Independent NNGP" = stlmm_color("contrast")
)

ggplot(annual_summaries_plot, aes(x = year)) +
  geom_vline(
    data = data_support_lines,
    aes(xintercept = year),
    inherit.aes = FALSE,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.4
  ) +
  geom_ribbon(aes(ymin = q025, ymax = q975, fill = prediction_method),
              alpha = 0.25,
              color = NA) +
  geom_line(aes(y = median,
                color = prediction_method,
                linetype = prediction_method),
            linewidth = 0.9) +
  geom_point(aes(y = median, color = prediction_method), size = 1.5) +
  geom_text(
    data = data_support_labels,
    aes(x = year, y = Inf, label = label),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 1.45,
    vjust = -0.35,
    size = 3,
    color = "grey25"
  ) +
  scale_color_manual(values = prediction_method_colors) +
  scale_fill_manual(values = prediction_method_colors) +
  scale_linetype_manual(values = c("Joint Vecchia" = "solid",
                                   "Independent NNGP" = "dashed")) +
  scale_x_continuous(breaks = seq(1999, 2026, by = 3)) +
  facet_wrap(vars(support), ncol = 1) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Year",
    y = "Live aboveground biomass density (Mg/ha)",
    title = "Posterior predictive biomass distributions",
    subtitle = "The ecoregion uses joint Vecchia only;\nthe National Forest also shows the independent NNGP approximation",
    color = NULL,
    fill = NULL,
    linetype = NULL
  )
