library(stLMM)
library(sf)
library(tidyverse)
library(knitr)
library(kableExtra)

source("article-utils.R")

set.seed(1)
theme_set(theme_bw(base_size = 12))

article_id <- "wa-unit-bhf-ladder"
data_dir <- "wa_data"

wa_counties <- read_rds(file.path(data_dir, "wa_counties.rds"))
unit_plots <- read_csv(file.path(data_dir, "wa_unit_plots.csv"), show_col_types = FALSE)
direct_estimates <- read_csv(file.path(data_dir, "wa_direct_estimates.csv"), show_col_types = FALSE)
prediction_grid <- read_csv(
  file.path(data_dir, "wa_unit_prediction_grid_10km.csv"),
  show_col_types = FALSE
)

wa_counties <- wa_counties |>
  mutate(county_fips = as.character(county_fips))

unit_plots <- unit_plots |>
  mutate(
    county_fips = as.character(county_fips),
    county_fips = factor(county_fips),
    year = as.integer(year),
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

prediction_grid <- prediction_grid |>
  mutate(
    county_fips = as.character(county_fips),
    county_fips = factor(county_fips, levels = county_levels),
    year = as.integer(year),
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
    "counties",
    "measurement years",
    "positive-biomass plot rows",
    "zero-biomass plot rows",
    "10 km prediction-grid rows",
    "plot TCC mean",
    "plot TCC SD"
  ),
  value = c(
    fmt_int(nrow(unit_plots)),
    fmt_int(n_distinct(unit_plots$county_fips)),
    paste(range(plot_years), collapse = "-"),
    fmt_int(sum(unit_plots$biomass_positive == 1)),
    fmt_int(sum(unit_plots$biomass_positive == 0)),
    fmt_int(nrow(prediction_grid)),
    fmt_num(tcc_center, 1),
    fmt_num(tcc_scale, 1)
  )
)

show_table(data_summary, caption = "Unit-response data used in the BHF ladder.")

tcc_bin_summary <- unit_plots |>
  filter(!is.na(tcc_mean)) |>
  mutate(tcc_bin = ntile(tcc_mean, 10)) |>
  group_by(tcc_bin) |>
  summarise(
    mean_tcc = mean(tcc_mean),
    positive_biomass_rate = mean(biomass_positive),
    mean_sqrt_positive_biomass = mean(
      sqrt(agb_live_mg_ha[biomass_positive == 1]),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    c(positive_biomass_rate, mean_sqrt_positive_biomass),
    names_to = "quantity",
    values_to = "value"
  ) |>
  mutate(
    quantity = recode(
      quantity,
      positive_biomass_rate = "positive-biomass rate (proportion)",
      mean_sqrt_positive_biomass = "mean sqrt positive biomass"
    )
  )

ggplot(tcc_bin_summary, aes(mean_tcc, value)) +
  geom_line(color = stlmm_color("primary"), linewidth = 0.75) +
  geom_point(color = stlmm_color("secondary"), size = 1.8) +
  facet_wrap(~ quantity, scales = "free_y", ncol = 1) +
  labs(
    x = "mean TCC within observed decile",
    y = "observed binned response"
  ) +
  theme(panel.spacing = grid::unit(0.45, "lines"))

year_summary <- unit_plots |>
  group_by(year) |>
  summarise(
    plots = n(),
    positive_biomass_rate = mean(biomass_positive),
    mean_response = mean(agb_live_mg_ha),
    median_response = median(agb_live_mg_ha),
    .groups = "drop"
  )

show_table(
  year_summary |>
    filter(year %in% c(2016, 2018, 2020, 2022, 2024)) |>
    mutate(
      positive_biomass_rate = fmt_num(positive_biomass_rate, 2),
      mean_response = fmt_num(mean_response, 1),
      median_response = fmt_num(median_response, 1)
    ),
  caption = "Selected-year summaries of the `agb_live_mg_ha` response."
)

n_samples <- 10000
burnin <- 5000
chains <- 3
n_keep <- 100
posterior_thin <- max(1L, floor((n_samples - burnin) / n_keep))
posterior_sub_sample <- list(start = burnin + 1, thin = posterior_thin)
chain_control <- list(seed = 7, dispersion = 1.5)
warmup_control <- list(batch_length = 25, min_batches = 10)

response_sd <- sd(unit_plots$agb_live_mg_ha, na.rm = TRUE)
sqrt_response_sd <- sd(positive_plots$biomass_sqrt, na.rm = TRUE)

fit_summary_parameters <- list(
  gaussian_ar1 = c(
    "(Intercept)", "tcc_mean_scaled", "tau_sq",
    "iid_1_sigma_sq", "ar1_1_sigma_sq", "ar1_1_phi"
  ),
  positive_ar1 = c(
    "(Intercept)", "tcc_mean_scaled",
    "iid_1_sigma_sq", "ar1_1_sigma_sq", "ar1_1_phi"
  ),
  positive_biomass_ar1 = c(
    "(Intercept)", "tcc_mean_scaled", "tau_sq",
    "iid_1_sigma_sq", "ar1_1_sigma_sq", "ar1_1_phi"
  )
)

select_summary_parameters <- function(x, parameters) {
  keep <- intersect(parameters, rownames(x$parameters))
  x$parameters <- x$parameters[keep, , drop = FALSE]
  if (!is.null(x$diagnostics)) {
    x$diagnostics <- x$diagnostics[x$diagnostics$parameter %in% keep, , drop = FALSE]
  }
  x
}



## gaussian_ar1_fit <- stLMM(
##   agb_live_mg_ha ~ tcc_mean_scaled + iid(county_fips) + ar1(year),
##   data = unit_plots,
##   priors = list(
##     beta = normal(mean = 0, sd = 2 * response_sd),
##     resid = list(tau_sq = half_t(df = 3, scale = response_sd)),
##     iid_1 = list(sigma_sq = ig(shape = 3, scale = 2 * response_sd^2)),
##     ar1_1 = list(
##       sigma_sq = half_t(df = 3, scale = response_sd),
##       phi = uniform(0.01, 0.99)
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   chain_control = chain_control,
##   warmup = warmup_control,
##   verbose = TRUE,
##   n_report = 500
## )
## 
## gaussian_ar1_summary <- summary(
##   gaussian_ar1_fit,
##   burn = burnin
## ) |>
##   select_summary_parameters(fit_summary_parameters$gaussian_ar1)

## positive_ar1_fit <- stLMM(
##   biomass_positive ~ tcc_mean_scaled + iid(county_fips) + ar1(year),
##   data = unit_plots,
##   family = "binomial",
##   priors = list(
##     beta = normal(mean = 0, sd = 2.5),
##     iid_1 = list(sigma_sq = ig(shape = 3, scale = 2)),
##     ar1_1 = list(
##       sigma_sq = half_t(df = 3, scale = 1),
##       phi = uniform(0.01, 0.99)
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   chain_control = chain_control,
##   warmup = warmup_control,
##   verbose = TRUE,
##   n_report = 500
## )
## 
## positive_ar1_summary <- summary(
##   positive_ar1_fit,
##   burn = burnin
## ) |>
##   select_summary_parameters(fit_summary_parameters$positive_ar1)
## 
## positive_biomass_ar1_fit <- stLMM(
##   biomass_sqrt ~ tcc_mean_scaled + iid(county_fips) + ar1(year),
##   data = positive_plots,
##   priors = list(
##     beta = normal(mean = 0, sd = 2 * sqrt_response_sd),
##     resid = list(tau_sq = half_t(df = 3, scale = sqrt_response_sd)),
##     iid_1 = list(sigma_sq = ig(shape = 3, scale = 2 * sqrt_response_sd^2)),
##     ar1_1 = list(
##       sigma_sq = half_t(df = 3, scale = sqrt_response_sd),
##       phi = uniform(0.01, 0.99)
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   chain_control = chain_control,
##   warmup = warmup_control,
##   verbose = TRUE,
##   n_report = 500
## )
## 
## positive_biomass_ar1_summary <- summary(
##   positive_biomass_ar1_fit,
##   burn = burnin
## ) |>
##   select_summary_parameters(fit_summary_parameters$positive_biomass_ar1)





## gaussian_ar1_summary
## positive_ar1_summary
## positive_biomass_ar1_summary

prediction_sample_matrix <- function(pred, sample) {
  as.matrix(as_samples(pred, sample = sample, metadata = FALSE))
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

## gaussian_ar1_rec <- recover(
##   gaussian_ar1_fit,
##   sub_sample = posterior_sub_sample
## )
## 
## positive_ar1_rec <- recover(
##   positive_ar1_fit,
##   sub_sample = posterior_sub_sample
## )
## 
## positive_biomass_ar1_rec <- recover(
##   positive_biomass_ar1_fit,
##   sub_sample = posterior_sub_sample
## )





## gaussian_ar1_pred <- predict(
##   gaussian_ar1_rec,
##   newdata = prediction_grid,
##   y_samples = TRUE
## )
## 
## positive_ar1_pred <- predict(
##   positive_ar1_rec,
##   newdata = prediction_grid,
##   y_samples = TRUE
## )
## 
## positive_biomass_ar1_pred <- predict(
##   positive_biomass_ar1_rec,
##   newdata = prediction_grid,
##   y_samples = TRUE
## )





## model_summaries <- bind_rows(
##   gaussian_ar1_county_summary |> mutate(model = "Gaussian BHF + TCC + AR(1)"),
##   two_stage_ar1_county_summary |> mutate(model = "two-stage + TCC + AR(1)")
## ) |>
##   mutate(
##     county_fips = as.character(county_fips),
##     model = factor(
##       model,
##       levels = c("Gaussian BHF + TCC + AR(1)", "two-stage + TCC + AR(1)")
##     )
##   )
## 
## direct_for_compare <- direct_estimates |>
##   filter(year %in% plot_years) |>
##   select(
##     county_fips, county, year, n,
##     direct_biomass_model, direct_biomass_se_model
##   )

## show_table(
##   model_summaries |>
##     filter(year == 2024) |>
##     left_join(
##       direct_for_compare |> filter(year == 2024),
##       by = c("county_fips", "county", "year")
##     ) |>
##     arrange(model, desc(theta_mean)) |>
##     group_by(model) |>
##     slice_head(n = 5) |>
##     ungroup() |>
##     select(model, county, n, direct_biomass_model, theta_mean, theta_lower, theta_upper) |>
##     mutate(
##       across(
##         c(direct_biomass_model, theta_mean, theta_lower, theta_upper),
##         ~ fmt_num(.x, 1)
##       )
##     ),
##   caption = "Highest 2024 county-year posterior predictive means by unit-response model."
## )

## panel_years <- c(2016, 2018, 2020, 2022, 2024)
## 
## map_dat <- wa_counties |>
##   left_join(
##     model_summaries |> filter(year %in% panel_years),
##     by = "county_fips"
##   )
## 
## ggplot(map_dat) +
##   geom_sf(aes(fill = theta_mean), color = "grey80", linewidth = 0.12) +
##   coord_sf(expand = FALSE) +
##   facet_grid(model ~ year) +
##   scale_fill_gradientn(
##     colors = stlmm_palette(),
##     name = "Mg/ha",
##     na.value = "grey92"
##   ) +
##   theme(
##     panel.spacing = grid::unit(0.03, "lines"),
##     strip.text = element_text(margin = margin(1, 1, 1, 1)),
##     plot.margin = margin(0, 0, 0, 0)
##   )

## profile_counties <- c("Clallam", "King", "Okanogan", "Yakima")
## 
## profile_dat <- model_summaries |>
##   filter(county %in% profile_counties) |>
##   mutate(county = factor(county, levels = profile_counties))
## 
## profile_direct <- direct_for_compare |>
##   filter(county %in% profile_counties) |>
##   mutate(county = factor(county, levels = profile_counties))
## 
## ggplot(profile_dat, aes(year, theta_mean, color = model, fill = model)) +
##   geom_ribbon(
##     aes(ymin = theta_lower, ymax = theta_upper),
##     alpha = 0.12,
##     color = NA
##   ) +
##   geom_line(linewidth = 0.75) +
##   geom_point(
##     data = profile_direct,
##     aes(year, direct_biomass_model),
##     inherit.aes = FALSE,
##     color = "grey20",
##     size = 1.4,
##     alpha = 0.75
##   ) +
##   facet_wrap(~ county, scales = "free_y", ncol = 2) +
##   scale_color_manual(values = stlmm_discrete_colors(2)) +
##   scale_fill_manual(values = stlmm_discrete_colors(2)) +
##   labs(
##     x = "year",
##     y = "mean agb_live_mg_ha (Mg/ha)",
##     color = NULL,
##     fill = NULL
##   ) +
##   theme(legend.position = "bottom")
