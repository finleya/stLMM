library(stLMM)
library(sf)
library(tidyverse)
library(knitr)
library(kableExtra)

source("article-utils.R")

set.seed(1)
theme_set(theme_bw(base_size = 12))

article_id <- "wa-unit-two-stage-car-time-tcc-fiaunit-var"
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
    "5 km prediction-grid rows"
  ),
  value = c(
    fmt_int(nrow(unit_plots)),
    fmt_int(nrow(positive_plots)),
    fmt_int(n_distinct(unit_plots$fia_unit)),
    fmt_int(nrow(prediction_grid))
  )
)

show_table(data_summary, caption = "Unit-response data for the FIA-unit variance model.")

g <- car_graph(
  wa_counties,
  id = "county_fips",
  island = "nearest",
  island_k = 4
)

g$island_added_edges

n_samples <- 10000
burnin <- 5000
chains <- 3
n_keep <- 100
posterior_thin <- max(1L, floor((n_samples - burnin) / n_keep))
posterior_sub_sample <- list(start = burnin + 1, thin = posterior_thin)
chain_control <- list(seed = 17, dispersion = 1.5)
warmup_control <- list(batch_length = 25, min_batches = 10)

sqrt_response_sd <- sd(positive_plots$biomass_sqrt, na.rm = TRUE)



## occurrence_fit <- stLMM(
##   biomass_positive ~
##     tcc_mean_scaled +
##     car_time(county_fips, year, graph = g, car_model = "leroux"),
##   data = unit_plots,
##   family = "binomial",
##   priors = list(
##     beta = normal(mean = 0, sd = 2.5),
##     car_time_1 = list(
##       sigma_sq = half_t(df = 3, scale = 1),
##       rho = uniform(0.01, 0.99),
##       phi = uniform(-0.99, 0.99)
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   chain_control = chain_control,
##   warmup = warmup_control,
##   verbose = TRUE,
##   n_report = 500
## )

## biomass_fit <- stLMM(
##   biomass_sqrt ~
##     tcc_mean_scaled +
##     car_time(county_fips, year, graph = g, car_model = "leroux") +
##     resid(model = "group", group = fia_unit),
##   data = positive_plots,
##   priors = list(
##     beta = normal(mean = 0, sd = 2 * sqrt_response_sd),
##     resid = list(tau_sq = half_t(df = 3, scale = sqrt_response_sd)),
##     car_time_1 = list(
##       sigma_sq = half_t(df = 3, scale = sqrt_response_sd),
##       rho = uniform(0.01, 0.99),
##       phi = uniform(-0.99, 0.99)
##     )
##   ),
##   n_samples = n_samples,
##   chains = chains,
##   chain_control = chain_control,
##   warmup = warmup_control,
##   verbose = TRUE,
##   n_report = 500
## )



## occurrence_summary
## biomass_summary

## biomass_draws <- as_samples(biomass_fit, burn = burnin, metadata = FALSE)
## resid_var_names <- paste0("tau_sq_", fia_unit_levels)
## if (!all(resid_var_names %in% colnames(biomass_draws))) {
##   resid_var_names <- fia_unit_levels
## }
## resid_var_draws <- biomass_draws[, resid_var_names, drop = FALSE]
## resid_sd_draws <- sqrt(resid_var_draws)
## 
## fia_resid_sd <- tibble(
##   fia_unit = factor(fia_unit_levels, levels = fia_unit_levels),
##   sd_mean = colMeans(resid_sd_draws),
##   sd_lower = apply(resid_sd_draws, 2, quantile, probs = 0.025),
##   sd_upper = apply(resid_sd_draws, 2, quantile, probs = 0.975)
## ) |>
##   left_join(fia_unit_names, by = "fia_unit")
## 
## fia_resid_map <- wa_counties |>
##   select(-fia_unit_name) |>
##   left_join(fia_resid_sd, by = "fia_unit")
## 
## fia_resid_labels <- fia_resid_map |>
##   group_by(fia_unit, fia_unit_name, sd_mean) |>
##   summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
##   mutate(label = paste0(fia_unit_name, "\n", fmt_num(sd_mean, 2))) |>
##   sf::st_point_on_surface()
## 
## ggplot(fia_resid_map) +
##   geom_sf(aes(fill = sd_mean), color = "grey80", linewidth = 0.15) +
##   geom_sf_label(
##     data = fia_resid_labels,
##     aes(label = label),
##     size = 3,
##     label.size = 0.12,
##     fill = "white",
##     alpha = 0.85
##   ) +
##   coord_sf(expand = FALSE) +
##   scale_fill_gradientn(
##     colors = stlmm_palette(),
##     name = "resid. SD",
##     na.value = "grey92"
##   ) +
##   labs(
##     title = "positive-stage residual SD by FIA unit",
##     x = NULL,
##     y = NULL
##   ) +
##   theme(
##     plot.title = element_text(size = 10, face = "bold", margin = margin(0, 0, 2, 0)),
##     plot.margin = margin(0, 0, 0, 0)
##   )

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
##   y_samples = TRUE
## )
## 
## biomass_pred <- predict(
##   biomass_rec,
##   newdata = prediction_grid,
##   y_samples = TRUE
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
##   caption = "Highest 2024 county-year posterior predictive means from the FIA-unit variance model."
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
