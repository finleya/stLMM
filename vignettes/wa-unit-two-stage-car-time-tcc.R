library(stLMM)
library(sf)
library(tidyverse)
library(knitr)
library(kableExtra)
library(patchwork)

source("article-utils.R")

set.seed(1)
theme_set(theme_bw(base_size = 12))

article_id <- "wa-unit-two-stage-car-time-tcc"
data_dir <- Sys.getenv(
  "STLMM_WA_DATA_DIR",
  unset = file.path(article_helper_dir, "wa_data")
)

wa_counties <- read_rds(file.path(data_dir, "wa_counties.rds"))
unit_plots <- read_csv(file.path(data_dir, "wa_unit_plots.csv"), show_col_types = FALSE)
direct_estimates <- read_csv(file.path(data_dir, "wa_direct_estimates.csv"), show_col_types = FALSE)
prediction_grid <- read_csv(
  file.path(data_dir, "wa_unit_prediction_grid_5km.csv"),
  show_col_types = FALSE
)

wa_counties <- wa_counties |>
  mutate(county_fips = as.character(county_fips))

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
    "counties",
    "measurement years",
    "positive-biomass plot rows",
    "zero-biomass plot rows",
    "5 km prediction-grid rows",
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

show_table(data_summary, caption = "Unit-response data for the two-stage CAR-time TCC model.")

g <- car_graph(
  wa_counties,
  id = "county_fips",
  island = "nearest",
  island_k = 4
)

g$island_added_edges

coord_dat <- data.frame(
  county_fips = wa_counties$county_fips,
  sf::st_coordinates(sf::st_point_on_surface(sf::st_geometry(wa_counties)))
)

edge_index <- which(as.matrix(g$adjacency) != 0, arr.ind = TRUE)
edge_index <- edge_index[edge_index[, "row"] < edge_index[, "col"], , drop = FALSE]
edge_dat <- tibble(
  from = g$ids[edge_index[, "row"]],
  to = g$ids[edge_index[, "col"]]
) |>
  mutate(
    key = paste(pmin(from, to), pmax(from, to), sep = "--"),
    x = coord_dat$X[match(from, coord_dat$county_fips)],
    y = coord_dat$Y[match(from, coord_dat$county_fips)],
    xend = coord_dat$X[match(to, coord_dat$county_fips)],
    yend = coord_dat$Y[match(to, coord_dat$county_fips)]
  )

island_key <- paste(
  pmin(g$island_added_edges$from, g$island_added_edges$to),
  pmax(g$island_added_edges$from, g$island_added_edges$to),
  sep = "--"
)
edge_dat$island <- edge_dat$key %in% island_key

ggplot(wa_counties) +
  geom_sf(fill = "grey96", color = "white", linewidth = 0.15) +
  geom_segment(
    data = edge_dat,
    aes(x = x, y = y, xend = xend, yend = yend, color = island, linewidth = island)
  ) +
  geom_point(data = coord_dat, aes(X, Y), color = stlmm_color("primary"), size = 1.3) +
  coord_sf(expand = FALSE) +
  scale_color_manual(
    values = c("FALSE" = "grey45", "TRUE" = stlmm_color("secondary")),
    guide = "none"
  ) +
  scale_linewidth_manual(values = c("FALSE" = 0.25, "TRUE" = 1), guide = "none") +
  theme(plot.margin = margin(0, 0, 0, 0))

n_samples <- 10000
burnin <- 5000
chains <- 3
n_keep <- 100
posterior_thin <- max(1L, floor((n_samples - burnin) / n_keep))
posterior_sub_sample <- list(start = burnin + 1, thin = posterior_thin)
chain_control <- list(seed = 13, dispersion = 1.5)
warmup_control <- list(batch_length = 25, min_batches = 10)

sqrt_response_sd <- sd(positive_plots$biomass_sqrt, na.rm = TRUE)

fit_summary_parameters <- list(
  occurrence = c(
    "(Intercept)", "tcc_mean_scaled",
    "car_time_1_sigma_sq", "car_time_1_rho", "car_time_1_phi"
  ),
  biomass = c(
    "(Intercept)", "tcc_mean_scaled", "tau_sq",
    "car_time_1_sigma_sq", "car_time_1_rho", "car_time_1_phi"
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
## 
## occurrence_summary <- summary(
##   occurrence_fit,
##   burn = burnin
## ) |>
##   select_summary_parameters(fit_summary_parameters$occurrence)

## biomass_fit <- stLMM(
##   biomass_sqrt ~
##     tcc_mean_scaled +
##     car_time(county_fips, year, graph = g, car_model = "leroux"),
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
## 
## biomass_summary <- summary(
##   biomass_fit,
##   burn = burnin
## ) |>
##   select_summary_parameters(fit_summary_parameters$biomass)



## occurrence_summary
## biomass_summary

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
##     select(
##       county, n, direct_biomass_model,
##       theta_mean, theta_lower, theta_upper,
##       p_positive_mean, positive_mean
##     ) |>
##     slice_head(n = 10) |>
##     mutate(
##       across(
##         c(
##           direct_biomass_model,
##           theta_mean, theta_lower, theta_upper,
##           p_positive_mean, positive_mean
##         ),
##         ~ fmt_num(.x, 1)
##       )
##     ),
##   caption = "Highest 2024 county-year posterior predictive means from the two-stage CAR-time TCC model."
## )

## panel_years <- c(2005, 2010, 2015, 2020, 2024)
## 
## map_dat <- wa_counties |>
##   left_join(
##     county_summary |> filter(year %in% panel_years),
##     by = "county_fips"
##   )
## 
## biomass_map_limits <- c(0, max(map_dat$positive_mean, na.rm = TRUE))
## 
## component_map <- function(fill_var, legend_title, title, limits = NULL) {
##   ggplot(map_dat) +
##     geom_sf(aes(fill = .data[[fill_var]]), color = "grey80", linewidth = 0.12) +
##     coord_sf(expand = FALSE) +
##     facet_wrap(~ year, nrow = 1) +
##     scale_fill_gradientn(
##       colors = stlmm_palette(),
##       name = legend_title,
##       limits = limits,
##       na.value = "grey92"
##     ) +
##     labs(title = title) +
##     theme(
##       panel.spacing = grid::unit(0.03, "lines"),
##       strip.text = element_text(margin = margin(1, 1, 1, 1)),
##       plot.title = element_text(size = 10, face = "bold", margin = margin(0, 0, 2, 0)),
##       plot.margin = margin(0, 0, 0, 0)
##     )
## }
## 
## component_map(
##   "p_positive_mean",
##   "prob.",
##   "positive-biomass probability"
## ) /
##   component_map(
##     "positive_mean",
##     "Mg/ha",
##     "positive-biomass magnitude",
##     limits = biomass_map_limits
##   ) /
##   component_map(
##     "theta_mean",
##     "Mg/ha",
##     "recombined county-year mean",
##     limits = biomass_map_limits
##   )

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
