library(stLMM)
library(sf)
library(tidyverse)
library(knitr)
library(kableExtra)

source("article-utils.R")

set.seed(1)
theme_set(theme_bw(base_size = 12))

article_id <- "wa-direct-car-time-scaled-variance"
data_dir <- Sys.getenv(
  "STLMM_WA_DATA_DIR",
  unset = file.path(article_helper_dir, "wa_data")
)

wa_counties <- read_rds(file.path(data_dir, "wa_counties.rds"))
direct_estimates <- read_csv(file.path(data_dir, "wa_direct_estimates.csv"), show_col_types = FALSE)

wa_counties <- wa_counties |>
  mutate(county_fips = as.character(county_fips))

direct_estimates <- direct_estimates |>
  mutate(
    county_fips = as.character(county_fips),
    year = as.integer(year)
  ) |>
  arrange(county_fips, year)

observed_direct <- direct_estimates |>
  filter(!is.na(direct_biomass_model))

direct_summary <- tibble(
  quantity = c(
    "counties",
    "years",
    "county-years",
    "model-ready direct estimates",
    "missing model responses",
    "single-plot rows retained but excluded",
    "variance-floor rows",
    "median plot count among modeled county-years",
    "median direct-estimate SE"
  ),
  value = c(
    n_distinct(direct_estimates$county_fips),
    n_distinct(direct_estimates$year),
    nrow(direct_estimates),
    sum(direct_estimates$direct_estimate_in_model),
    sum(is.na(direct_estimates$direct_biomass_model)),
    sum(direct_estimates$direct_estimate_status == "single_plot_no_variance"),
    sum(direct_estimates$direct_estimate_in_model &
          direct_estimates$direct_biomass_vhat_source == "floor", na.rm = TRUE),
    median(observed_direct$n, na.rm = TRUE),
    median(observed_direct$direct_biomass_se_model, na.rm = TRUE)
  )
) |>
  mutate(value = fmt_num(as.numeric(value), 1))

show_table(direct_summary, caption = "County-year direct-estimate data used in the scaled-variance CAR-time model.")

year_coverage <- direct_estimates |>
  group_by(year) |>
  summarise(
    counties_with_direct_estimates = sum(direct_estimate_in_model),
    median_plot_n = median(n[n > 0], na.rm = TRUE),
    median_direct_biomass = median(direct_biomass_model, na.rm = TRUE),
    median_direct_se = median(direct_biomass_se_model, na.rm = TRUE),
    .groups = "drop"
  )

show_table(
  year_coverage |>
    filter(year %in% c(1999, 2005, 2010, 2015, 2020, 2024, 2025)) |>
    mutate(
      median_plot_n = fmt_num(median_plot_n, 0),
      median_direct_biomass = fmt_num(median_direct_biomass, 1),
      median_direct_se = fmt_num(median_direct_se, 1)
    ),
  caption = "Selected-year direct-estimate coverage and precision."
)

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
  scale_linewidth_manual(values = c("FALSE" = 0.25, "TRUE" = 1), guide = "none")

n_samples <- 12000
burnin <- 6000
chains <- 3
n_keep <- 100
posterior_thin <- max(1L, floor((n_samples - burnin) / n_keep))
posterior_sub_sample <- list(start = burnin + 1, thin = posterior_thin)
chain_control <- list(seed = 1, dispersion = 1.5)
warmup_control <- list(batch_length = 25, min_batches = 10)
summary_parameters <- c(
  "(Intercept)",
  "kappa",
  "tau0_sq",
  "car_time_1_sigma_sq",
  "car_time_1_rho",
  "car_time_1_phi"
)



## fit <- stLMM(
##   direct_biomass_model ~
##     car_time(county_fips, year, graph = g, car_model = "leroux") +
##     resid(
##       model = "scaled",
##       variance = direct_biomass_vhat_model,
##       n = n_eff_model,
##       shrinkage = 10,
##       kappa_log_prior = c(mean = 0, sd = 1)
##     ),
##   data = direct_estimates,
##   priors = list(
##     car_time_1 = list(
##       sigma_sq = half_t(
##         df = 3,
##         scale = sd(observed_direct$direct_biomass_model, na.rm = TRUE)
##       ),
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
## fit_summary <- summary(fit, burn = burnin, parameters = summary_parameters)



fit_summary

variance_draws <- as_samples(fit, burn = burnin, metadata = FALSE)

scaled_variance_summary <- observed_direct |>
  mutate(
    variance_weight = n_eff_model / (n_eff_model + 10),
    tau_sq_mean = purrr::map2_dbl(
      direct_biomass_vhat_model,
      variance_weight,
      \(vhat_i, weight_i) {
        mean(
          variance_draws$kappa *
            exp(
              weight_i * log(vhat_i) +
                (1 - weight_i) * log(variance_draws$tau0_sq)
            )
        )
      }
    )
  )

variance_axis_limits <- range(
  c(
    scaled_variance_summary$direct_biomass_vhat_model,
    scaled_variance_summary$tau_sq_mean
  ),
  na.rm = TRUE
)
variance_axis_limits <- 10^c(
  floor(log10(variance_axis_limits[1])),
  ceiling(log10(variance_axis_limits[2]))
)

ggplot(
  scaled_variance_summary,
  aes(x = direct_biomass_vhat_model, y = tau_sq_mean)
) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.45) +
  geom_point(aes(color = n_eff_model), alpha = 0.65, size = 1.7) +
  scale_x_log10(
    limits = variance_axis_limits,
    labels = scales::label_number(big.mark = ",", trim = TRUE)
  ) +
  scale_y_log10(
    limits = variance_axis_limits,
    labels = scales::label_number(big.mark = ",", trim = TRUE)
  ) +
  scale_color_viridis_c(option = "C", name = "Plot n") +
  coord_fixed() +
  labs(
    x = "Supplied direct-estimate variance",
    y = "Posterior mean scaled observation variance"
  )

## rec <- recover(
##   fit,
##   sub_sample = posterior_sub_sample
## )



## fitted_draws <- as.matrix(as_samples(fitted(rec, summary = FALSE), metadata = FALSE))
## 
## county_year_summary <- direct_estimates |>
##   bind_cols(
##     summarize_draw_matrix(fitted_draws, prefix = "theta_") |>
##       select(-prediction_row)
##   )



show_table(
  county_year_summary |>
    filter(year == 2024) |>
    arrange(desc(theta_mean)) |>
    select(
      county, n, direct_biomass_model, direct_biomass_se_model,
      theta_mean, theta_lower, theta_upper
    ) |>
    slice_head(n = 10) |>
    mutate(
      across(
        c(direct_biomass_model, direct_biomass_se_model, theta_mean, theta_lower, theta_upper),
        ~ fmt_num(.x, 1)
      )
    ),
  caption = "Highest 2024 model-smoothed county-year biomass means, with 95% posterior credible intervals."
)

panel_years <- c(2005, 2010, 2015, 2020, 2024)

map_dat <- wa_counties |>
  left_join(
    county_year_summary |> filter(year %in% panel_years),
    by = "county_fips"
  )

map_long <- bind_rows(
  map_dat |> mutate(quantity = "direct estimate", value = direct_biomass_model),
  map_dat |> mutate(quantity = "model-smoothed\nmean", value = theta_mean)
) |>
  mutate(quantity = factor(quantity, levels = c("direct estimate", "model-smoothed\nmean")))

ggplot(map_long) +
  geom_sf(aes(fill = value), color = "grey80", linewidth = 0.12) +
  coord_sf(expand = FALSE) +
  facet_grid(quantity ~ year) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "Mg/ha",
    na.value = "grey92"
  ) +
  theme(
    panel.spacing = grid::unit(0.02, "lines"),
    strip.text = element_text(margin = margin(1, 1, 1, 1)),
    plot.margin = margin(0, 0, 0, 0)
  )

scatter_dat <- county_year_summary |>
  filter(!is.na(direct_biomass_model))

fit_axis_limits <- range(
  c(scatter_dat$direct_biomass_model, scatter_dat$theta_mean, 0),
  na.rm = TRUE
)

ggplot(scatter_dat, aes(x = direct_biomass_model, y = theta_mean)) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.45) +
  geom_point(
    aes(size = 1 / direct_biomass_vhat_model, color = year),
    alpha = 0.55
  ) +
  scale_size_continuous(range = c(0.5, 3.5), guide = "none") +
  scale_color_viridis_c(option = "C") +
  coord_equal(xlim = fit_axis_limits, ylim = fit_axis_limits) +
  labs(
    x = "Direct estimate Mg/ha",
    y = "CAR-time posterior mean Mg/ha",
    color = "Year"
  )

selected_counties <- c("Clallam", "King", "Okanogan", "Yakima")

series_dat <- county_year_summary |>
  filter(county %in% selected_counties) |>
  mutate(county = factor(county, levels = selected_counties))

ggplot(series_dat, aes(x = year)) +
  geom_ribbon(
    aes(ymin = theta_lower, ymax = theta_upper),
    fill = "#9ecae1",
    alpha = 0.45
  ) +
  geom_line(aes(y = theta_mean), color = stlmm_color("primary"), linewidth = 0.7) +
  geom_point(
    aes(y = direct_biomass_model, size = n),
    color = stlmm_color("secondary"),
    alpha = 0.75,
    na.rm = TRUE
  ) +
  facet_wrap(~ county, ncol = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(series_dat$year), max(series_dat$year), by = 5)) +
  scale_size_continuous(range = c(1, 3.5), name = "Plot n") +
  labs(
    x = "Year",
    y = "Live aboveground biomass Mg/ha"
  )
