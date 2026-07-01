library(sf)
library(tidyverse)
library(knitr)
library(kableExtra)

source("article-utils.R")

set.seed(1)
theme_set(theme_bw(base_size = 12))

data_dir <- "wa_data"

wa_counties <- read_rds(file.path(data_dir, "wa_counties.rds"))
unit_plots <- read_csv(file.path(data_dir, "wa_unit_plots.csv"), show_col_types = FALSE)
direct_estimates <- read_csv(file.path(data_dir, "wa_direct_estimates.csv"), show_col_types = FALSE)
grid_5km <- read_csv(file.path(data_dir, "wa_unit_prediction_grid_5km.csv"), show_col_types = FALSE)
grid_10km <- read_csv(file.path(data_dir, "wa_unit_prediction_grid_10km.csv"), show_col_types = FALSE)
manifest <- read_csv(file.path(data_dir, "DATA_MANIFEST.csv"), show_col_types = FALSE)

unit_plots <- unit_plots |>
  mutate(
    forested = as.logical(forested),
    county_fips = as.character(county_fips)
  )

direct_estimates <- direct_estimates |>
  mutate(county_fips = as.character(county_fips))

grid_5km <- grid_5km |>
  mutate(county_fips = as.character(county_fips))


data_pack_summary <- tibble(
  quantity = c(
    "counties",
    "FIA units",
    "direct-estimate years",
    "unit plot rows",
    "direct-estimate rows",
    "model-ready direct-estimate rows",
    "single-plot rows retained but excluded",
    "direct-estimate rows using variance floor",
    "5 km prediction-grid rows",
    "10 km prediction-grid rows",
    "plot rows dropped for missing TCC"
  ),
  value = c(
    nrow(wa_counties),
    n_distinct(unit_plots$fia_unit),
    n_distinct(direct_estimates$year),
    nrow(unit_plots),
    nrow(direct_estimates),
    manifest_value("direct_estimate_model_rows"),
    manifest_value("direct_estimate_single_plot_excluded_rows"),
    manifest_value("direct_estimate_variance_floor_rows"),
    manifest_value("prediction_grid_5km_full_rows"),
    nrow(grid_10km),
    manifest_value("unit_plot_dropped_missing_tcc_rows")
  )
) |>
  mutate(value = fmt_int(as.numeric(value)))

show_table(data_pack_summary, caption = "Washington county biomass article data bundle.")

plot_status_summary <- unit_plots |>
  mutate(status = if_else(forested, "forested", "nonforested")) |>
  count(status, name = "plot_rows") |>
  mutate(percent = 100 * plot_rows / sum(plot_rows))

fia_unit_summary <- unit_plots |>
  mutate(status = if_else(forested, "forested", "nonforested")) |>
  count(fia_unit, fia_unit_name, status, name = "plot_rows") |>
  pivot_wider(names_from = status, values_from = plot_rows, values_fill = 0) |>
  mutate(
    total = forested + nonforested,
    forested_percent = 100 * forested / total
  ) |>
  arrange(fia_unit)

show_table(
  plot_status_summary |>
    mutate(
      plot_rows = fmt_int(plot_rows),
      percent = fmt_num(percent, 1)
    ),
  caption = "Forested and nonforested FIA plot rows."
)

show_table(
  fia_unit_summary |>
    mutate(
      across(c(forested, nonforested, total), fmt_int),
      forested_percent = fmt_num(forested_percent, 1)
    ),
  caption = "FIA plot rows by Washington FIA unit."
)

direct_status_summary <- direct_estimates |>
  count(direct_estimate_status, name = "county_years") |>
  mutate(
    status = recode(
      direct_estimate_status,
      modeled = "model-ready direct estimate",
      no_plots = "no FIA plots",
      single_plot_no_variance = "single plot; no design variance",
      missing_model_variance = "missing model variance",
      missing = "missing"
    ),
    percent = 100 * county_years / sum(county_years)
  ) |>
  select(status, county_years, percent)

show_table(
  direct_status_summary |>
    mutate(
      county_years = fmt_int(county_years),
      percent = fmt_num(percent, 1)
    ),
  caption = "County-year direct-estimate status."
)

direct_year_summary <- direct_estimates |>
  group_by(year) |>
  summarise(
    counties_with_plots = sum(n > 0),
    counties_with_model_direct_estimates = sum(direct_estimate_in_model),
    single_plot_county_years = sum(direct_estimate_status == "single_plot_no_variance"),
    median_plot_n = median(n[n > 0], na.rm = TRUE),
    median_direct_biomass = median(direct_biomass_model, na.rm = TRUE),
    median_direct_se = median(direct_biomass_se_model, na.rm = TRUE),
    .groups = "drop"
  )

show_table(
  direct_year_summary |>
    filter(year %in% c(2016, 2018, 2020, 2022, 2024, 2025)) |>
    mutate(
      counties_with_model_direct_estimates = fmt_int(counties_with_model_direct_estimates),
      single_plot_county_years = fmt_int(single_plot_county_years),
      median_plot_n = fmt_num(median_plot_n, 0),
      median_direct_biomass = fmt_num(median_direct_biomass, 1),
      median_direct_se = fmt_num(median_direct_se, 1)
    ),
  caption = "Selected-year county direct-estimate characteristics."
)

panel_years <- c(2016, 2018, 2020, 2022, 2024)

direct_map_dat <- wa_counties |>
  left_join(
    direct_estimates |> filter(year %in% panel_years),
    by = "county_fips"
  )

ggplot(direct_map_dat) +
  geom_sf(aes(fill = direct_biomass_model), color = "grey80", linewidth = 0.12) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "Mg/ha",
    na.value = "grey92"
  ) +
  labs(fill = "Mg/ha") +
  theme(plot.margin = margin(0, 0, 0, 0))

plot_map_dat <- unit_plots |>
  filter(year %in% panel_years)

ggplot() +
  geom_sf(data = wa_counties, fill = "grey96", color = "grey70", linewidth = 0.12) +
  geom_point(
    data = plot_map_dat,
    aes(x = x, y = y, color = agb_live_mg_ha),
    alpha = 0.65,
    size = 0.55
  ) +
  coord_sf(crs = st_crs(wa_counties), expand = FALSE) +
  facet_wrap(~ year, ncol = 3) +
  scale_color_gradientn(
    colors = stlmm_palette(),
    name = "Mg/ha",
    trans = "sqrt"
  ) +
  theme(plot.margin = margin(0, 0, 0, 0))

tcc_map_dat <- wa_counties |>
  left_join(
    direct_estimates |> filter(year %in% panel_years),
    by = "county_fips"
  )

ggplot(tcc_map_dat) +
  geom_sf(aes(fill = county_mean_tcc), color = "grey80", linewidth = 0.12) +
  coord_sf(expand = FALSE) +
  facet_wrap(~ year, ncol = 3) +
  scale_fill_gradientn(
    colors = stlmm_palette(),
    name = "% TCC",
    limits = c(0, 100),
    na.value = "grey92"
  ) +
  theme(plot.margin = margin(0, 0, 0, 0))

grid_plot_dat <- grid_5km |>
  filter(year %in% panel_years)

ggplot() +
  geom_sf(data = wa_counties, fill = "grey96", color = "grey75", linewidth = 0.12) +
  geom_point(
    data = grid_plot_dat,
    aes(x = x, y = y, color = tcc_mean),
    size = 0.18,
    alpha = 0.75
  ) +
  coord_sf(crs = st_crs(wa_counties), expand = FALSE) +
  facet_wrap(~ year, ncol = 3) +
  scale_color_gradientn(
    colors = stlmm_palette(),
    name = "% TCC",
    limits = c(0, 100)
  ) +
  theme(plot.margin = margin(0, 0, 0, 0))
