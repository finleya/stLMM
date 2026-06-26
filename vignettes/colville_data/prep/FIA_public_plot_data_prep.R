## Acknowledgment: This script is based primarily on code and workflow
## developed by Brian Walters for preparing FIA public plot data.

library(tidyverse)
library(rFIA)
library(sf)

fia_states <- "WA"
raw_fia_dir <- file.path("raw", "fia")
derived_dir <- "derived"
dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

decimal_year <- function(year, month, day) {
  date <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", year, month, day)))
  year_start <- suppressWarnings(as.Date(sprintf("%04d-01-01", year)))
  next_year_start <- suppressWarnings(as.Date(sprintf("%04d-01-01", year + 1)))
  out <- year + as.numeric(date - year_start) / as.numeric(next_year_start - year_start)
  ifelse(is.na(date) | is.na(year_start) | is.na(next_year_start), NA_real_, out)
}

## Read the FIA tables needed for plot, condition, tree, and survey attributes.
dat <- readFIA(dir = raw_fia_dir, states = fia_states,
               tables = c("COND", "PLOT", "TREE", "SURVEY"))

## Restrict the data to annual-design plots and remove ozone-only surveys.
dat$SURVEY <- dat$SURVEY %>%
  rename_with(tolower) %>%
  filter(ann_inventory == "Y" & p3_ozone_ind == "N")

dat$PLOT <- dat$PLOT %>%
  rename_with(tolower) %>%
  inner_join(dat$SURVEY %>% select(srv_cn = cn), by = join_by(srv_cn))

dat$COND <- dat$COND %>%
  rename_with(tolower) %>%
  inner_join(dat$PLOT %>% select(plt_cn = cn), by = join_by(plt_cn))

dat$TREE <- dat$TREE %>%
  rename_with(tolower) %>%
  inner_join(dat$PLOT %>% select(plt_cn = cn), by = join_by(plt_cn)) %>%
  filter(between(subp, 1, 4)) %>% # Retain trees measured on FIA subplots 1 through 4.
  filter(!is.na(tpa_unadj)) # Retain trees with an expansion factor.


## Build a plot-level table with location, inventory timing, condition summaries,
## and live/dead aboveground tree carbon and biomass.
plot_data <- dat$PLOT %>%
  select(plt_cn = cn, statecd, unitcd, countycd, plot,
         invyr, measyear, measmon, measday, intensity, lat, lon) %>%
  mutate(
    meas_date = suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", measyear, measmon, measday))),
    meas_year_decimal = coalesce(decimal_year(measyear, measmon, measday), as.numeric(measyear), as.numeric(invyr))
  ) %>%
  ## Summarize forested conditions, forest proportion, and non-sampled proportion.
  inner_join(dat$COND %>%
               summarise(forest_cond_count = sum(ifelse(cond_status_cd == 1, 1, 0)),
                         forest_prop = sum(ifelse(cond_status_cd == 1, condprop_unadj, 0)),
                         nonsamp_prop = sum(ifelse(cond_status_cd == 5, condprop_unadj, 0)),
                         .by = plt_cn),
             by = join_by(plt_cn)) %>%
  ## Add live aboveground tree carbon and biomass on forested conditions.
  left_join(dat$TREE %>%
              filter(statuscd == 1) %>% # Live trees only.
              inner_join(dat$COND %>% filter(cond_status_cd == 1) %>% select(plt_cn, condid),
                         by = join_by(plt_cn, condid)) %>% # Forested conditions only.
              summarise(ag_live_tree_c_tons_per_acre = sum(tpa_unadj * coalesce(carbon_ag / 2000, 0)),
                        ag_live_tree_biomass_tons_per_acre = sum(tpa_unadj * coalesce(drybio_ag / 2000, 0)),
                        .by = plt_cn),
            by = join_by(plt_cn)) %>%
  ## Add standing dead aboveground tree carbon and biomass on forested conditions.
  left_join(dat$TREE %>%
              ## Standing dead saplings (DBH 1.0-4.9 in) were not measured in all FIA years.
              ## For consistent comparisons across years, consider filtering to dia >= 5.0.
              filter(statuscd == 2 & standing_dead_cd == 1) %>%
              inner_join(dat$COND %>% filter(cond_status_cd == 1) %>% select(plt_cn, condid),
                         by = join_by(plt_cn, condid)) %>% # Forested conditions only.
              summarise(ag_dead_tree_c_tons_per_acre = sum(tpa_unadj * coalesce(carbon_ag / 2000, 0)),
                        ag_dead_tree_biomass_tons_per_acre = sum(tpa_unadj * coalesce(drybio_ag / 2000, 0)),
                        .by = plt_cn),
            by = join_by(plt_cn)) %>%
  ## Convert missing carbon and biomass values to true zeroes for forested plots.
  mutate(across(ag_live_tree_c_tons_per_acre:ag_dead_tree_biomass_tons_per_acre,
                ~ ifelse(forest_cond_count > 0, coalesce(.x, 0), .x)))


## TODO: Before finalizing the vignette, remove columns that should not be shared
## and rename columns for reader-facing clarity.
plot_data_for_export <- plot_data

## FIA public plot coordinates are latitude/longitude in NAD83 (EPSG:4269).
plot_data_sf <- plot_data_for_export %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4269, remove = FALSE)

## Write the public plot-level table for use in the vignette.
write_csv(plot_data_for_export,
          file.path(derived_dir, "FIA_public_plot_data.csv"))

## Write the same plot data as a spatial GeoPackage.
st_write(plot_data_sf,
         file.path(derived_dir, "FIA_public_plot_data.gpkg"),
         layer = "fia_public_plot_data",
         delete_dsn = TRUE,
         quiet = TRUE)
