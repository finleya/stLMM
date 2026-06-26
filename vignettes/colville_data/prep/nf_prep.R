rm(list = ls())

library(tidyverse)
library(sf)

raw_shape_dir <- file.path("raw", "usda_eco_shapes")
derived_dir <- "derived"
dir.create(derived_dir, showWarnings = FALSE, recursive = TRUE)

## Read National Forest boundaries and subset to the forest of interest.
nfs <- st_read(file.path(raw_shape_dir, "S_USA.BdyAdm_LSRS_AdministrativeForest.shp"),
               quiet = TRUE) %>%
  st_make_valid()

nf <- nfs %>%
  filter(FORESTNAME == "Colville National Forest")

## Read ecoregion sections and subset to the section of interest.
ecos <- st_read(file.path(raw_shape_dir, "S_USA.ECOSYS_ECOMAPSECTIONS_2025.shp"),
                quiet = TRUE) %>%
  st_make_valid()

eco <- ecos %>%
  filter(SECTION_NA == "Northern Cascades")

## Reproject the forest boundary to match the ecoregion layer.
nf <- nf %>%
  st_transform(st_crs(eco))

## Clip the National Forest boundary to the selected ecoregion.
nf <- st_intersection(nf, eco)

## Read FIA public plot data produced by FIA_public_plot_data_prep.R.
fia_plots <- st_read(file.path(derived_dir, "FIA_public_plot_data.gpkg"),
                     layer = "fia_public_plot_data",
                     quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(st_crs(eco))

if (!"measmon" %in% names(fia_plots)) {
  fia_plots$measmon <- NA_integer_
}

if (!"measday" %in% names(fia_plots)) {
  fia_plots$measday <- NA_integer_
}

if (!"meas_date" %in% names(fia_plots)) {
  fia_plots$meas_date <- as.Date(NA)
}

if (!"meas_year_decimal" %in% names(fia_plots)) {
  fia_plots$meas_year_decimal <- coalesce(as.numeric(fia_plots$measyear),
                                          as.numeric(fia_plots$invyr))
}

## Retain FIA plots that fall within the selected ecoregion.
fia_plots <- fia_plots %>%
  st_filter(eco, .predicate = st_intersects)

## Reproject to match TCC grids.
tcc_crs <- st_crs(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

eco <- eco %>%
    st_transform(tcc_crs)

nf <- nf %>%
    st_transform(tcc_crs)

fia_plots <- fia_plots %>%
    st_transform(tcc_crs)

## Overlay the ecoregion, clipped National Forest boundary, and FIA plots.
ggplot() +
  geom_sf(data = eco, fill = "grey92", color = "grey45", linewidth = 0.4) +
  geom_sf(data = nf, fill = NA, color = "red3", linewidth = 0.8) +
  geom_sf(data = fia_plots,
          aes(color = meas_year_decimal),
          alpha = 0.75,
          size = 1.1) +
  scale_color_viridis_c(name = "Measurement year", option = "C") +
  coord_sf(datum = NA) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(color = "grey88", linewidth = 0.2),
        legend.position = "right")

## Write out for modeling.

fia_coords <- st_coordinates(fia_plots)
fia_plots$x <- fia_coords[, 1]
fia_plots$y <- fia_coords[, 2]
fia_plots$x_km <- fia_plots$x / 1000
fia_plots$y_km <- fia_plots$y / 1000

fia_plots_final <- fia_plots %>%
    select(x, y, x_km, y_km, invyr, measyear, measmon, measday,
           meas_date, meas_year_decimal,
           forest_prop,
           ag_live_tree_biomass_tons_per_acre,
           ag_live_tree_c_tons_per_acre) %>%
    transmute(x, y, x_km, y_km,
              year = as.integer(coalesce(measyear, invyr)),
              time = meas_year_decimal,
              invyr, measyear, measmon, measday, meas_date,
              forested = ifelse(forest_prop > 0, 1, 0),
              biomass_Mg_ha = ifelse(is.na(ag_live_tree_biomass_tons_per_acre), 0, 2.241702 * ag_live_tree_biomass_tons_per_acre),
              c_Mg_ha = ifelse(is.na(ag_live_tree_c_tons_per_acre), 0, 2.241702 * ag_live_tree_c_tons_per_acre))

write_csv(st_drop_geometry(fia_plots_final),
          file.path(derived_dir, "fia_plots_model_data.csv"))
st_write(fia_plots_final,
         file.path(derived_dir, "fia_plots_model_data.gpkg"),
         layer = "fia_plots_model_data",
         delete_dsn = TRUE,
         quiet = TRUE)

## Create annual prediction grids for the full ecoregion and the clipped
## National Forest. The ecoregion grid is coarser to keep vignette prediction
## runtime manageable. The National Forest gets both a 2 km vignette grid and
## a 1 km grid for finer prediction runs.
make_prediction_grid <- function(boundary, cellsize_m, years) {
  grid_points <- st_sf(
    geometry = st_make_grid(boundary, cellsize = cellsize_m, what = "centers")
  ) %>%
    st_set_crs(st_crs(boundary)) %>%
    st_filter(boundary, .predicate = st_intersects)

  grid_coords <- st_coordinates(grid_points)

  grid_points <- grid_points %>%
    mutate(x = grid_coords[, 1],
           y = grid_coords[, 2],
           x_km = x / 1000,
           y_km = y / 1000) %>%
    select(x, y, x_km, y_km)

  grid_index <- rep(seq_len(nrow(grid_points)), each = length(years))

  st_sf(
    st_drop_geometry(grid_points)[grid_index, ],
    year = rep(years, times = nrow(grid_points)),
    time = as.numeric(rep(years, times = nrow(grid_points))),
    geometry = st_geometry(grid_points)[grid_index],
    crs = st_crs(grid_points)
  )
}

prediction_grid_eco_5km <- make_prediction_grid(eco, 5000, 1999:2026)
prediction_grid_eco_10km <- make_prediction_grid(eco, 10000, 1999:2026)
prediction_grid_nf_2km <- make_prediction_grid(nf, 2000, 1999:2026)
prediction_grid_nf_1km <- make_prediction_grid(nf, 1000, 1999:2026)

write_csv(st_drop_geometry(prediction_grid_eco_5km),
          file.path(derived_dir, "prediction_grid_eco_5km_1999_2026.csv"))
st_write(prediction_grid_eco_5km,
         file.path(derived_dir, "prediction_grid_eco_5km_1999_2026.gpkg"),
         layer = "prediction_grid_eco_5km_1999_2026",
         delete_dsn = TRUE,
         quiet = TRUE)

write_csv(st_drop_geometry(prediction_grid_eco_10km),
          file.path(derived_dir, "prediction_grid_eco_10km_1999_2026.csv"))
st_write(prediction_grid_eco_10km,
         file.path(derived_dir, "prediction_grid_eco_10km_1999_2026.gpkg"),
         layer = "prediction_grid_eco_10km_1999_2026",
         delete_dsn = TRUE,
         quiet = TRUE)

write_csv(st_drop_geometry(prediction_grid_nf_2km),
          file.path(derived_dir, "prediction_grid_nf_2km_1999_2026.csv"))
st_write(prediction_grid_nf_2km,
         file.path(derived_dir, "prediction_grid_nf_2km_1999_2026.gpkg"),
         layer = "prediction_grid_nf_2km_1999_2026",
         delete_dsn = TRUE,
         quiet = TRUE)

write_csv(st_drop_geometry(prediction_grid_nf_1km),
          file.path(derived_dir, "prediction_grid_nf_1km_1999_2026.csv"))
st_write(prediction_grid_nf_1km,
         file.path(derived_dir, "prediction_grid_nf_1km_1999_2026.gpkg"),
         layer = "prediction_grid_nf_1km_1999_2026",
         delete_dsn = TRUE,
         quiet = TRUE)

## Write spatial boundaries used by the modeling vignette.
st_write(nf,
         file.path(derived_dir, "colville_national_forest_northern_cascades.gpkg"),
         layer = "colville_national_forest_northern_cascades",
         delete_dsn = TRUE,
         quiet = TRUE)

st_write(eco,
         file.path(derived_dir, "northern_cascades_ecoregion.gpkg"),
         layer = "northern_cascades_ecoregion",
         delete_dsn = TRUE,
         quiet = TRUE)
