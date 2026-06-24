## Build Washington county biomass vignette data from raw FIA tables.
##
## This script is intentionally written as a readable, auditable workflow.  It
## starts from FIA CSV tables, prepares plot-level biomass responses, builds a
## Washington county support from stLMM's county polygons, creates annual
## prediction grids, attaches annual TCC covariates, and writes the small files
## expected by the WA county vignette series.
##
## The TCC predictor follows the PNW vignette convention: annual percent tree
## canopy cover is averaged over a square block with about a one-mile half-width
## support, then sampled at plot and grid locations.  Cached PNW block rasters
## can be reused to save time; set WA_TCC_REUSE_PNW_BLOCKS=false to rebuild the
## WA-clipped block rasters from the raw TCC rasters.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(rFIA)
  library(sf)
  library(terra)
  library(tidyr)
})

## ---- configuration ---------------------------------------------------------

env_path <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_int <- function(name, default) {
  as.integer(env_path(name, as.character(default)))
}

env_num <- function(name, default) {
  as.numeric(env_path(name, as.character(default)))
}

env_bool <- function(name, default = FALSE) {
  value <- tolower(env_path(name, if (default) "true" else "false"))
  value %in% c("1", "true", "t", "yes", "y")
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else
  file.path(getwd(), "scripts", "build_wa_vignette_data.R")
script_dir <- dirname(normalizePath(script_file, mustWork = FALSE))
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

input_fia_dir <- env_path("WA_FIA_DATA_DIR", file.path(project_dir, "..", "data"))
output_dir <- env_path("WA_VIGNETTE_DATA_DIR", file.path(project_dir, "wa_data"))
cache_dir <- env_path("WA_PREP_CACHE_DIR", file.path(project_dir, "cache"))

stunitco_file <- env_path(
  "WA_STUNITCO_GPKG",
  "/home/andy/Rdevel/stLMM-public/inst/extdata/stunitco.gpkg"
)

tcc_input_dir <- env_path("WA_TCC_INPUT_DIR", "/mnt/disk1/TCC/tcc_conus_v2025-6-tifs")
tcc_output_dir <- env_path("WA_TCC_OUTPUT_DIR", "/mnt/disk2/wa_vignettes")
pnw_block_dir <- env_path("WA_PNW_TCC_BLOCK_DIR", "/mnt/disk2/pnw/tcc_mean_1mi_block_pnw")
reuse_pnw_blocks <- env_bool("WA_TCC_REUSE_PNW_BLOCKS", TRUE)
clip_reused_pnw_blocks <- env_bool("WA_TCC_CLIP_REUSED_PNW_BLOCKS", TRUE)

prediction_years <- seq(env_int("WA_START_YEAR", 1999), env_int("WA_END_YEAR", 2025))
grid_resolutions <- c("1km" = 1000, "2km" = 2000, "5km" = 5000, "10km" = 10000)
grid_cells_per_county_year <- env_int("WA_GRID_CELLS_PER_COUNTY_YEAR", 400)
grid_sample_seed <- env_int("WA_GRID_SAMPLE_SEED", 23)
window_radius_m <- env_num("WA_TCC_WINDOW_RADIUS_M", 1609.344)
force <- env_bool("WA_FORCE", FALSE)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tcc_output_dir, "tcc_clipped_wa"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tcc_output_dir, "tcc_mean_1mi_block_wa"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tcc_output_dir, "terra_tmp"), recursive = TRUE, showWarnings = FALSE)

set.seed(grid_sample_seed)
Sys.setenv(GDAL_NUM_THREADS = Sys.getenv("GDAL_NUM_THREADS", "ALL_CPUS"))
terraOptions(tempdir = file.path(tcc_output_dir, "terra_tmp"), memfrac = 0.75, threads = TRUE)

log_msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
  flush.console()
}

write_csv_logged <- function(x, path) {
  log_msg("Writing", path)
  write_csv(x, path)
}

write_rds_logged <- function(x, path) {
  log_msg("Writing", path)
  write_rds(x, path)
}

## ---- helpers ---------------------------------------------------------------

decimal_year <- function(year, month, day) {
  date <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", year, month, day)))
  year_start <- suppressWarnings(as.Date(sprintf("%04d-01-01", year)))
  next_year_start <- suppressWarnings(as.Date(sprintf("%04d-01-01", year + 1)))
  out <- year + as.numeric(date - year_start) / as.numeric(next_year_start - year_start)
  ifelse(is.na(date) | is.na(year_start) | is.na(next_year_start), NA_real_, out)
}

mg_ha_from_tons_acre <- function(x) {
  2.241702 * x
}

albers_conus_crs <- st_crs(
  "+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
)

## ---- Washington county support --------------------------------------------

wa_counties_path <- file.path(output_dir, "wa_counties.rds")

build_wa_counties <- function() {
  log_msg("Reading stLMM county polygons:", stunitco_file)
  counties <- st_read(stunitco_file, quiet = TRUE) |>
    st_make_valid() |>
    st_transform(albers_conus_crs)

  wa <- counties |>
    filter(STATENM == "Washington") |>
    group_by(COUNTYFIPS, COUNTYNM, STATECD, STATENM, UNITCD, UNITNM) |>
    summarise(.groups = "drop") |>
    mutate(
      county_fips = as.character(COUNTYFIPS),
      county = as.character(COUNTYNM),
      fia_unit = as.integer(UNITCD),
      fia_unit_name = as.character(UNITNM)
    )

  names(wa)[names(wa) == attr(wa, "sf_column")] <- "geometry"
  st_geometry(wa) <- "geometry"

  wa |>
    mutate(county_area_ha = as.numeric(st_area(geometry)) / 10000) |>
    select(county_fips, county, fia_unit, fia_unit_name, county_area_ha, geometry) |>
    arrange(county_fips)
}

if (file.exists(wa_counties_path) && !force) {
  wa_counties <- read_rds(wa_counties_path)
} else {
  wa_counties <- build_wa_counties()
  write_rds_logged(wa_counties, wa_counties_path)
}

wa_boundary <- wa_counties |>
  summarise(geometry = st_union(geometry), .groups = "drop")

county_lookup <- wa_counties |>
  st_drop_geometry() |>
  select(county_fips, county, fia_unit, fia_unit_name, county_area_ha)

## ---- raw FIA plot table ----------------------------------------------------

public_plot_csv <- file.path(cache_dir, "wa_fia_public_plot_data.csv")
public_plot_gpkg <- file.path(cache_dir, "wa_fia_public_plot_data.gpkg")

build_public_plot_data <- function() {
  log_msg("Reading raw FIA tables from:", input_fia_dir)
  dat <- readFIA(
    dir = input_fia_dir,
    states = "WA",
    tables = c("COND", "PLOT", "TREE", "SURVEY")
  )

  log_msg("Restricting to annual-design, non-ozone surveys.")
  dat$SURVEY <- dat$SURVEY |>
    rename_with(tolower) |>
    filter(ann_inventory == "Y", p3_ozone_ind == "N")

  dat$PLOT <- dat$PLOT |>
    rename_with(tolower) |>
    inner_join(dat$SURVEY |> select(srv_cn = cn), by = join_by(srv_cn))

  dat$COND <- dat$COND |>
    rename_with(tolower) |>
    inner_join(dat$PLOT |> select(plt_cn = cn), by = join_by(plt_cn))

  dat$TREE <- dat$TREE |>
    rename_with(tolower) |>
    inner_join(dat$PLOT |> select(plt_cn = cn), by = join_by(plt_cn)) |>
    filter(between(subp, 1, 4), !is.na(tpa_unadj))

  log_msg("Building plot-level condition and tree biomass summaries.")
  plot_data <- dat$PLOT |>
    select(
      plt_cn = cn, statecd, unitcd, countycd, plot,
      invyr, measyear, measmon, measday, intensity, lat, lon
    ) |>
    mutate(
      meas_date = suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", measyear, measmon, measday))),
      meas_year_decimal = coalesce(
        decimal_year(measyear, measmon, measday),
        as.numeric(measyear),
        as.numeric(invyr)
      )
    ) |>
    inner_join(
      dat$COND |>
        summarise(
          forest_cond_count = sum(ifelse(cond_status_cd == 1, 1, 0)),
          forest_prop = sum(ifelse(cond_status_cd == 1, condprop_unadj, 0)),
          nonsamp_prop = sum(ifelse(cond_status_cd == 5, condprop_unadj, 0)),
          .by = plt_cn
        ),
      by = join_by(plt_cn)
    ) |>
    left_join(
      dat$TREE |>
        filter(statuscd == 1) |>
        inner_join(
          dat$COND |> filter(cond_status_cd == 1) |> select(plt_cn, condid),
          by = join_by(plt_cn, condid)
        ) |>
        summarise(
          ag_live_tree_c_tons_per_acre = sum(tpa_unadj * coalesce(carbon_ag / 2000, 0)),
          ag_live_tree_biomass_tons_per_acre = sum(tpa_unadj * coalesce(drybio_ag / 2000, 0)),
          .by = plt_cn
        ),
      by = join_by(plt_cn)
    ) |>
    left_join(
      dat$TREE |>
        filter(statuscd == 2, standing_dead_cd == 1) |>
        inner_join(
          dat$COND |> filter(cond_status_cd == 1) |> select(plt_cn, condid),
          by = join_by(plt_cn, condid)
        ) |>
        summarise(
          ag_dead_tree_c_tons_per_acre = sum(tpa_unadj * coalesce(carbon_ag / 2000, 0)),
          ag_dead_tree_biomass_tons_per_acre = sum(tpa_unadj * coalesce(drybio_ag / 2000, 0)),
          .by = plt_cn
        ),
      by = join_by(plt_cn)
    ) |>
    mutate(across(
      ag_live_tree_c_tons_per_acre:ag_dead_tree_biomass_tons_per_acre,
      ~ ifelse(forest_cond_count > 0, coalesce(.x, 0), .x)
    ))

  plot_sf <- plot_data |>
    filter(!is.na(lon), !is.na(lat)) |>
    st_as_sf(coords = c("lon", "lat"), crs = 4269, remove = FALSE)

  write_csv_logged(plot_data, public_plot_csv)
  log_msg("Writing", public_plot_gpkg)
  st_write(plot_sf, public_plot_gpkg, layer = "wa_fia_public_plot_data",
           delete_dsn = TRUE, quiet = TRUE)

  plot_sf
}

if (file.exists(public_plot_gpkg) && !force) {
  public_plot_sf <- st_read(public_plot_gpkg, quiet = TRUE)
} else {
  public_plot_sf <- build_public_plot_data()
}

## ---- model plot table ------------------------------------------------------

plot_model_csv <- file.path(cache_dir, "wa_fia_plots_model_data.csv")
plot_model_gpkg <- file.path(cache_dir, "wa_fia_plots_model_data.gpkg")

build_model_plot_data <- function(public_plot_sf) {
  log_msg("Projecting FIA plots and retaining plots inside Washington counties.")
  plots <- public_plot_sf |>
    st_make_valid() |>
    st_transform(st_crs(wa_counties)) |>
    st_filter(wa_boundary, .predicate = st_intersects)

  county_idx <- st_intersects(plots, wa_counties)
  has_county <- lengths(county_idx) > 0L
  if (any(!has_county)) {
    log_msg("Dropping", sum(!has_county), "plots without a county intersection.")
  }

  plots <- plots[has_county, ]
  county_rows <- vapply(county_idx[has_county], `[`, integer(1), 1L)
  plot_county <- wa_counties[county_rows, ] |>
    st_drop_geometry() |>
    select(county_fips, county, fia_unit, fia_unit_name)

  coords <- st_coordinates(plots)

  plot_model <- plots |>
    st_drop_geometry() |>
    bind_cols(plot_county) |>
    mutate(
      x = coords[, 1],
      y = coords[, 2],
      x_km = x / 1000,
      y_km = y / 1000,
      year = as.integer(coalesce(measyear, invyr)),
      time = meas_year_decimal,
      forested = forest_prop > 0,
      proportion_forested = forest_prop,
      agb_live_mg_ha = ifelse(
        is.na(ag_live_tree_biomass_tons_per_acre),
        0,
        mg_ha_from_tons_acre(ag_live_tree_biomass_tons_per_acre)
      ),
      agb_live_c_mg_ha = ifelse(
        is.na(ag_live_tree_c_tons_per_acre),
        0,
        mg_ha_from_tons_acre(ag_live_tree_c_tons_per_acre)
      ),
      agb_dead_mg_ha = ifelse(
        is.na(ag_dead_tree_biomass_tons_per_acre),
        0,
        mg_ha_from_tons_acre(ag_dead_tree_biomass_tons_per_acre)
      )
    ) |>
    select(
      plt_cn, statecd, unitcd, countycd, plot,
      x, y, x_km, y_km, year, time,
      invyr, measyear, measmon, measday, meas_date,
      county_fips, county, fia_unit, fia_unit_name,
      forested, proportion_forested, nonsamp_prop,
      agb_live_mg_ha, agb_live_c_mg_ha, agb_dead_mg_ha
    ) |>
    filter(year %in% prediction_years) |>
    arrange(county_fips, year, plt_cn)

  plot_model_sf <- st_as_sf(plot_model, coords = c("x", "y"), crs = st_crs(wa_counties), remove = FALSE)

  write_csv_logged(plot_model, plot_model_csv)
  log_msg("Writing", plot_model_gpkg)
  st_write(plot_model_sf, plot_model_gpkg, layer = "wa_fia_plots_model_data",
           delete_dsn = TRUE, quiet = TRUE)

  plot_model
}

if (file.exists(plot_model_csv) && !force) {
  plot_model <- read_csv(plot_model_csv, show_col_types = FALSE)
} else {
  plot_model <- build_model_plot_data(public_plot_sf)
}

## ---- annual prediction grids ----------------------------------------------

grid_cache_path <- function(resolution) {
  file.path(cache_dir, sprintf("wa_prediction_grid_%s_%s_%s.csv",
                              resolution, min(prediction_years), max(prediction_years)))
}

make_prediction_grid <- function(cellsize_m, years) {
  grid_points <- st_sf(
    geometry = st_make_grid(wa_boundary, cellsize = cellsize_m, what = "centers")
  ) |>
    st_set_crs(st_crs(wa_counties)) |>
    st_filter(wa_boundary, .predicate = st_intersects)

  county_idx <- st_intersects(grid_points, wa_counties)
  has_county <- lengths(county_idx) > 0L
  grid_points <- grid_points[has_county, ]
  county_rows <- vapply(county_idx[has_county], `[`, integer(1), 1L)

  coords <- st_coordinates(grid_points)
  base_grid <- grid_points |>
    st_drop_geometry() |>
    mutate(x = coords[, 1], y = coords[, 2], x_km = x / 1000, y_km = y / 1000) |>
    bind_cols(
      wa_counties[county_rows, ] |>
        st_drop_geometry() |>
        select(county_fips, county, fia_unit, fia_unit_name)
    )

  row_index <- rep(seq_len(nrow(base_grid)), each = length(years))
  base_grid[row_index, ] |>
    mutate(year = rep(years, times = nrow(base_grid))) |>
    select(x, y, x_km, y_km, year, county_fips, county, fia_unit, fia_unit_name) |>
    arrange(county_fips, year, x, y)
}

grid_tables <- list()
for (resolution in names(grid_resolutions)) {
  grid_file <- grid_cache_path(resolution)
  if (file.exists(grid_file) && !force) {
    grid_tables[[resolution]] <- fread(grid_file)
  } else {
    log_msg("Building", resolution, "Washington prediction grid.")
    grid_tables[[resolution]] <- make_prediction_grid(grid_resolutions[[resolution]], prediction_years)
    fwrite(grid_tables[[resolution]], grid_file)
  }
}

## ---- TCC raster catalog and sampling --------------------------------------

tcc_files <- list.files(tcc_input_dir, pattern = "\\.tif$", full.names = TRUE)
tcc_year <- as.integer(sub(".*_([0-9]{4})0101_[0-9]{4}1231\\.tif$", "\\1", basename(tcc_files)))
tcc_catalog <- data.table(year = tcc_year, file = tcc_files)[!is.na(year)]
setorder(tcc_catalog, year)

missing_tcc_years <- setdiff(prediction_years, tcc_catalog$year)
if (length(missing_tcc_years)) {
  log_msg("TCC unavailable for years:", paste(missing_tcc_years, collapse = ", "))
}

template_raster <- rast(tcc_catalog$file[1])
wa_boundary_vect <- vect(st_transform(wa_boundary, crs(template_raster)))

kernel_radius_cells <- ceiling(window_radius_m / mean(res(template_raster)))
kernel_size <- 2 * kernel_radius_cells + 1
kernel <- matrix(1, nrow = kernel_size, ncol = kernel_size)

log_msg(
  "Using", window_radius_m, "m half-width TCC support;",
  kernel_size, "x", kernel_size, "pixel block."
)

wa_clipped_tcc_path <- function(year) {
  file.path(tcc_output_dir, "tcc_clipped_wa", sprintf("tcc_wa_%d.tif", year))
}

wa_block_tcc_path <- function(year) {
  file.path(tcc_output_dir, "tcc_mean_1mi_block_wa", sprintf("tcc_wa_mean_1mi_block_%d.tif", year))
}

pnw_block_tcc_path <- function(year) {
  file.path(pnw_block_dir, sprintf("tcc_pnw_mean_1mi_block_%d.tif", year))
}

tcc_block_cache <- new.env(parent = emptyenv())

make_wa_clipped_tcc <- function(year) {
  out <- wa_clipped_tcc_path(year)
  if (file.exists(out) && !force) return(out)

  src <- tcc_catalog[["file"]][tcc_catalog[["year"]] == year]
  if (length(src) != 1L) {
    stop("Expected exactly one TCC raster for year ", year, ".")
  }

  log_msg("Clipping/masking raw TCC for Washington:", year)
  r <- rast(src)
  NAflag(r) <- 255
  r <- crop(r, wa_boundary_vect, snap = "out")
  r <- mask(r, wa_boundary_vect)
  r <- classify(r, rbind(c(101, 255, NA)), include.lowest = TRUE)
  names(r) <- "tcc"
  writeRaster(
    r, out, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES",
                         "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS"))
  )
  out
}

make_wa_block_tcc <- function(year) {
  out <- wa_block_tcc_path(year)
  if (file.exists(out) && !force) return(out)

  log_msg("Computing WA one-mile half-width block mean TCC:", year)
  clipped <- make_wa_clipped_tcc(year)
  r <- rast(clipped)
  block_r <- aggregate(r, fact = kernel_size, fun = "mean", na.rm = TRUE)
  names(block_r) <- "tcc_mean_1mi"
  writeRaster(
    block_r, out, overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES",
                         "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS"))
  )
  out
}

tcc_block_for_year <- function(year) {
  cache_key <- as.character(year)
  if (exists(cache_key, envir = tcc_block_cache, inherits = FALSE)) {
    return(get(cache_key, envir = tcc_block_cache, inherits = FALSE))
  }

  cached <- pnw_block_tcc_path(year)
  if (reuse_pnw_blocks && file.exists(cached)) {
    if (!clip_reused_pnw_blocks) {
      assign(cache_key, cached, envir = tcc_block_cache)
      return(cached)
    }

    out <- wa_block_tcc_path(year)
    if (file.exists(out) && !force) {
      return(out)
    }

    log_msg("Clipping cached PNW block-mean TCC to Washington:", year)
    r <- rast(cached)
    wa_vect <- vect(st_transform(wa_boundary, crs(r)))
    r <- crop(r, wa_vect, snap = "out")
    r <- mask(r, wa_vect)
    names(r) <- "tcc_mean_1mi"
    writeRaster(
      r, out, overwrite = TRUE,
      wopt = list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES",
                           "BIGTIFF=IF_SAFER", "NUM_THREADS=ALL_CPUS"))
    )
    assign(cache_key, out, envir = tcc_block_cache)
    return(out)
  }

  out <- make_wa_block_tcc(year)
  assign(cache_key, out, envir = tcc_block_cache)
  out
}

attach_tcc <- function(dt, label) {
  dt <- as.data.table(dt)
  dt[, tcc_mean := NA_real_]

  years_to_sample <- sort(intersect(unique(dt$year), tcc_catalog$year))
  for (yr in years_to_sample) {
    log_msg("Sampling TCC for", label, "year", yr)
    idx <- which(dt$year == yr)
    r <- rast(tcc_block_for_year(yr))
    vals <- terra::extract(r, as.matrix(dt[idx, .(x, y)]))[[1]]
    dt[idx, tcc_mean := vals]
  }

  as_tibble(dt)
}

## ---- final unit plot file --------------------------------------------------

unit_plots_file <- file.path(output_dir, "wa_unit_plots.csv")

if (file.exists(unit_plots_file) && !force) {
  unit_plots <- read_csv(unit_plots_file, show_col_types = FALSE)
  unit_plot_dropped_missing_tcc_rows <- NA_integer_
} else {
  unit_plots_with_tcc <- plot_model |>
    attach_tcc("FIA plots") |>
    mutate(plot_id = row_number()) |>
    select(
      plot_id, plt_cn, x, y, year, time,
      county_fips, county, fia_unit, fia_unit_name,
      forested, proportion_forested, nonsamp_prop,
      agb_live_mg_ha, agb_live_c_mg_ha, agb_dead_mg_ha,
      tcc_mean
    ) |>
    arrange(county_fips, year, plot_id)

  unit_plot_dropped_missing_tcc_rows <- sum(is.na(unit_plots_with_tcc$tcc_mean))
  if (unit_plot_dropped_missing_tcc_rows > 0L) {
    log_msg(
      "Dropping", unit_plot_dropped_missing_tcc_rows,
      "FIA plot rows with missing plot-level TCC."
    )
  }

  unit_plots <- unit_plots_with_tcc |>
    filter(!is.na(tcc_mean))

  write_csv_logged(unit_plots, unit_plots_file)
}

## ---- final prediction grid files ------------------------------------------

grid_outputs <- list()
grid_sample_outputs <- list()
grid_full_rows <- setNames(integer(length(grid_tables)), names(grid_tables))
grid_sample_rows <- setNames(integer(length(grid_tables)), names(grid_tables))

count_csv_rows <- function(path) {
  if (!file.exists(path)) {
    return(NA_integer_)
  }
  nrow(fread(path, select = 1L, showProgress = FALSE))
}

for (resolution in names(grid_tables)) {
  output_file <- file.path(output_dir, sprintf("wa_unit_prediction_grid_%s.csv", resolution))
  sample_output_file <- file.path(output_dir, sprintf("wa_unit_prediction_grid_%s_sample.csv", resolution))

  if (file.exists(output_file) && file.exists(sample_output_file) && !force) {
    grid_full_rows[[resolution]] <- count_csv_rows(output_file)
    grid_sample_rows[[resolution]] <- count_csv_rows(sample_output_file)
    grid_outputs[[resolution]] <- NULL
    grid_sample_outputs[[resolution]] <- NULL
    next
  }

  grid_with_tcc <- grid_tables[[resolution]] |>
    attach_tcc(paste(resolution, "prediction grid")) |>
    filter(!is.na(tcc_mean))

  write_csv_logged(grid_with_tcc, output_file)
  grid_full_rows[[resolution]] <- nrow(grid_with_tcc)

  if (!is.na(grid_cells_per_county_year) && grid_cells_per_county_year > 0) {
    log_msg(
      "Sampling", resolution, "grid to at most",
      grid_cells_per_county_year, "cells per county-year."
    )
    grid_sample <- grid_with_tcc |>
      group_by(county_fips, year) |>
      group_modify(function(.x, .y) {
        if (nrow(.x) > grid_cells_per_county_year) {
          slice_sample(.x, n = grid_cells_per_county_year)
        } else {
          .x
        }
      }) |>
      ungroup() |>
      arrange(county_fips, year, x, y)
  } else {
    grid_sample <- grid_with_tcc
  }

  write_csv_logged(grid_sample, sample_output_file)
  grid_sample_rows[[resolution]] <- nrow(grid_sample)
  grid_outputs[[resolution]] <- grid_with_tcc
  grid_sample_outputs[[resolution]] <- grid_sample
}

## ---- direct county-year estimates -----------------------------------------

direct_estimates_file <- file.path(output_dir, "wa_direct_estimates.csv")

summarize_county_tcc <- function(years) {
  county_tcc_rows <- lapply(years, function(yr) {
    if (!yr %in% tcc_catalog$year) {
      return(county_lookup |>
               mutate(year = yr, county_mean_tcc = NA_real_, n_tcc_pixels = NA_integer_))
    }

    log_msg("Computing county zonal mean TCC from full raster for year", yr)
    r <- rast(tcc_block_for_year(yr))
    counties_raster_crs <- st_transform(wa_counties, crs(r))

    stats <- exactextractr::exact_extract(
      r,
      counties_raster_crs,
      function(values, coverage_fraction) {
        valid <- is.finite(values) & coverage_fraction > 0
        data.frame(
          county_mean_tcc = if (any(valid)) {
            weighted.mean(values[valid], coverage_fraction[valid])
          } else {
            NA_real_
          },
          n_tcc_pixels = sum(valid)
        )
      }
    )

    bind_cols(
      county_lookup |> select(county_fips, county, fia_unit, fia_unit_name),
      year = yr,
      bind_rows(stats)
    )
  })

  bind_rows(county_tcc_rows) |>
    arrange(county_fips, year)
}

county_year_tcc <- summarize_county_tcc(prediction_years)

if (file.exists(direct_estimates_file) && !force) {
  direct_estimates <- read_csv(direct_estimates_file, show_col_types = FALSE)
} else {
  direct_estimates <- unit_plots |>
    group_by(county_fips, county, fia_unit, fia_unit_name, year) |>
    summarise(
      n = n(),
      direct_biomass = mean(agb_live_mg_ha, na.rm = TRUE),
      direct_biomass_se = sd(agb_live_mg_ha, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) |>
    right_join(
      expand_grid(
        county_lookup |> select(county_fips, county, fia_unit, fia_unit_name),
        year = prediction_years
      ),
      by = c("county_fips", "county", "fia_unit", "fia_unit_name", "year")
    ) |>
    mutate(
      n = replace_na(n, 0L),
      direct_biomass = ifelse(n > 0, direct_biomass, NA_real_),
      direct_biomass_se = ifelse(n > 1, direct_biomass_se, NA_real_),
      direct_biomass_var = direct_biomass_se^2,
      n_eff = n
    ) |>
    left_join(
      county_year_tcc |> select(county_fips, year, county_mean_tcc, n_tcc_pixels),
      by = c("county_fips", "year")
    ) |>
    arrange(county_fips, year)

  positive_var <- direct_estimates$direct_biomass_var[
    is.finite(direct_estimates$direct_biomass_var) &
      direct_estimates$direct_biomass_var > 0
  ]
  vhat_floor <- if (length(positive_var)) {
    as.numeric(quantile(positive_var, probs = 0.05, na.rm = TRUE))
  } else {
    NA_real_
  }

  direct_estimates <- direct_estimates |>
    mutate(
      direct_biomass_vhat = case_when(
        is.na(direct_biomass) ~ NA_real_,
        is.finite(direct_biomass_var) & direct_biomass_var > 0 ~ direct_biomass_var,
        is.finite(vhat_floor) & vhat_floor > 0 ~ vhat_floor,
        TRUE ~ NA_real_
      ),
      direct_biomass_vhat_floor = vhat_floor,
      direct_biomass_vhat_source = case_when(
        is.na(direct_biomass) ~ NA_character_,
        is.finite(direct_biomass_var) & direct_biomass_var > 0 ~ "direct_biomass_var",
        is.finite(vhat_floor) & vhat_floor > 0 ~ "floor",
        TRUE ~ NA_character_
      ),
      direct_estimate_status = case_when(
        n == 0 ~ "no_plots",
        n == 1 ~ "single_plot_no_variance",
        !is.na(direct_biomass) & !is.na(direct_biomass_vhat) ~ "modeled",
        !is.na(direct_biomass) ~ "missing_model_variance",
        TRUE ~ "missing"
      ),
      direct_estimate_in_model = direct_estimate_status == "modeled",
      direct_biomass_model = if_else(direct_estimate_in_model, direct_biomass, NA_real_),
      direct_biomass_se_model = if_else(direct_estimate_in_model, direct_biomass_se, NA_real_),
      direct_biomass_vhat_model = if_else(direct_estimate_in_model, direct_biomass_vhat, NA_real_),
      n_eff_model = if_else(direct_estimate_in_model, as.numeric(n_eff), NA_real_)
    )

  write_csv_logged(direct_estimates, direct_estimates_file)
}

## ---- manifest --------------------------------------------------------------

manifest_file <- file.path(output_dir, "DATA_MANIFEST.csv")

manifest <- tibble(
  key = c(
    "created_at",
    "fia_data_dir",
    "stunitco_gpkg",
    "tcc_input_dir",
    "tcc_output_dir",
    "reuse_pnw_tcc_blocks",
    "clip_reused_pnw_tcc_blocks",
    "pnw_block_dir",
    "prediction_years",
    "tcc_window_radius_m",
    "tcc_kernel_size_pixels",
    "county_mean_tcc_source",
    "grid_cells_per_county_year",
    "counties",
    "fia_units",
    "unit_plot_rows",
    "unit_plot_missing_tcc_rows",
    "unit_plot_dropped_missing_tcc_rows",
    "direct_estimate_rows",
    "direct_estimate_model_rows",
    "direct_estimate_single_plot_excluded_rows",
    "direct_estimate_variance_floor_rows",
    "direct_estimate_missing_county_tcc_rows",
    "direct_estimate_missing_vhat_rows",
    "direct_estimate_missing_model_vhat_rows",
    paste0("prediction_grid_", names(grid_full_rows), "_full_rows"),
    paste0("prediction_grid_", names(grid_sample_rows), "_sample_rows")
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    input_fia_dir,
    stunitco_file,
    tcc_input_dir,
    tcc_output_dir,
    as.character(reuse_pnw_blocks),
    as.character(clip_reused_pnw_blocks),
    pnw_block_dir,
    paste(range(prediction_years), collapse = "-"),
    as.character(window_radius_m),
    as.character(kernel_size),
    "exactextractr_zonal_mean_full_annual_block_raster",
    as.character(grid_cells_per_county_year),
    as.character(nrow(wa_counties)),
    paste(sort(unique(county_lookup$fia_unit_name)), collapse = ";"),
    as.character(nrow(unit_plots)),
    as.character(sum(is.na(unit_plots$tcc_mean))),
    as.character(unit_plot_dropped_missing_tcc_rows),
    as.character(nrow(direct_estimates)),
    as.character(sum(direct_estimates$direct_estimate_in_model)),
    as.character(sum(direct_estimates$direct_estimate_status == "single_plot_no_variance")),
    as.character(sum(
      direct_estimates$direct_estimate_in_model &
        direct_estimates$direct_biomass_vhat_source == "floor",
      na.rm = TRUE
    )),
    as.character(sum(is.na(direct_estimates$county_mean_tcc))),
    as.character(sum(!is.na(direct_estimates$direct_biomass) & is.na(direct_estimates$direct_biomass_vhat))),
    as.character(sum(!is.na(direct_estimates$direct_biomass_model) & is.na(direct_estimates$direct_biomass_vhat_model))),
    as.character(grid_full_rows),
    as.character(grid_sample_rows)
  )
)

write_csv_logged(manifest, manifest_file)

log_msg("Finished WA vignette data build.")
