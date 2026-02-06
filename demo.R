## gproj Phase 1 demo
## Run this to exercise all the new functionality

# source("R/gproj.R")  ## or devtools::load_all()

# =============================================================================
# 1. Basic gproj: forward/inverse
# =============================================================================

cat("=== 1. gproj forward/inverse ===\n")

## Identity
g0 <- gproj()
g0@proj_xy(cbind(147, -42))

## LAEA centred on Hobart
g1 <- gproj("+proj=laea +lon_0=147 +lat_0=-42")
xy <- g1@proj_xy(cbind(147, -42))
print(xy)  # should be ~0,0
g1@proj_ll(xy)  # should be 147, -42

## Stereographic for Antarctica
g2 <- gproj("+proj=stere +lat_0=-90 +lon_0=0 +lat_ts=-71 +datum=WGS84")


# =============================================================================
# 2. params() on various geometries
# =============================================================================

cat("\n=== 2. params() examples ===\n")

## A. Small region around Hobart
hobart_pts <- cbind(
  lon = 147 + runif(100, -1, 1),
  lat = -42 + runif(100, -0.5, 0.5)
)
(p_hobart <- params(hobart_pts))
cat("Auto projection:", proj(p_hobart), "\n")
cat("Equal area:", proj(p_hobart, "equal_area"), "\n")

## B. Antarctic polar cap
polar_pts <- cbind(
  lon = runif(200, -180, 180),
  lat = runif(200, -90, -60)
)
(p_polar <- params(polar_pts))
cat("Auto projection:", proj(p_polar), "\n")

## C. Mid-latitude wide band (like Australia)
aus_pts <- cbind(
  lon = seq(112, 154, length.out = 100),
  lat = seq(-10, -44, length.out = 100)
)
(p_aus <- params(aus_pts))
cat("Auto projection:", proj(p_aus), "\n")

## D. Elongated track (simulated satellite swath, oblique)
## A great circle from (-113, -20) heading NE
t_seq <- seq(0, 0.3, length.out = 200)
track_lon <- -113 + 30 * t_seq + rnorm(200, 0, 0.5)
track_lat <- -20 + 50 * t_seq + rnorm(200, 0, 0.5)
track_pts <- cbind(track_lon, track_lat)
(p_track <- params(track_pts))
cat("Auto projection:", proj(p_track), "\n")
cat("Oblique Mercator:", proj(p_track, "omerc"), "\n")
cat("  -> axis azimuth:", round(p_track@axis$azimuth, 1), "degrees\n")
cat("  -> variance ratio:", round(p_track@axis$variance_ratio[1], 3), "\n")


# =============================================================================
# 3. gproj_region: quick "give me a projection here"
# =============================================================================

cat("\n=== 3. gproj_region() ===\n")

## Mawson Station
mawson <- gproj_region(lon = 62.87, lat = -67.6, radius_km = 500, property = "equal_area")
cat("Mawson LAEA:", mawson$crs, "\n")
cat("Extent:", mawson$extent, "\n")

## Grid spec for 1km resolution
gs <- grid_spec(mawson, resolution_m = 1000)
cat("Grid:", gs$dimension["ncol"], "x", gs$dimension["nrow"], "at", gs$resolution["dx"], "m\n")

## Davis Station, conformal
davis <- gproj_region(lon = 77.97, lat = -68.58, radius_km = 200, property = "conformal")
cat("Davis stereo:", davis$crs, "\n")


# =============================================================================
# 4. Mercator devolution detection
# =============================================================================

cat("\n=== 4. Mercator devolution detection ===\n")

## Simulate the AVISO MADT grid:
## Uniform longitude, Mercator-spaced latitude
## The real file has 1/3 degree Mercator resolution

R <- 6378137  # WGS84 semi-major axis

## Build Mercator y values at uniform spacing
merc_res <- R * (1/3) * pi / 180  # ~37km at equator
merc_y_seq <- seq(-merc_res * 100, merc_res * 100, by = merc_res)

## Invert to get latitudes: phi = 2*atan(exp(y/R)) - pi/2
lat_from_merc <- (2 * atan(exp(merc_y_seq / R)) - pi/2) * 180 / pi

## Uniform longitudes
lon_seq <- seq(0, 359.6667, by = 1/3)

cat("Simulated MADT-like grid:\n")
cat("  Lon:", length(lon_seq), "values, range:", range(lon_seq), "\n")
cat("  Lat:", length(lat_from_merc), "values, range:", round(range(lat_from_merc), 2), "\n")

## Test the latitude sequence
merc_test <- test_mercator_spacing(lat_from_merc)
cat("  Mercator test:", merc_test$is_mercator, "\n")
cat("  Max relative deviation:", merc_test$max_relative_deviation, "\n")
cat("  Resolution:", round(merc_test$merc_resolution_m), "m\n")

## Test longitude regularity
lon_test <- test_regular_spacing(lon_seq)
cat("  Lon regular:", lon_test$is_regular, "\n")
cat("  Lon resolution:", round(lon_test$resolution, 4), "degrees\n")

## Now build 2D arrays and run the full rectilinear detection
lon_2d <- matrix(rep(lon_seq, each = length(lat_from_merc)),
                 nrow = length(lat_from_merc), ncol = length(lon_seq))
lat_2d <- matrix(rep(lat_from_merc, times = length(lon_seq)),
                 nrow = length(lat_from_merc), ncol = length(lon_seq))

result <- detect_rectilinear(lon_2d, lat_2d)
cat("\ndetect_rectilinear() result:\n")
cat("  is_rectilinear:", result$is_rectilinear, "\n")
cat("  type:", result$type, "\n")
cat("  orientation:", result$orientation, "\n")
cat("  Mercator CRS:", result$mercator_crs, "\n")
cat("  Mercator resolution:", round(result$mercator_resolution), "m\n")

## Contrast with a genuinely regular longlat grid
cat("\n--- Regular longlat grid ---\n")
lon_reg <- seq(-180, 179, by = 1)
lat_reg <- seq(-90, 89, by = 1)
lon_2d_r <- matrix(rep(lon_reg, each = length(lat_reg)),
                   nrow = length(lat_reg), ncol = length(lon_reg))
lat_2d_r <- matrix(rep(lat_reg, times = length(lon_reg)),
                   nrow = length(lat_reg), ncol = length(lon_reg))
result_r <- detect_rectilinear(lon_2d_r, lat_2d_r)
cat("  type:", result_r$type, "\n")  # should be "regular_lonlat"


# =============================================================================
# 5. PCA axis detection with the VIIRS-like example
# =============================================================================

cat("\n=== 5. PCA axis on oblique swath ===\n")

## Simulate a VIIRS-like swath: narrow, oblique, roughly NW-SE
## (like the example in the README)
n_along <- 300
n_across <- 40
along_t <- seq(0, 1, length.out = n_along)
across_t <- seq(-0.02, 0.02, length.out = n_across)

swath_lon <- outer(across_t * 50, along_t * 60 - 143, "+")
swath_lat <- outer(across_t * 30, along_t * 50 - 45, "+")

swath_xy <- cbind(as.vector(swath_lon), as.vector(swath_lat))

p_swath <- params(swath_xy)
cat("Swath centre:", round(p_swath@centre, 2), "\n")
cat("Axis azimuth:", round(p_swath@axis$azimuth, 1), "\n")
cat("Axis tilt:", round(p_swath@axis$tilt, 1), "\n")
cat("Variance ratio:", round(p_swath@axis$variance_ratio, 3), "\n")

## The omerc projection should align with the swath
omerc_crs <- proj(p_swath, "omerc")
cat("Oblique Mercator:", omerc_crs, "\n")

## Project and check: the swath should now be roughly axis-aligned
g_omerc <- gproj(omerc_crs)
swath_proj <- g_omerc@proj_xy(swath_xy)
cat("Projected extent X:", round(range(swath_proj[,1])), "\n")
cat("Projected extent Y:", round(range(swath_proj[,2])), "\n")
cat("X span / Y span:", round(diff(range(swath_proj[,1])) / diff(range(swath_proj[,2])), 2), "\n")
cat("  (should be >> 1 if omerc aligned the swath along X)\n")

cat("\n=== Done ===\n")

