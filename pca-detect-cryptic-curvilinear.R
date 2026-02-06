## detect_grid demo: CryoSat NSIDC "cryptic curvilinear" → EPSG:3412
##
## Requires: the CS2WFA_25km_201007.nc file from the hypertidy blog post
## source("R/gproj.R")

library(terra)

## netcdf file from a blog post critiquing its coordinate representation
## https://zenodo.org/records/7327711/files/CS2WFA_25km_201007.nc?download=1
## https://www.hypertidy.org/posts/2025-09-04_broken-netcdf/
#system("wget https://zenodo.org/records/7327711/files/CS2WFA_25km_201007.nc?download=1 -O CS2WFA_25km_201007.nc")
#piggyback::pb_upload("CS2WFA_25km_201007.nc")


dsn <- "/vsicurl/https://github.com/mdsumner/gproj/releases/download/v0.0.1/CS2WFA_25km_201007.nc"

## 1. Read the lon/lat coordinate arrays
lon <- values(r0 <- rast(sprintf("NETCDF:%s:lon", dsn)))
lat <- values(rast(sprintf("NETCDF:%s:lat", dsn)))

## reshape to 2D (316 cols x 332 rows from the file metadata)
nc <- dim(r0)[2L]; nr <- dim(r0)[1L]

lon_2d <- matrix(lon, nrow = nr, ncol = nc, byrow = TRUE)
lat_2d <- matrix(lat, nrow = nr, ncol = nc, byrow = TRUE)

cat("Grid dimensions:", nr, "x", nc, "\n")
cat("Total coordinate values:", nr * nc * 2, "\n\n")

## 2. Is it rectilinear in longlat?
rect <- detect_rectilinear(lon_2d, lat_2d)
cat("Rectilinear in lonlat?", rect$is_rectilinear, "\n\n")

## 3. Run params on the coordinates
p <- params(lon = lon_2d, lat = lat_2d)
print(p)

## 4. Project to candidate CRS and test regularity
## params sees polar_south → try polar stereographic
candidate <- "EPSG:3412"
xy_proj <- reproj::reproj_xy(cbind(as.vector(lon_2d), as.vector(lat_2d)),
                             target = candidate, source = "EPSG:4326")

x_2d <- matrix(xy_proj[, 1], nrow = nr, ncol = nc)
y_2d <- matrix(xy_proj[, 2], nrow = nr, ncol = nc)

## 5. Test: are rows horizontal in projected space? (y constant along each row)
row_y_sd <- apply(y_2d, 1, sd, na.rm = TRUE)
cat("Max row-y SD:", round(max(row_y_sd, na.rm = TRUE), 2), "m\n")
cat("  (should be near 0 if rows are horizontal)\n\n")

## 6. Test: is x-spacing constant along each row?
row_dx <- t(apply(x_2d, 1, diff))
dx_median <- median(row_dx, na.rm = TRUE)
dx_dev <- max(abs(row_dx - dx_median), na.rm = TRUE)
cat("Median dx:", round(dx_median), "m\n")
cat("Max dx deviation:", round(dx_dev, 2), "m\n\n")

## 7. Test: is y-spacing constant along each column?
col_dy <- apply(y_2d, 2, diff)
dy_median <- median(col_dy, na.rm = TRUE)
dy_dev <- max(abs(col_dy - dy_median), na.rm = TRUE)
cat("Median dy:", round(dy_median), "m\n")
cat("Max dy deviation:", round(dy_dev, 2), "m\n\n")

## 8. Snap resolution
snap <- function(x, nice = c(1000, 5000, 10000, 12500, 25000, 50000, 100000)) {
  nice[which.min(abs(nice - abs(x)))]
}
res_x <- snap(dx_median)
res_y <- snap(abs(dy_median))
cat("Snapped resolution:", res_x, "x", res_y, "m\n\n")

## 9. Recover extent (cell centres → grid edges)
x_centres <- sort(unique(round(x_2d[1, ])))  # one row suffices
y_centres <- sort(unique(round(y_2d[, 1])))   # one column suffices

xmin <- min(x_centres) - res_x / 2
xmax <- max(x_centres) + res_x / 2
ymin <- min(y_centres) - res_y / 2
ymax <- max(y_centres) + res_y / 2

cat("Recovered extent:\n")
cat("  xmin:", xmin, "xmax:", xmax, "\n")
cat("  ymin:", ymin, "ymax:", ymax, "\n\n")

## 10. Compare to known NSIDC grid
cat("Known NSIDC extent: -3950000 3950000 -3950000 4350000\n")
cat("Difference: ", xmin - (-3950000), xmax - 3950000,
    ymin - (-3950000), ymax - 4350000, "\n\n")

## 11. The punchline
cat("=== Summary ===\n")
cat("CRS:         ", candidate, "\n")
cat("Resolution:  ", res_x, "x", res_y, "m\n")
cat("Dimensions:  ", nc, "x", nr, "\n")
cat("Extent:      ", xmin, xmax, ymin, ymax, "\n")
cat("Replaced:    ", nr * nc * 2, "coordinate values with 4 numbers + 1 CRS\n")
cat("Entropy ratio:", round(nr * nc * 2 / 4), ":1\n")


## ============================================================================
## EXPECTED OUTPUT (approximate — values depend on float precision in the file)
## ============================================================================
#
# Grid dimensions: 332 x 316
# Total coordinate values: 209824
#
# Rectilinear in lonlat? FALSE
#
# <params: 104912 points>
#   centre:     0.XX, -YY.YY           ← should be near south pole
#   extent:     [-180, 180] x [-90, -39]
#   span:       ~360° lon x ~50° lat
#   hemisphere: polar_south
#   aspect:     wide
#   axis:       azimuth=...°, tilt=...°, var_ratio=0.50/0.50/0.00
#                ↑ polar data has no dominant axis — variance split evenly in x/y
#
# Max row-y SD: ~2-10 m
#   (should be near 0 if rows are horizontal)
#
# Median dx: ~25000 m
# Max dx deviation: ~1-10 m          ← floating point noise from trig roundtrip
#
# Median dy: ~-25000 m               ← negative because y decreases with row index
# Max dy deviation: ~1-10 m
#
# Snapped resolution: 25000 x 25000 m
#
# Recovered extent:
#   xmin: ~-3950000 xmax: ~3950000
#   ymin: ~-3950000 ymax: ~4350000
#
# Known NSIDC extent: -3950000 3950000 -3950000 4350000
# Difference:  ~0 ~0 ~0 ~0           ← if snapping works, exact match
#
# === Summary ===
# CRS:          EPSG:3412
# Resolution:   25000 x 25000 m
# Dimensions:   316 x 332
# Extent:       -3950000 3950000 -3950000 4350000
# Replaced:     209824 coordinate values with 4 numbers + 1 CRS
# Entropy ratio: 52456 :1
