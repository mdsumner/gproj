  nc <- tidync::tidync("dt_ref_global_merged_madt_uv_19921014_19921014_20060315.nc")
lat_from_merc <- nc$transforms$NbLatitudes$NbLatitudes
lon_seq <- nc$transforms$NbLongitudes$NbLongitudes


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

