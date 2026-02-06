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
swath_xy <- viirs_boundary
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
plot(swath_proj, pch = 19, cex = .3, col = "firebrick", asp = 1)
lines(reproj::reproj_xy(do.call(cbind, maps::map(plot = F)[1:2]), omerc_crs, source = "EPSG:4326"))
