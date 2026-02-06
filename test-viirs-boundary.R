## gproj Phase 1: params() and proj() with VIIRS swath boundary
##
## Run from the gproj package root, or adjust the source/load paths

#source("R/gproj.R")
load("data/viirs_boundary.rda")

cat("viirs_boundary:", nrow(viirs_boundary), "points\n")
cat("  lon range:", round(range(viirs_boundary[,1]), 2), "\n")
cat("  lat range:", round(range(viirs_boundary[,2]), 2), "\n\n")

# =============================================================================
# 1. params() — what does it see?
# =============================================================================

p <- params(viirs_boundary)
print(p)

cat("\nAxis details:\n")
cat("  azimuth:", round(p@axis$azimuth, 2), "°\n")
cat("  tilt:", round(p@axis$tilt, 2), "°\n")
cat("  variance ratio:", round(p@axis$variance_ratio, 4), "\n")
cat("  PC1 lonlat:", round(p@axis$pc1_lonlat, 2), "\n")
cat("  centroid (spherical):", round(p@axis$lonc, 2), round(p@axis$latc, 2), "\n")

# =============================================================================
# 2. proj() — auto and explicit families
# =============================================================================

cat("\n=== Projection strings ===\n")
cat("auto:       ", proj(p, "auto"), "\n")
cat("equal_area: ", proj(p, "equal_area"), "\n")
cat("conformal:  ", proj(p, "conformal"), "\n")
cat("equidistant:", proj(p, "equidistant"), "\n")
cat("stere:      ", proj(p, "stere"), "\n")
cat("laea:       ", proj(p, "laea"), "\n")
cat("ortho:      ", proj(p, "ortho"), "\n")
cat("omerc:      ", proj(p, "omerc"), "\n")
cat("lcc:        ", proj(p, "lcc"), "\n")

# =============================================================================
# 3. Compare: original README omerc vs PCA-derived omerc
# =============================================================================

cat("\n=== Oblique Mercator comparison ===\n")

## The hand-tuned omerc from the README
readme_omerc <- "+proj=omerc +lonc=-113 +lat_0=-20 +alpha=15 +gamma=1"
cat("README (hand-tuned):", readme_omerc, "\n")

## PCA-derived omerc
pca_omerc <- proj(p, "omerc")
cat("PCA-derived:        ", pca_omerc, "\n")

cat("\nDifference tells us how well PCA recovers the swath alignment.\n")
cat("The README alpha=15 was eyeballed. PCA should find something close.\n")

# =============================================================================
# 4. Visual comparison — project and plot
# =============================================================================

mp <- do.call(cbind, maps::map(plot = FALSE)[1:2])

families <- c("stere", "laea", "ortho", "omerc")
par(mfrow = c(2, 2), mar = c(1, 1, 2, 1))

for (fam in families) {
  prj <- proj(p, fam)
  g <- gproj(prj)

  xy_viirs <- g@proj_xy(viirs_boundary)
  xy_coast <- g@proj_xy(mp)

  ## Clip coastline to reasonable range (projections can blow up far from centre)
  ## Use 3x the viirs extent as clip window
  vx <- range(xy_viirs[,1], na.rm = TRUE)
  vy <- range(xy_viirs[,2], na.rm = TRUE)
  pad <- 0.5
  xlim <- vx + diff(vx) * c(-pad, pad)
  ylim <- vy + diff(vy) * c(-pad, pad)

  ok <- xy_coast[,1] > xlim[1] & xy_coast[,1] < xlim[2] &
    xy_coast[,2] > ylim[1] & xy_coast[,2] < ylim[2]
  ok[is.na(ok)] <- FALSE

  plot(xy_viirs, cex = 0.2, col = "firebrick", asp = 1,
       xlim = xlim, ylim = ylim, xlab = "", ylab = "",
       main = sprintf("%s (auto from params)", fam))
  points(xy_coast[ok,], pch = ".", col = "grey40")
}

# =============================================================================
# 5. The omerc alignment test
# =============================================================================

cat("\n=== Omerc alignment quality ===\n")

g_omerc <- gproj(pca_omerc)
xy_omerc <- g_omerc@proj_xy(viirs_boundary)

x_range <- diff(range(xy_omerc[,1], na.rm = TRUE))
y_range <- diff(range(xy_omerc[,2], na.rm = TRUE))
cat("In omerc coordinates:\n")
cat("  X span:", round(x_range / 1000), "km\n")
cat("  Y span:", round(y_range / 1000), "km\n")
cat("  Aspect (X/Y):", round(x_range / y_range, 2), "\n")
cat("  (>> 1 means the swath is well-aligned along omerc X axis)\n")

## Compare with README omerc
g_readme <- gproj(readme_omerc)
xy_readme <- g_readme@proj_xy(viirs_boundary)
x_range_r <- diff(range(xy_readme[,1], na.rm = TRUE))
y_range_r <- diff(range(xy_readme[,2], na.rm = TRUE))
cat("\nREADME omerc aspect (X/Y):", round(x_range_r / y_range_r, 2), "\n")
cat("PCA omerc aspect (X/Y):   ", round(x_range / y_range, 2), "\n")
cat("(closer to the true swath aspect = better alignment)\n")

# =============================================================================
# 6. Does auto pick omerc for this swath?
# =============================================================================

cat("\n=== Auto family selection ===\n")
cat("auto_family picks:", proj(p, "auto"), "\n")
cat("Variance ratio[1]:", round(p@axis$variance_ratio[1], 3), "\n")
cat("Axis tilt:", round(p@axis$tilt, 1), "°\n")
cat("\nFor omerc to be auto-selected, need variance_ratio > 0.85 AND tilt > 20.\n")
cat("If auto didn't pick omerc, we may need to tune those thresholds.\n")

# =============================================================================
# 7. Quick gproj_region for context
# =============================================================================

cat("\n=== gproj_region at swath centre ===\n")
reg <- gproj_region(
  lon = p@centre["lon_0"],
  lat = p@centre["lat_0"],
  radius_km = 3000,
  property = "equal_area"
)
cat("Region CRS:", reg$crs, "\n")
gs <- grid_spec(reg, resolution_m = 10000)
cat("Grid:", gs$dimension["ncol"], "x", gs$dimension["nrow"],
    "at", gs$resolution["dx"]/1000, "km\n")

# =============================================================================
# 8. Secant params for LCC
# =============================================================================

cat("\n=== LCC with secant latitudes ===\n")
p_sec <- params(viirs_boundary, secant = TRUE)
lcc_crs <- proj(p_sec, "lcc")
cat("LCC:", lcc_crs, "\n")
cat("Secant lats:", round(p_sec@secant, 2), "\n")

