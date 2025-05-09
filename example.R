file <- "/vsicurl/https://github.com/mdsumner/gproj/releases/download/v0.0.1/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc"

file <- normalizePath("~/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc")
dsn <- sprintf("NetCDF:\"%s\":/geophysical_data/bbp_443", file)
library(gdalraster)
s <- 1.2e7
lon <- seq(-160, -40, by = 5)
i <- 12
source("R/gproj.R")

w <- warp(dsn, tf <- tempfile(tmpdir = "/vsimem", fileext = ".tif"), prj, cl_arg = c("-te", -s, -s, s, s, "-ts",  1024, 0))
elev_pal <- hcl.colors(256)
ramp <- scales::colour_ramp(elev_pal, alpha=FALSE)

plot_raster(new(GDALRaster, tf), col_map_fn = ramp)


mp <- do.call(cbind, maps::map(plot = F)[1:2])


plot(proj_xy(viirs_boundary), pch = ".", col = "firebrick", asp = 1)
lines(proj_xy(viirs_boundary), col = "firebrick")
points(proj_xy(midll(viirs_boundary)))
points(proj_xy(mp), pch = ".", col = "grey")


ex <- c(apply(proj_xy(m), 2, range))
library(gdalraster)
w <- warp(dsn, tf <- tempfile(tmpdir = "/vsimem", fileext = ".tif"), prj, cl_arg = c("-te", ex[1], ex[3], ex[2], ex[4], "-ts",  1024, 0))
plot_raster(new(GDALRaster, tf), col_map_fn = ramp)

