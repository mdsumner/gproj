
bound <- vaster::boundary_cell

f <- "~/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc"
library(terra)
ll <- flip(c(rast(f, "/navigation_data/longitude"),
  rast(f, "/navigation_data/latitude")), "horizontal")

plot(ll[bound(dim(ll)[2:1])], pch = ".")
maps::map(add = T)
r0 <- rast(f, "/geophysical_data/bbp_443")
set.crs(r0, "EPSG:4326")
plot(project(r0, rast(ext(160, 180, -90, -60), res = .2, crs = "EPSG:4326"), by_util = TRUE, method = "sum"))
plot(project(r0, rast(ext(-180, -140, -80, 0), res = .2, crs = "EPSG:4326"), by_util = TRUE, method = "sum"))

dsn <- sprintf("NetCDF:\"%s\":/geophysical_data/bbp_443", normalizePath(f))
library(gdalraster)
s <- 1.2e7
lon <- seq(-160, -40, by = 5)
#i <- 0
i <- 12
prj <- sprintf("+proj=laea +lon_0=%i +lat_0=-10", lon[i])
w <- warp(dsn, tf <- tempfile(tmpdir = "/vsimem", fileext = ".tif"), prj, cl_arg = c("-te", -s, -s, s, s, "-ts",  1024, 0))
plot_raster(read_ds(new(GDALRaster, tf)))
# plot(reproj::reproj_xy(ll[bound(dim(ll)[2:1])], prj), pch = ".", col = "firebrick")
points(reproj::reproj_xy(ll[bound(dim(ll)[2:1])], prj), pch = ".", col = "firebrick")

m <- reproj::reproj_xy(do.call(cbind, maps::map(plot = F)[1:2]), prj)
lines(m)




plot(proj_xy(m), pch = ".", col = "firebrick", asp = 1)
lines(proj_xy(m), col = "firebrick")
points(proj_xy(midll(m)))
mp <- do.call(cbind, maps::map(plot = F)[1:2])
points(proj_xy(mp), pch = ".", col = "grey")


ex <- c(apply(proj_xy(m), 2, range))
dsn <- sprintf("NetCDF:\"%s\":/geophysical_data/bbp_443", normalizePath(f))
library(gdalraster)
w <- warp(dsn, tf <- tempfile(tmpdir = "/vsimem", fileext = ".tif"), prj, cl_arg = c("-te", ex[1], ex[3], ex[2], ex[4], "-ts",  1024, 0))
plot_raster(read_ds(new(GDALRaster, tf)))
image(sqrt(rast(tf)), col = hcl.colors(256))
