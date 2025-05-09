Sys.setenv("GDAL_HTTP_HEADER_FILE" = normalizePath("~/earthdata"),
          "GDAL_DISABLE_READDIR_ON_OPEN" = "EMPTY_DIR")

base <- "http://oceandata.sci.gsfc.nasa.gov/getfile"
#http://oceandata.sci.gsfc.nasa.gov/getfile/S3B_OLCI_EFRNT.20250101T000134.L2.OC.nc

allh5 <- readLines("data-raw/h5.tx")
dsn <- sprintf("/vsicurl/%s", file.path(base, allh5))


viirs_swaths <- tibble::tibble(source = dsn))



## biggest file is
# viirs_swaths |> dplyr::filter(stringr::str_detect(source, "ERRNT")) |> dplyr::slice(which.max(size))
# url <- "/vsicurl/https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1.002/2025.05.06/VNP10A1.A2025126.h23v02.002.2025127112146.h5"
url <- sprintf("/vsicurl/%s", file.path(base, "S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc"))

header <- readLines("~/earthdata")
tf <- tempfile(fileext = ".h5")
cmd <- sprintf("wget --header=\"%s\" %s -O %s", header, gsub("/vsicurl/", "", url), tf)
system(cmd)

## didn't exactly go like that above so use
## https://github.com/mdsumner/gproj/releases/download/v0.0.1/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc
library(piggyback)
pb_new_release("mdsumner/gproj", "v0.0.1")
pb_upload("~/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc",  repo = "mdsumner/gproj", tag = "v0.0.1")


file <- "/vsicurl/https://github.com/mdsumner/gproj/releases/download/v0.0.1/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc"
#library(terra)
#ll <- terra::flip(c(rast(file, "/navigation_data/longitude"),
#  rast(file, "/navigation_data/latitude")), "horizontal")

lon <- new(gdalraster::GDALRaster, sprintf("vrt://%s?sd_name=/navigation_data/longitude", file))
lat <- new(gdalraster::GDALRaster, file)
cell <- vaster:::vaster_boundary_cell(lon$dim()[2:1])
library(gdalraster)
## wasteful but whatever, terra definitly good at this if the file is not remote
viirs_boundary <- cbind(read_ds(lon)[cell], read_ds(lat)[cell])


plot(viirs_boundary)
usethis::use_data(viirs_boundary)
