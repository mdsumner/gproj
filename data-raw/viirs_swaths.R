Sys.setenv("GDAL_HTTP_HEADER_FILE" = normalizePath("~/earthdata"),
          "GDAL_DISABLE_READDIR_ON_OPEN" = "EMPTY_DIR")

base <- "http://oceandata.sci.gsfc.nasa.gov/getfile"
#http://oceandata.sci.gsfc.nasa.gov/getfile/S3B_OLCI_EFRNT.20250101T000134.L2.OC.nc

allh5 <- readLines("data-raw/h5.tx")
dsn <- sprintf("/vsicurl/%s", file.path(base, allh5))


viirs_swaths <- tibble::tibble(source = dsn))



## biggest file is
viirs_swaths |> dplyr::filter(stringr::str_detect(source, "ERRNT")) |> dplyr::slice(which.max(size))
url <- "/vsicurl/https://n5eil01u.ecs.nsidc.org/VIIRS/VNP10A1.002/2025.05.06/VNP10A1.A2025126.h23v02.002.2025127112146.h5"
header <- readLines("~/earthdata")
tf <- tempfile(fileext = ".h5")
cmd <- sprintf("wget --header=\"%s\" %s -O %s", header, gsub("/vsicurl/", "", url), tf)
system(cmd)
Sys.setenv("GDAL_HTTP_HEADER_FILE" = normalizePath("~/earthdata"),
          "GDAL_DISABLE_READDIR_ON_OPEN" = "EMPTY_DIR")

library(terra)
ll <- terra::flip(c(rast(tf, "/navigation_data/longitude"),
  rast(file, "/navigation_data/latitude")), "horizontal")

viirs_boundary <- ll[vaster:::vaster_boundary_cell(dim(ll)[2:1])]


