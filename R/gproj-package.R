#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


# VIIRS swath sources
#
# Data frame of 'source' urls to VIIRs swaths in HDF4, along with 'size' the file size in bytes.
#
# To read these directly we must set Earthdata authorization (i.e. GDAL_HTTP_HEADER_FILE, GDAL_HTTP_HEADERS, or perhaps a .netrc file).
#
# Created in data-raw/
# @name viirs_swaths
# @docType data
#'viirs_swaths'

# An example swath source
#
# Only for use in examples in this package.
#
# @returns string to a url with vsicurl protocol
# @export
#
# @examples
# viirs_eg() ## open this with 'new(gdalraster::GDALRaster, <>)' or terra::rast
