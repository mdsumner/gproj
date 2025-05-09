#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


# VIIRS swath sources
#
# Data frame of 'source' urls to VIIRs swaths in HDF5.
#
# To read these directly we can't. The oc nasa getfile mech is not streamable. See viirs_eg for an example.
#
# Created in data-raw/
# @name viirs_swaths
# @docType data
#'viirs_swaths'

#' An example swath source
#'
#' Only for use in examples in this package.
#'
#' @returns string to a url with vsicurl protocol
#' @export
#'
#' @examples
#' viirs_eg() ## open this with 'new(gdalraster::GDALRaster, <>)' or terra::rast
viirs_eg <- function() {
  "/vsicurl/https://github.com/mdsumner/gproj/releases/download/v0.0.1/S3B_OLCI_ERRNT.20250101T164657.L2.IOP.NRT.nc"
}
