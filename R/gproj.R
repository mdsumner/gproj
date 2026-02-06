## gproj: projection concierge
## S7 classes and methods for easy projection selection and coordinate transformation
##
## Phase 1: core gproj class, enriched params(), proj(), gproj_region()

library(S7)

# =============================================================================
# gproj class: stores a CRS and provides forward/inverse projection
# =============================================================================

gproj <- new_class("gproj",
                   properties = list(
                     crs = class_character,
                     proj_xy = class_function,
                     proj_ll = class_function
                   ),
                   constructor = function(crs = "") {
                     if (nchar(crs) == 0) {
                       ## identity: no projection
                       new_object(S7_object(),
                                  crs = crs,
                                  proj_xy = function(xy) { colnames(xy) <- c("x", "y"); xy },
                                  proj_ll = function(xy) { colnames(xy) <- c("x", "y"); xy }
                       )
                     } else {
                       ll <- getOption("reproj.default.longlat", "OGC:CRS84")
                       new_object(S7_object(),
                                  crs = crs,
                                  proj_xy = function(xy) {
                                    reproj::reproj_xy(xy, target = crs, source = ll)
                                  },
                                  proj_ll = function(xy) {
                                    reproj::reproj_xy(xy, target = ll, source = crs)
                                  }
                       )
                     }
                   }
)


# =============================================================================
# params class: characterises the geometry of a set of coordinates
# =============================================================================

params <- new_class("params",
                    properties = list(
                      centre = class_numeric,       # c(lon_0, lat_0)
                      extent_ll = class_numeric,    # c(xmin, xmax, ymin, ymax) in longlat
                      span = class_numeric,         # c(lon_span, lat_span) in degrees
                      axis = class_list,            # PCA-derived great circle axis (azimuth, components)
                      hemisphere = class_character,  # "polar_south", "polar_north", "mid_south", "mid_north", "equatorial"
                      aspect = class_character,     # "wide", "tall", "square" — shape of the region
                      secant = class_numeric,       # lat_1, lat_2 for conic projections (or empty)
                      n = class_integer             # number of input points
                    ),
                    constructor = function(xy = NULL, lon = NULL, lat = NULL, secant = FALSE) {
                      ## accept either a 2-col matrix or separate lon/lat
                      if (is.null(xy) && !is.null(lon) && !is.null(lat)) {
                        ## lon/lat might be 2D arrays (curvilinear grids) — flatten
                        xy <- cbind(as.vector(lon), as.vector(lat))
                      }
                      if (is.null(xy)) {
                        stop("provide either 'xy' (2-col matrix of lon,lat) or 'lon' and 'lat'")
                      }
                      xy <- as.matrix(xy[, 1:2])
                      ## remove NAs
                      ok <- complete.cases(xy)
                      xy <- xy[ok, , drop = FALSE]

                      n <- nrow(xy)
                      lon_vals <- xy[, 1]
                      lat_vals <- xy[, 2]

                      ## --- Centre ---
                      ## Use column-wise mean (fine for most cases, breaks near ±180)
                      ## TODO: circular mean for longitude when data spans antimeridian
                      centre <- c(lon_0 = mean(lon_vals), lat_0 = mean(lat_vals))

                      ## --- Extent ---
                      extent_ll <- c(
                        xmin = min(lon_vals), xmax = max(lon_vals),
                        ymin = min(lat_vals), ymax = max(lat_vals)
                      )

                      ## --- Span ---
                      lon_span <- diff(range(lon_vals))
                      lat_span <- diff(range(lat_vals))
                      span <- c(lon_span = lon_span, lat_span = lat_span)

                      ## --- Aspect ---
                      ## Compare angular spans (rough, but useful for projection family choice)
                      ratio <- lon_span / max(lat_span, 1e-6)
                      aspect <- if (ratio > 1.5) "wide" else if (ratio < 0.67) "tall" else "square"

                      ## --- Hemisphere / latitude band ---
                      mean_lat <- centre[2]
                      hemisphere <- if (mean_lat < -60) {
                        "polar_south"
                      } else if (mean_lat > 60) {
                        "polar_north"
                      } else if (mean_lat < -23.5) {
                        "mid_south"
                      } else if (mean_lat > 23.5) {
                        "mid_north"
                      } else {
                        "equatorial"
                      }

                      ## --- PCA axis (great circle alignment) ---
                      axis <- pca_axis(xy)

                      ## --- Secant latitudes ---
                      sec <- numeric(0)
                      if (isTRUE(secant)) {
                        sec <- c(lat_1 = extent_ll["ymin"], lat_2 = extent_ll["ymax"])
                        names(sec) <- c("lat_1", "lat_2")
                      }

                      new_object(S7_object(),
                                 centre = centre,
                                 extent_ll = extent_ll,
                                 span = span,
                                 axis = axis,
                                 hemisphere = hemisphere,
                                 aspect = aspect,
                                 secant = sec,
                                 n = as.integer(n)
                      )
                    }
)


# =============================================================================
# PCA axis detection: find the dominant great circle through a set of lon/lat
# =============================================================================

#' Compute the dominant axis of a set of lon/lat points via PCA on 3D cartesian
#'
#' Returns a list with:
#'   - azimuth: bearing of the dominant axis from the centroid (degrees from north)
#'   - tilt: how far the dominant axis deviates from E-W (0 = pure E-W, 90 = pure N-S)
#'   - eigenvectors: the 3x3 rotation matrix
#'   - variance_ratio: proportion of variance along first component
#'   - lonc, lat_c: centre of the great circle fit
#'   - alpha: oblique Mercator alpha parameter (azimuth of the axis)
#'
pca_axis <- function(xy) {
  lon <- xy[, 1] * pi / 180
  lat <- xy[, 2] * pi / 180

  ## Convert to 3D cartesian on unit sphere
  x <- cos(lat) * cos(lon)
  y <- cos(lat) * sin(lon)
  z <- sin(lat)

  xyz <- cbind(x, y, z)

  ## Centre and do PCA
  xyz_c <- scale(xyz, scale = FALSE)
  sv <- svd(xyz_c)

  ## Variance explained by each component
  var_total <- sum(sv$d^2)
  var_ratio <- sv$d^2 / var_total

  ## The first principal component direction (in 3D cartesian)
  pc1 <- sv$v[, 1]

  ## Convert PC1 direction back to lon/lat (as a point on the sphere)
  ## This gives the "pole" of the great circle perpendicular to the data spread
  ## Actually PC1 is the direction of maximum variance, so it IS the data axis
  pc1_lon <- atan2(pc1[2], pc1[1]) * 180 / pi
  pc1_lat <- asin(pc1[3]) * 180 / pi

  ## Centroid in lonlat
  centroid_xyz <- colMeans(xyz)
  centroid_xyz <- centroid_xyz / sqrt(sum(centroid_xyz^2))  # normalize
  lonc <- atan2(centroid_xyz[2], centroid_xyz[1]) * 180 / pi
  latc <- asin(centroid_xyz[3]) * 180 / pi

  ## Azimuth: bearing from centroid toward the PC1 direction
  ## We compute initial bearing from centroid to pc1 point
  azimuth <- initial_bearing(lonc, latc, pc1_lon, pc1_lat)

  ## Tilt: 0 = axis runs E-W, 90 = axis runs N-S
  tilt <- abs(90 - abs(azimuth %% 180))

  list(
    azimuth = azimuth,
    tilt = tilt,
    variance_ratio = var_ratio,
    pc1_lonlat = c(pc1_lon, pc1_lat),
    lonc = lonc,
    latc = latc,
    ## omerc parameters
    alpha = azimuth
  )
}


#' Initial bearing (forward azimuth) from point 1 to point 2
#' All inputs in degrees, output in degrees [0, 360)
initial_bearing <- function(lon1, lat1, lon2, lat2) {
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  dlon <- lon2 - lon1
  x <- sin(dlon) * cos(lat2)
  y <- cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
  bearing <- atan2(x, y) * 180 / pi
  (bearing + 360) %% 360
}


# =============================================================================
# proj(): generate a PROJ string from params + projection family/property
# =============================================================================

#' Build a PROJ string from params and a projection family or property
#'
#' @param p a params object
#' @param family projection family name (e.g. "laea", "stere", "omerc", "lcc", "ortho")
#'        or a property: "equal_area", "equidistant", "conformal", "auto"
#'
proj <- function(p, family = "auto") {
  lon_0 <- round(p@centre["lon_0"], 6)
  lat_0 <- round(p@centre["lat_0"], 6)

  ## Map properties to families based on data geometry

  if (family == "auto") {
    family <- auto_family(p)
  } else if (family == "equal_area") {
    family <- auto_family(p, property = "equal_area")
  } else if (family == "equidistant") {
    family <- auto_family(p, property = "equidistant")
  } else if (family == "conformal") {
    family <- auto_family(p, property = "conformal")
  }

  ## Build the PROJ string
  base <- sprintf("+proj=%s +lon_0=%s +lat_0=%s", family, lon_0, lat_0)

  ## Add family-specific parameters
  if (family == "lcc") {
    if (length(p@secant) == 2) {
      base <- sprintf("%s +lat_1=%s +lat_2=%s", base, p@secant["lat_1"], p@secant["lat_2"])
    } else {
      ## default secants at 1/6 and 5/6 of the latitude range
      lat_range <- p@extent_ll[c("ymin", "ymax")]
      lat_1 <- lat_range[1] + diff(lat_range) / 6
      lat_2 <- lat_range[2] - diff(lat_range) / 6
      base <- sprintf("%s +lat_1=%s +lat_2=%s", base, round(lat_1, 2), round(lat_2, 2))
    }
  }

  if (family == "omerc") {
    alpha <- p@axis$alpha
    if (is.null(alpha) || !is.finite(alpha)) alpha <- 0
    ## omerc uses lonc not lon_0
    base <- sprintf("+proj=omerc +lonc=%s +lat_0=%s +alpha=%s +gamma=0",
                    lon_0, lat_0, round(alpha, 2))
  }

  if (family == "stere") {
    ## add lat_ts (true scale latitude) for polar cases
    if (grepl("polar", p@hemisphere)) {
      lat_ts <- if (p@hemisphere == "polar_south") -71 else 71
      base <- sprintf("%s +lat_ts=%s", base, lat_ts)
    }
  }

  base
}


#' Choose a projection family based on data geometry
#'
auto_family <- function(p, property = NULL) {

  is_polar <- grepl("polar", p@hemisphere)
  is_mid <- grepl("mid", p@hemisphere)
  is_small <- p@span["lon_span"] < 30 && p@span["lat_span"] < 30

  ## If a specific property is requested, pick the best family for it
  if (!is.null(property)) {
    return(switch(property,
                  equal_area = if (is_polar || is_small) "laea" else if (is_mid && p@aspect == "wide") "aea" else "laea",
                  equidistant = if (is_polar) "aeqd" else if (is_mid && p@aspect == "wide") "eqdc" else "aeqd",
                  conformal = if (is_polar) "stere" else if (is_mid && p@aspect == "wide") "lcc" else "stere",
                  "laea"  # fallback
    ))
  }

  ## Auto mode: pick based on geometry
  ## Small region → stereographic (conformal, good for local work)
  if (is_small) return("stere")

  ## Polar → azimuthal
  if (is_polar) return("stere")

  ## Mid-latitude wide band → conic
  if (is_mid && p@aspect == "wide") return("lcc")

  ## Elongated along a great circle → oblique Mercator
  if (!is.null(p@axis$variance_ratio) &&
      length(p@axis$variance_ratio) >= 2 &&
      p@axis$variance_ratio[1] > 0.85 &&
      p@axis$tilt > 20) {
    return("omerc")
  }

  ## Default: Lambert Azimuthal Equal-Area (good all-rounder)
  "laea"
}


# =============================================================================
# gproj_region(): quick projection from a point + radius
# =============================================================================

#' Create a gproj object centred on a point with a given radius
#'
#' @param lon centre longitude
#' @param lat centre latitude
#' @param radius_km radius in kilometres
#' @param property "equal_area", "conformal", "equidistant", or a PROJ family name
#' @return a list with gproj object, extent in projected coords, and grid_spec
#'
gproj_region <- function(lon, lat, radius_km = 500, property = "equal_area") {
  ## Create params from a single point (with hemisphere info from lat)
  ## We fake a small extent around the point to get hemisphere classification
  fake_xy <- cbind(
    lon + c(-1, 1, 0, 0),
    lat + c(0, 0, -1, 1)
  )
  p <- params(fake_xy)
  ## Override centre with exact values
  p@centre <- c(lon_0 = lon, lat_0 = lat)

  ## Get CRS
  crs_string <- proj(p, property)
  g <- gproj(crs_string)

  ## Compute extent in projected coordinates
  radius_m <- radius_km * 1000
  extent_xy <- c(
    xmin = -radius_m, xmax = radius_m,
    ymin = -radius_m, ymax = radius_m
  )

  list(
    gproj = g,
    crs = crs_string,
    extent = extent_xy,
    centre_xy = c(x = 0, y = 0),
    radius_m = radius_m
  )
}


#' Convenience: generate a grid specification from a gproj_region result
#'
#' @param region result from gproj_region()
#' @param resolution_m grid cell size in metres
#' @return list with crs, extent, dimension, resolution
#'
grid_spec <- function(region, resolution_m = 1000) {
  ex <- region$extent
  nx <- ceiling((ex["xmax"] - ex["xmin"]) / resolution_m)
  ny <- ceiling((ex["ymax"] - ex["ymin"]) / resolution_m)
  list(
    crs = region$crs,
    extent = ex,
    dimension = c(ncol = nx, nrow = ny),
    resolution = c(dx = resolution_m, dy = resolution_m)
  )
}


# =============================================================================
# Mercator devolution detection (foundation for Phase 2 detect_grid)
# =============================================================================

#' Test whether a 1D latitude sequence follows Mercator spacing
#'
#' Given a vector of latitude values (e.g. from a rectilinear grid),
#' check if the spacing matches the Mercator y formula.
#' Returns a list with:
#'   - is_mercator: logical
#'   - residual: max deviation from perfect Mercator spacing
#'   - merc_resolution: the resolution in Mercator y-units (if detected)
#'
#' @param lat_values numeric vector of latitude values (sorted)
#' @param tol_relative relative tolerance for spacing uniformity (default 0.001)
#'
test_mercator_spacing <- function(lat_values, tol_relative = 0.001) {
  lat_values <- sort(unique(lat_values))
  n <- length(lat_values)
  if (n < 3) return(list(is_mercator = FALSE, reason = "too few values"))

  ## Convert to Mercator y (on unit sphere, then scale doesn't matter for uniformity test)
  ## y = ln(tan(pi/4 + phi/2))  [Web Mercator / Spherical Mercator]
  phi <- lat_values * pi / 180
  merc_y <- log(tan(pi/4 + phi/2))

  ## Check uniformity of Mercator y spacing
  dy <- diff(merc_y)
  median_dy <- median(dy)

  if (abs(median_dy) < 1e-12) {
    return(list(is_mercator = FALSE, reason = "zero spacing"))
  }

  ## Relative deviation from uniform spacing
  rel_dev <- abs(dy - median_dy) / abs(median_dy)
  max_dev <- max(rel_dev)

  is_merc <- max_dev < tol_relative

  ## If Mercator, compute the resolution in metres (assuming WGS84 sphere approx)
  merc_res <- NULL
  if (is_merc) {
    ## Earth radius ~6378137 for WGS84
    R <- 6378137
    merc_res <- abs(median_dy) * R
  }

  list(
    is_mercator = is_merc,
    max_relative_deviation = max_dev,
    merc_resolution_m = merc_res,
    n_values = n
  )
}

#' Test whether a 1D coordinate sequence is uniformly spaced (regular)
#'
#' @param values numeric vector
#' @param tol_relative relative tolerance
#' @return list with is_regular, resolution, max_deviation
#'
test_regular_spacing <- function(values, tol_relative = 0.001) {
  values <- sort(unique(values))
  n <- length(values)
  if (n < 3) return(list(is_regular = FALSE, reason = "too few values"))

  dv <- diff(values)
  median_dv <- median(dv)

  if (abs(median_dv) < 1e-12) {
    return(list(is_regular = FALSE, reason = "zero spacing"))
  }

  rel_dev <- abs(dv - median_dv) / abs(median_dv)
  max_dev <- max(rel_dev)

  list(
    is_regular = max_dev < tol_relative,
    resolution = median_dv,
    max_relative_deviation = max_dev,
    n_values = n
  )
}


# =============================================================================
# detect_rectilinear(): test if 2D lon/lat arrays are rectilinear
# =============================================================================

#' Test if 2D coordinate arrays are rectilinear (lon varies only along cols,
#' lat varies only along rows)
#'
#' @param lon 2D matrix of longitudes (nrow x ncol)
#' @param lat 2D matrix of latitudes (nrow x ncol)
#' @param tol tolerance for "constant along an axis" test (in degrees)
#' @return list describing the rectilinear structure
#'
detect_rectilinear <- function(lon, lat, tol = 0.01) {
  stopifnot(is.matrix(lon), is.matrix(lat))
  stopifnot(identical(dim(lon), dim(lat)))

  nr <- nrow(lon)
  nc <- ncol(lon)

  ## Test: does lon vary only along columns (i.e. each row has the same lon sequence)?
  ## Check SD of each column of lon — should be near-zero
  lon_col_sd <- apply(lon, 2, sd, na.rm = TRUE)
  lon_row_varies <- max(lon_col_sd, na.rm = TRUE) < tol

  ## Alternatively: does lon vary only along rows? (transposed convention)
  lon_row_sd <- apply(lon, 1, sd, na.rm = TRUE)
  lon_col_varies <- max(lon_row_sd, na.rm = TRUE) < tol

  ## Same for lat

  lat_row_sd <- apply(lat, 1, sd, na.rm = TRUE)
  lat_row_const <- max(lat_row_sd, na.rm = TRUE) < tol

  lat_col_sd <- apply(lat, 2, sd, na.rm = TRUE)
  lat_col_const <- max(lat_col_sd, na.rm = TRUE) < tol

  ## Standard orientation: lon varies along columns (x-axis), lat along rows (y-axis)
  ## i.e. lon is constant within each row, lat is constant within each column
  ## Wait — NetCDF convention: lon[y,x] so lon should be constant within each column
  ## and lat constant within each row? No...
  ## Actually the standard case: lon[i,j] depends only on j, lat[i,j] depends only on i
  ## So: for lon, within each row (fixed i), values change with j → lon varies along cols
  ##     for lat, within each column (fixed j), values change with i → lat varies along rows
  ## Test: lon constant within columns means lon[,j] has low SD → lon_col_sd is low? No,
  ## lon[,j] has all the values for column j across all rows — if lon depends only on j
  ## then lon[,j] should all be the same value → sd ≈ 0 → that's lon_col_sd

  is_rectilinear <- lon_col_varies && lat_col_const  ## lon[i,j]~f(i), lat[i,j]~g(j) — transposed
  orientation <- "standard"

  if (!is_rectilinear) {
    ## Try the other orientation
    is_rectilinear <- lon_row_varies && lat_row_const
    ## Actually let me be more careful. The two cases:
    ## Case A: lon depends on column index (j), lat depends on row index (i)
    ##   → lon constant along rows of the lon matrix (each col of lon has same value)
    ##   → lat constant along cols of the lat matrix (each row of lat has same value)
    ##   → lon_col_sd ≈ 0 and lat_row_sd ≈ 0
    is_rectilinear_A <- all(lon_col_sd < tol) && all(lat_row_sd < tol)

    ## Case B: transposed — lon depends on row index, lat depends on column index
    is_rectilinear_B <- all(lon_row_sd < tol) && all(lat_col_sd < tol)

    is_rectilinear <- is_rectilinear_A || is_rectilinear_B
    orientation <- if (is_rectilinear_A) "standard" else if (is_rectilinear_B) "transposed" else "none"
  }

  if (!is_rectilinear) {
    return(list(
      is_rectilinear = FALSE,
      reason = "coordinates vary in both dimensions"
    ))
  }

  ## Extract marginal coordinate sequences
  if (orientation == "standard") {
    lon_seq <- lon[1, ]    # lon varies along columns
    lat_seq <- lat[, 1]    # lat varies along rows
  } else {
    lon_seq <- lon[, 1]    # lon varies along rows
    lat_seq <- lat[1, ]    # lat varies along columns
  }

  ## Test regularity of each marginal
  lon_reg <- test_regular_spacing(lon_seq)
  lat_reg <- test_regular_spacing(lat_seq)

  ## Test Mercator hypothesis on latitude
  lat_merc <- test_mercator_spacing(lat_seq)

  ## Classify
  type <- if (lon_reg$is_regular && lat_reg$is_regular) {
    "regular_lonlat"
  } else if (lon_reg$is_regular && lat_merc$is_mercator) {
    "mercator_devolved"
  } else {
    "rectilinear_nonuniform"
  }

  list(
    is_rectilinear = TRUE,
    orientation = orientation,
    type = type,
    dimensions = dim(lon),
    lon_marginal = lon_seq,
    lat_marginal = lat_seq,
    lon_regular = lon_reg,
    lat_regular = lat_reg,
    lat_mercator = lat_merc,
    ## If Mercator, provide the reconstruction
    mercator_crs = if (type == "mercator_devolved") {
      sprintf("+proj=merc +lon_0=%s +datum=WGS84",
              round(mean(range(lon_seq)), 2))
    } else NULL,
    mercator_resolution = if (type == "mercator_devolved") {
      lat_merc$merc_resolution_m
    } else NULL
  )
}


# =============================================================================
# Printing / show methods
# =============================================================================

method(format, gproj) <- function(x, ...) {
  if (nchar(x@crs) == 0) {
    "<gproj: identity (no projection)>"
  } else {
    sprintf("<gproj: %s>", x@crs)
  }
}

method(print, gproj) <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}

method(format, params) <- function(x, ...) {
  lines <- c(
    sprintf("<params: %d points>", x@n),
    sprintf("  centre:     %.4f, %.4f", x@centre[1], x@centre[2]),
    sprintf("  extent:     [%.4f, %.4f] x [%.4f, %.4f]",
            x@extent_ll[1], x@extent_ll[2], x@extent_ll[3], x@extent_ll[4]),
    sprintf("  span:       %.2f° lon x %.2f° lat", x@span[1], x@span[2]),
    sprintf("  hemisphere: %s", x@hemisphere),
    sprintf("  aspect:     %s", x@aspect)
  )
  if (!is.null(x@axis$azimuth) && is.finite(x@axis$azimuth)) {
    lines <- c(lines,
               sprintf("  axis:       azimuth=%.1f°, tilt=%.1f°, var_ratio=%.2f/%.2f/%.2f",
                       x@axis$azimuth, x@axis$tilt,
                       x@axis$variance_ratio[1],
                       x@axis$variance_ratio[2],
                       x@axis$variance_ratio[3])
    )
  }
  if (length(x@secant) == 2) {
    lines <- c(lines,
               sprintf("  secant:     lat_1=%.2f, lat_2=%.2f", x@secant[1], x@secant[2])
    )
  }
  paste(lines, collapse = "\n")
}

method(print, params) <- function(x, ...) {
  cat(format(x), "\n")
  invisible(x)
}
