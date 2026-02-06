# gproj <- S7::new_class(name = "gproj", package = "gproj",
#   properties = list(
#     crs = S7::new_property(class = S7::class_character, default = "EPSG:4326"),
#     proj_xy = S7::new_property(
#       class = S7::class_function, getter = function(self) {mk_proj_xy(self@crs)}),
#     proj_ll = S7::new_property(
#       class = S7::class_function, getter = function(self) {mk_proj_ll(self@crs)}
#     ))
# )
#
#
# mk_proj_xy <- function(crs) {
#   function(x) {
#     reproj::reproj_xy(x, crs, source = "EPSG:4326")
#  }
# }
#
# mk_proj_xyz <- function(crs) {
#   function(x) {
#     reproj::reproj_xyz(x, "+proj=geocent", source = "EPSG:4326")
#   }
# }
#
# mk_proj_llz <- function(crs)  {
#   function(x) {
#   reproj::reproj_xyz(x, "EPSG:4326", source = "+proj=geocent")
#   }
# }
#
# mk_proj_ll <- function(crs) {
#   function(x) {
#     reproj::reproj_xy(x, "EPSG:4326", source = crs)
#   }
# }
# area <- function(x) {
#     area::polygon_area(x)
# }
#
# mid <- function(x) {
#     matrix(apply(x, 2, mean), ncol = 2L)
# }
# proj_xyz <- function(x) {
#     reproj::reproj_xyz(x, "+proj=geocent", source = "EPSG:4326")
# }
# proj_llz <- function(x) {
#     reproj::reproj_xyz(x, "EPSG:4326", source = "+proj=geocent")
# }
# midll <- function(p) {
#     proj_llz(matrix(apply(proj_xyz(p), 2L, mean), ncol = 3L))[,1:2, drop = F]
# }
# ## treat these as a closure for now
# params <- function(x, secant = FALSE) {
#     lon <- as.integer(midll(x)[1])
#     lat <- as.integer(midll(x)[2])
#     if (secant) {
#         sec <-  as.integer(range(x[,2]))
#         return(list(lon_0 = lon, lat_0 = lat, lat_1 = sec[1], lat_2 = sec[2]))
#     }
#     list(lon_0 = lon, lat_0 = lat)
# }
# params_fold <- function(x) {
#     unname(apply(cbind(names(x), unlist(x)), 1, paste0, collapse = "="))
# }
#
# proj <- function(params, family   = NULL) {
#     if (is.null(family)) stop("family must be named (i.e. laea, stere, lcc")
#     params$proj <- family
#     nms <- setdiff(names(params), "proj")
#     params <- params[c("proj", nms)]
#
#  paste0( sprintf("+%s", params_fold(params)), collapse = " ")
# }
# # # alpha <- seq(-50, 50, by = 5)
# # # i <- 0
# # # i <- i + 1
# # prj <- sprintf("+proj=omerc +lonc=-113 +lat_0=-20 +alpha=%f +gamma=1", 15)
# # proj_xy <- function(x) {
# #     reproj::reproj_xy(x, prj, source = "EPSG:4326")
# # }
# # proj_ll <- function(x) {
# #     reproj::reproj_xy(x, "EPSG:4326", source = prj)
# # }
