# gproj design notes

## Identity

gproj is a **projection concierge**: given coordinates or a region of interest,
it tells you what projection to use and gives you the transform. It's not
spatial infrastructure — it sits between "I have lon/lat data" and "I need a
sensible projected view".

## Dependencies

- **S7**: class system
- **reproj**: coordinate transformation (wraps PROJ)
- **geographiclib** (optional, Phase 3): geodesic area/distance

## Phase 1: Core (current)

### Classes

**`gproj(crs)`** — S7 object storing a CRS string, with `@proj_xy` (forward)
and `@proj_ll` (inverse) methods. Identity when crs is empty.

**`params(xy)` / `params(lon, lat)`** — Characterises the geometry of a
coordinate set:
- `@centre` — lon_0, lat_0 (column mean; TODO: circular mean for antimeridian)
- `@extent_ll` — bounding box in longlat
- `@span` — angular extent in degrees
- `@hemisphere` — polar_south/north, mid_south/north, equatorial
- `@aspect` — wide/tall/square (lon_span vs lat_span ratio)
- `@axis` — PCA on 3D cartesian: azimuth, tilt, variance_ratio
- `@secant` — lat_1/lat_2 for conic projections

### Functions

**`proj(params, family)`** — Generates a PROJ string. `family` can be:
- A PROJ name: "laea", "stere", "omerc", "lcc", "ortho", "aeqd", etc.
- A property: "equal_area", "equidistant", "conformal"
- "auto" — picks family from data geometry (the textbook decision tree)

**`gproj_region(lon, lat, radius_km, property)`** — Quick entry point:
"equal area centred here, 500km radius". Returns gproj + extent.

**`grid_spec(region, resolution_m)`** — Target grid specification from a region.

### Mercator devolution detection (foundation for Phase 2)

**`test_mercator_spacing(lat_values)`** — Tests if 1D latitude sequence follows
Mercator y-spacing: `y = R * ln(tan(π/4 + φ/2))`. Converts lat to Mercator y,
checks if diff(merc_y) is constant.

**`test_regular_spacing(values)`** — Tests if a 1D coordinate sequence has
uniform spacing.

**`detect_rectilinear(lon_2d, lat_2d)`** — Tests if 2D coordinate arrays are
rectilinear (lon varies only along one axis, lat along the other). Classifies as:
- `regular_lonlat` — uniform spacing in both lon and lat
- `mercator_devolved` — uniform lon, Mercator-spaced lat
- `rectilinear_nonuniform` — rectilinear but neither regular nor Mercator

## Phase 2: detect_grid (planned)

The full grid detection pipeline for "cryptic curvilinear" data.

### Algorithm

1. **Rectilinear test** (from Phase 1): check if lon/lat arrays are separable
2. **Candidate CRS generation**: use params() centre + hemisphere to pick candidates
3. **Regularity test in projected space**:
   - Project lon/lat to candidate CRS (subsample: edges + cross suffice)
   - For each row: is diff(x) constant?
   - For each column: is diff(y) constant?
   - Are row-y values constant (rows truly horizontal)?
4. **Resolution snapping**: 24999.84 → 25000
5. **Half-cell extent recovery**: cell centres + resolution → true grid extent
6. **EPSG matching**: compare detected CRS + resolution against known grid database

### The "centre trick"

Even without a candidate database, project to `+proj=stere +lat_0={detected}
+lon_0={detected}` (or lcc, etc. based on hemisphere) and test regularity.
If regular, you've found the family. Then search EPSG for a match.

### Output: detect_grid(lon, lat) returns

```r
result@type        # "regular_projected" | "rectilinear_lonlat" | "mercator_devolved" | "curvilinear"
result@crs         # "EPSG:3412" or PROJ string
result@resolution  # c(25000, 25000)
result@extent      # c(-3950000, 3950000, -3950000, 4350000)
result@dimensions  # c(316, 332)
result@confidence  # 0-1
result@residuals   # max deviation from perfect grid (CRS units)
result@gproj       # ready-to-use gproj object
result@grid_spec   # list(crs, extent, res, dim)
result@gdal_args   # c("-a_srs", "EPSG:3412", "-a_ullr", ...) for immediate use
```

### Known grid database (internal data)

Common gridded product CRSs with typical resolutions:

| CRS | Res (m) | Domain | Products |
|-----|---------|--------|----------|
| EPSG:3412 | 25000 | Antarctic sea ice | NSIDC Bootstrap/NASA Team |
| EPSG:3413 | 25000 | Arctic sea ice | NSIDC |
| EPSG:3031 | various | Antarctic | BAS, AAD, many |
| EPSG:3976 | 25000 | Antarctic sea ice | NSIDC alt |
| EPSG:3995 | various | Arctic | various |
| +proj=merc | various | Global ocean | AVISO, GHRSST, Copernicus |
| +proj=lcc (Euro) | 12000-50000 | Europe | CORDEX, EURO-CORDEX, WRF |
| +proj=stere +lat_0=90 | various | Arctic | ERA5 reduced |

## Phase 3: Geodesic operations (planned)

If geographiclib is available:
- `geodesic_area(lonlat_polygon)` — exact ellipsoidal area
- `geodesic_distance(lon1, lat1, lon2, lat2)` — exact ellipsoidal distance
- `geodesic_line(lon1, lat1, lon2, lat2, n)` — points along a geodesic

Otherwise fallback to spherical approximations.

## Key insight: the entropy problem

A "cryptic curvilinear" grid stores N×M×2 coordinate values that can be
replaced by 4 numbers (extent) + 1 string (CRS). That's a compression ratio
of ~35,000:1 for a 316×332 grid, and it's *lossless* — better than the
original because the lon/lat arrays introduce floating-point noise.

The goal of detect_grid() is to discover this latent structure automatically.

## Antimeridian handling (TODO)

params() currently uses simple column mean for centre. Needs circular mean
for datasets spanning ±180. Detection: if max(lon) - min(lon) > 350 and there
are values near both -180 and 180, shift to 0-360 before computing centre.
