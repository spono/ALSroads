#' Computation of the least cost path
#'
#' The function computes the conductivity layers, compute the global conductivity, mask the waterbodies,
#' add ending caps, generates location of point A and B, computes the least cost path and returns
#' the new centerline geometry.
#'
#' @param las LAS. the point cloud extracted with a buffer from the centerline
#' @param centerline sf, the road
#' @param DTM,water DTM and water body mask
#' @param list. parameters
#'
#' @noRd
least_cost_path = function(las, centerline, dtm, conductivity, water, param)
{
  display <- getOption("ALSroads.debug.finding")

  if (is.null(dtm) && is.null(conductivity)) stop("'dtm' and 'conductivity' cannot be both NULL", call. = FALSE)
  if (!is.null(dtm) && !is.null(conductivity) && conductivity != "v2") stop("'dtm' or 'conductivity' must be NULL", call. = FALSE)

  # Compute the limit polygon where we actually need to work
  # This allows to crop the raster (DTM or conductivty) and extract
  # only the part we are working with
  hold <- sf::st_buffer(centerline, param[["extraction"]][["road_buffer"]] + 2)

  # Clip and mask the conductivity (if provided)
  # If the conductivity has a resolution more than 2 m, aggregate to 2 m.
  # If the resolution is less than 2 meter, abort.
  if (inherits(conductivity, "RasterLayer"))
  {
    conductivity <- raster::crop(conductivity, hold)
    conductivity <- raster::mask(conductivity, hold)

    res <- round(raster::res(conductivity)[1], 2)

    if (res > 2)
      stop("The conductivity must have a resolution of 1 m or less.")

    if (res < 2)
      conductivity <- raster::aggregate(conductivity, fact = 2/res, fun = mean, na.rm = TRUE)
  }

  if (is.character(conductivity) && conductivity == "v2")
  {
    verbose("Computing conductivity maps...\n")
    conductivity <- rasterize_conductivity2(las, dtm = dtm, param = param, return_all = FALSE, return_stack = FALSE)
  }

  # If no conductivity is provided it means that the DTM is. We compute the conductivity
  if (is.null(conductivity))
  {
    verbose("Computing conductivity maps...\n")
    conductivity <- rasterize_conductivity(las, dtm = dtm, param = param, return_all = FALSE, return_stack = FALSE)
  }

  # Handle bridge case:
  # If a road crosses a water body it builds a bridge, i.e. a polygon in which we will
  # force a conductivity of 1 later
  bridge = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    water <- sf::st_crop(water, bbox)
    bridge <- sf::st_intersection(sf::st_geometry(centerline), water)
    if (length(bridge) > 0) bridge <- sf::st_buffer(bridge, 5)
    if (length(water)  > 0) conductivity <- raster::mask(conductivity, sf::as_Spatial(water), inverse = TRUE)
  }

  if (length(bridge) > 0)
  {
    cells <- raster::cellFromPolygon(conductivity, sf::as_Spatial(bridge))
    cells <- unlist(cells)
    tmp   <- conductivity
    tmp[cells] <- 1
    conductivity <- tmp

    if (display) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 1m with bridge")
  }

  # Compute low resolution conductivity with mask
  conductivity <- mask_conductivity(conductivity, centerline, param)

  # Compute start and end points
  AB <- start_end_points(centerline, param)
  A  <- AB$A
  B  <- AB$B

  # Compute the transition
  trans <- transition(conductivity)

  # Find the path
  verbose("Computing least cost path...\n")
  path <- find_path(trans, centerline, A, B, param)

  return(path)
}

#' Mask the conductivity and modify some pixels
#'
#' Update some pixels by masking the map with the bounding polygon, multiplying by a distance to
#' road cost factor (not described in the paper) and add terminal caps conductive pixels
#'  to allow driving from the point A and B further apart the road
#' @noRd
mask_conductivity <- function(conductivity, centerline, param)
{
  verbose("Computing conductivity masks...\n")

  # Boundary masking
  hull <- sf::st_buffer(centerline, param$extraction$road_buffer)

  # Fix an issue for road at the very edge of a catalog
  bb_hull <- sf::st_bbox(hull)
  bb_cond <- sf::st_bbox(conductivity)
  if (bb_hull[1] < bb_cond[1] | bb_hull[2] < bb_cond[2] | bb_hull[3] > bb_cond[3] | bb_hull[4] > bb_cond[4])
    conductivity = raster::extend(conductivity, bb_hull)

  conductivity <- raster::mask(conductivity, hull)

  if (getOption("ALSroads.debug.finding")) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity 2m")
  verbose("   - Masking\n")

  # Penalty factor based on distance-to-road.
  # Actually small penalty since this part was not very successful. Not described in the paper.
  # Could maybe be removed. Yet it might be useful if putting more constraints
  p <- sf::st_buffer(centerline, 1)
  f <- fasterize::fasterize(p, conductivity)
  f <- raster::distance(f)
  f <- raster::mask(f, hull)
  fmin <- min(f[], na.rm = T)
  fmax <- max(f[], na.rm = T)
  target_min <- 1-param[["constraint"]][["confidence"]]
  f <- (1-(((f - fmin) * (1 - target_min)) / (fmax - fmin)))
  conductivity <- f*conductivity

  if (getOption("ALSroads.debug.finding")) raster::plot(f, col = viridis::viridis(25), main = "Distance factor")
  verbose("   - Road rasterization and distance factor map\n")

  # Set a conductivity of 1 in the caps and 0 on the outer half ring link in figure 6
  # We could use raster::cellFromPolygon but it is slow. This workaround using lidR is complex but fast.
  # Maybe using terra we could simplify the code.
  caps <- make_caps(centerline, param)
  xy <- raster::xyFromCell(conductivity, 1: raster::ncell(conductivity))
  xy <- as.data.frame(xy)
  xy$z <- 0
  names(xy) <- c("X", "Y", "Z")
  xy <- lidR::LAS(xy, lidR::LASheader(xy))
  xy@header$`Global Encoding`$WKT = TRUE
  lidR::st_crs(xy) <- sf::st_crs(caps$caps)

  res <- !is.na(lidR:::point_in_polygons(xy, caps$caps))
  conductivity[res] <- 1
  res <- !is.na(lidR:::point_in_polygons(xy, caps$shields))
  conductivity[res] <- 0
  conductivity[is.nan(conductivity)] <- NA

  if (getOption("ALSroads.debug.finding")) raster::plot(conductivity, col = viridis::inferno(15), main = "Conductivity with end caps")
  verbose("   - Add full conductivity end blocks\n")

  return(conductivity)
}

start_end_points = function(centerline, param)
{
  caps <- make_caps(centerline, param)$caps
  #P <- sf::st_cast(caps, "POLYGON")
  C <- sf::st_centroid(caps)
  A <- sf::st_coordinates(C[1])
  B <- sf::st_coordinates(C[2])
  return(list(A = A, B = B))
}

#' @importClassesFrom  gdistance TransitionLayer
transition <- function(conductivity, directions = 8, geocorrection = TRUE)
{
  verbose("Computing graph map...\n")

  x = conductivity
  use_terra = FALSE

  if (methods::is(x, "RasterLayer"))
  {
    bb = raster::extent(x)
    ncells = raster::ncell(x)
    val = raster::values(x)
  }
  else if (methods::is(x, "SpatRaster"))
  {
    bb = terra::ext(x)
    bb = raster::extent(bb[1:4])
    ncells = terra::ncell(x)
    val = terra::values(x)[,1]
    use_terra = TRUE
  }
  else
  {
    stop("Only RasterLayer and SpatRaster are supported")
  }

  symm = TRUE

  tr <- methods::new("TransitionLayer",
            nrows=as.integer(nrow(x)),
            ncols=as.integer(ncol(x)),
            extent=bb,
            crs=sp::CRS(),
            transitionMatrix = Matrix::Matrix(0, ncells,ncells),
            transitionCells = 1:ncells)

  transitionMatr <- gdistance::transitionMatrix(tr)
  Cells <- which(!is.na(val))

  if (use_terra)
    adj <- terra::adjacent(x, cells=Cells, pairs=TRUE, directions=directions)
  else
    adj = raster::adjacent(x, cells=Cells, pairs=TRUE, directions=directions)

  if(symm)
    adj <- adj[adj[,1] < adj[,2],]

  dataVals <- cbind(val[adj[,1]], val[adj[,2]])
  transition.values <- rowMeans(dataVals)
  transition.values[is.na(transition.values)] <- 0

  if(!all(transition.values >= 0))
    warning("transition function gives negative values")

  transitionMatr[adj] <- as.vector(transition.values)
  if(symm)
    transitionMatr <- Matrix::forceSymmetric(transitionMatr)

  gdistance::transitionMatrix(tr) <- transitionMatr
  gdistance::matrixValues(tr) <- "conductance"


  #trans <- gdistance::transition(conductivity, transitionFunction = mean, directions = 8)
  trans = tr
  verbose("   - Transition graph\n")

  if (isTRUE(geocorrection))
  {
    trans <- gdistance::geoCorrection(trans)
    verbose("   - Geocorrection graph\n")
  }

  return(trans)
}


find_path = function(trans, centerline, A, B, param)
{
  caps <- make_caps(centerline, param)$caps
  caps <- sf::st_union(caps)
  trans@crs <- methods::as(sf::NA_crs_, "CRS") # workaround to get rid of rgdal warning

  cost <- gdistance::costDistance(trans, A, B)

  if (is.infinite(cost))
  {
    verbose("    - Impossible to reach the end of the road\n")
    path <- sf::st_geometry(centerline)
    path <- sf::st_as_sf(path)
    path$CONDUCTIVITY <- 0
    return(path)
  }

  path <- gdistance::shortestPath(trans, A, B, output = "SpatialLines") |> suppressWarnings()
  path <- sf::st_as_sf(path)
  len  <- sf::st_length(path)
  path <- sf::st_simplify(path, dTolerance = 3)
  path <- sf::st_set_crs(path, sf::NA_crs_)
  path <- sf::st_set_crs(path, sf::st_crs(centerline))
  path <- sf::st_difference(path, caps)
  path$CONDUCTIVITY <- round(as.numeric(len/cost),2)

  if (getOption("ALSroads.debug.finding")) plot(sf::st_geometry(path), col = "red", add = T, lwd = 2)

  return(path)
}

sobel <- function(img, ker = 3) UseMethod("sobel", img)

sobel.RasterLayer <- function(img, ker = 3)
{
  slop <- raster::as.matrix(img)
  img[] <- sobel.matrix(slop, ker)
  img
}

sobel.matrix <- function(img, ker = 3)
{
  # define horizontal and vertical Sobel kernel
  if (ker == 3)
    Shoriz <- matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)

  if (ker == 5)
  {
    A = 2*sqrt(2)
    B = sqrt(5)
    C = sqrt(2)
    Shoriz <- matrix(c(2,2,4,2,2,
                       1,1,2,1,1,
                       0,0,0,0,0,
                       -1,-1,-2,-1,-1,
                       -2,-2,-4,-2,-2), nrow = 5, byrow = TRUE)
  }
  Svert <- t(Shoriz)
  nas <- is.na(img)
  img[nas] <- 0

  # get horizontal and vertical edges
  imgH <- EBImage::filter2(img, Shoriz)
  imgV <- EBImage::filter2(img, Svert)

  # combine edge pixel data to get overall edge data
  hdata <- EBImage::imageData(imgH)
  vdata <- EBImage::imageData(imgV)
  edata <- sqrt(hdata^2 + vdata^2)
  edata[nas] <- NA
  edata
}

make_caps <- function(centerline, param)
{

  buf <- param[["extraction"]][["road_buffer"]]

  len <- as.numeric(20/sf::st_length(centerline))
  start <- lwgeom::st_linesubstring(centerline, 0, len)
  end <- lwgeom::st_linesubstring(centerline, 1-len, 1)
  s1 <- lwgeom::st_startpoint(start)
  e1 <- lwgeom::st_endpoint(start)
  s2 <- lwgeom::st_startpoint(end)
  e2 <- lwgeom::st_endpoint(end)
  l1 <- sf::st_coordinates(c(s1, e1))
  l2 <- sf::st_coordinates(c(s2, e2))

  start <- sf::st_sfc(sf::st_linestring(l1), crs = sf::st_crs(centerline))
  end <- sf::st_sfc(sf::st_linestring(l2), crs = sf::st_crs(centerline))

  poly1 <- sf::st_geometry(sf::st_buffer(start, buf, endCapStyle = "FLAT"))
  poly2 <- sf::st_geometry(sf::st_buffer(end,   buf, endCapStyle = "FLAT"))
  poly3 <- sf::st_geometry(sf::st_buffer(centerline, buf))

  A <- lwgeom::st_startpoint(centerline)
  B <- lwgeom::st_endpoint(centerline)

  cap_A <- sf::st_buffer(A, buf)
  cap_B <- sf::st_buffer(B, buf)

  shield_A <- sf::st_buffer(A, buf - 6)
  shield_B <- sf::st_buffer(B, buf - 6)
  shield <- sf::st_union(shield_A, shield_B)

  caps_A <- sf::st_difference(cap_A, poly1)
  caps_A <- sf::st_cast(caps_A, "POLYGON")
  if (length(caps_A) > 1)
    caps_A <- caps_A[which.max(sf::st_area(caps_A))]

  caps_B <- sf::st_difference(cap_B, poly2)
  caps_B <- sf::st_cast(caps_B, "POLYGON")
  if (length(caps_B) > 1)
    caps_B <- caps_B[which.max(sf::st_area(caps_B))]

  caps <- c(caps_A, caps_B)
  #caps <- sf::st_union(caps)
  shield <- sf::st_difference(sf::st_buffer(sf::st_union(caps, sf::st_union(poly1, poly2)), 0.01), shield)

  sf::st_crs(caps) <- sf::st_crs(centerline)
  sf::st_crs(shield) <- sf::st_crs(centerline)

  return(list(caps = caps, shields = shield))
}
