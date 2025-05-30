#' @export
#' @rdname rasterize_conductivity
rasterize_conductivity2 <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  UseMethod("rasterize_conductivity2", las)
}

#' @export
rasterize_conductivity2.LAS <- function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  use_intensity <- "Intensity" %in% names(las)
  display <- getOption("ALSroads.debug.finding")
  pkg <- if (is.null(dtm)) getOption("lidR.raster.default") else lidR:::raster_pkg(dtm)

  if (is.null(dtm))
  {
    dtm <- lidR::rasterize_terrain(las, 1, lidR::tin(), pkg = "raster")
  }
  else if (lidR:::is_raster(dtm))
  {
    res <- round(raster::res(dtm)[1], 2)
    if (res > 1) stop("The DTM must have a resolution of 1 m or less.")

    bb = lidR::st_bbox(las)
    dtm <- raster::crop(dtm, raster::extent(bb))

    if (res < 1)
      dtm <- raster::aggregate(dtm, fact = 1/res, fun = mean)
  }
  else
  {
    stop("dtm must be a RasterLayer or must be NULL")
  }

  # plot(dtm, col = gray(1:30/30))

  mask = NULL
  if (!is.null(water) && length(sf::st_geometry(water)) > 0)
  {
    id <- NULL
    water <- sf::st_geometry(water)
    bbox <- suppressWarnings(sf::st_bbox(las))
    bbox <- sf::st_set_crs(bbox, sf::st_crs(water))
    mask <- sf::st_crop(water, bbox)
    if (length(mask) > 0)
    {
      las = lidR::classify_poi(las, lidR::LASWATER, roi = water)
      mask <- sf::as_Spatial(mask)
    }
    else
    {
      mask = NULL
    }
  }

  # Force to use raster
  if (lidR:::raster_pkg(dtm) == "terra")
    dtm <- raster::raster(dtm)

  nlas <- lidR::normalize_height(las, dtm) |> suppressMessages() |> suppressWarnings()
  nlas@data[Classification == LASWATER & Z > 2, Classification := LASBRIGDE]
  bridge = lidR::filter_poi(nlas, Classification == LASBRIGDE)
  bridge = sf::st_coordinates(bridge, z = FALSE)

  # Terrain metrics using the raster package (slope, roughness)
  slope <- terra::terrain(dtm, opt = c("slope"), unit = "degrees")
  if (!is.null(mask)) slope <- raster::mask(slope, mask, inverse = T)
  if(display) raster::plot(slope, col = gray(1:30/30), main = "slope")

  smoothdtm = raster::focal(dtm, matrix(1,5,5), mean)
  roughdtm = dtm - smoothdtm
  roughness <- terra::terrain(roughdtm, opt = c("roughness"), unit = "degrees")

  if(display) {
    raster::plot(roughdtm, col = gray(1:30/30), main = "Residual roughness")
    raster::plot(roughness, col = gray(1:30/30), main = "Roughness")
  }

  # Slope-based conductivity
  s <- param$conductivity$s
  sigma_s <- activation(slope, s, "piecewise-linear", asc = FALSE)
  sigma_s <- raster::aggregate(sigma_s, fact = 2, fun = mean)

  if (display) raster::plot(sigma_s, col = viridis::viridis(25), main = "Conductivity slope")
  verbose("   - Slope conductivity map \n")

  # Roughness-based conductivity
  r <- param$conductivity$r
  sigma_r <- activation(roughness, r, "piecewise-linear", asc = FALSE)
  sigma_r <- raster::aggregate(sigma_r, fact = 2, fun = mean)

  if (display) raster::plot(sigma_r, col = viridis::viridis(25), main = "Conductivity roughness")
  verbose("   - Roughness conductivity map\n")

  # Edge-based conductivity
  e    <- param$conductivity$e
  sobl <- sobel.RasterLayer(slope)
  #plot(sigma_e, col = gray(1:30/30))
  sigma_e <- activation(sobl, e, "piecewise-linear", asc = FALSE)
  sigma_e <- raster::aggregate(sigma_e, fact = 2, fun = mean)

  if (display) raster::plot(sigma_e, col = viridis::viridis(25), main = "Conductivity Sobel edges")
  verbose("   - Sobel conductivity map\n")

  # Intensity-based conductivity
  sigma_i   <- dtm
  sigma_i[] <- 0
  if (use_intensity)
  {
    template2m <- raster::aggregate(dtm, fact = 2)
    template2m[] <- 0
    irange = intensity_range_by_flightline(las, template2m)
    #irange = terra::focal(irange, matrix(1,3,3), mean, na.rm = T)
    if (!is.null(mask)) irange <- raster::mask(irange, mask, inverse = T)

    if (display) raster::plot(irange, col = heat.colors(20), main = "Intensity range")

    #th <- stats::quantile(irange[], probs = q, na.rm = TRUE)
    th <- c(0.25,0.35)
    #sigma_i <- template2m
    sigma_i <- activation(irange, th, "piecewise-linear", asc = FALSE)
    #sigma_i = activation2(irange)

    if (display) raster::plot(sigma_i, col = viridis::viridis(20), main = "Conductivity intensity")
    verbose("   - Intensity conductivity map\n")
  }


  # CHM-based conductivity
  h <- param$conductivity$h
  chm <- lidR::grid_canopy(nlas, dtm, lidR::p2r())
  if (display) raster::plot(chm, col = height.colors(25),  main = "CHM")
  sigma_h <- dtm
  sigma_h <- activation(chm, h, "piecewise-linear", asc = FALSE)
  sigma_h <- raster::aggregate(sigma_h, fact = 2, fun = mean)

  if (display) raster::plot(sigma_h, col = viridis::inferno(25), main = "Conductivity CHM")
  verbose("   - CHM conductivity map\n")

  # Lowpoints-based conductivity
  # Check the presence/absence of lowpoints
  z1 <- 0.5
  z2 <- 3
  th <- 0.01
  lp <- density_lp_by_flightline(nlas, dtm)
  lp[is.na(lp)] = 0
  sigma_lp <- activation(lp, th, "thresholds", asc = FALSE)
  sigma_lp <- raster::aggregate(sigma_lp, fact = 2, fun = min)

  if (display) raster::plot(lp, col = viridis::inferno(20), main = "Number of low point")
  if (display) raster::plot(sigma_lp, col = viridis::inferno(2), main = "Bottom layer")
  verbose("   - Bottom layer conductivity map\n")

  # Density-based conductivity
  # Notice that the paper makes no mention of smoothing
  q <- param$conductivity$d

  d <- density_gnd_by_flightline(las, sigma_lp, drop_angles = 0)
  d[is.na(d)] = 0
  M   <- matrix(1,3,3)
  d = terra::focal(d, M, mean, na.rm = T)
  th <- mean(d[d>0], na.rm = T)
  d[is.na(d)] = 0

  if (!is.null(mask)) d <- raster::mask(d, mask, inverse = T, updatevalue = 0)
  if (display) raster::plot(d, col = viridis::inferno(15), main = "Density of ground points")

  sigma_d <- activation(d, c(0.33, 0.66), "piecewise-linear")

  if (display)  raster::plot(sigma_d, col = viridis::inferno(25), main = "Conductivity density")
  verbose("   - Density conductivity map\n")

  hard_slope = slope < 25
  hard_slope <- raster::aggregate(hard_slope, fact = 2, fun = min)


  # Final conductivity sigma
  alpha = param$conductivity$alpha
  alpha$i = alpha$i * as.numeric(use_intensity)
  max_coductivity <- sum(unlist(alpha))+1
  sigma <- (sigma_e + alpha$d * sigma_d + alpha$h * sigma_h + alpha$r * sigma_r + alpha$i * sigma_i)
  sigma <- sigma/max_coductivity

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Conductivity")

  sigma <- edge_enhancement(sigma, interation = 50, lambda = 0.05, k = 20, pass = 2)

  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[sigma < 0.1] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")
  sigma[hard_slope == 0] = 0.1
  #raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity")

  sigma[is.na(sigma)] <- 0.1 # lakes
  cells = raster::cellFromXY(sigma, bridge)
  sigma[cells] = 0.75
  if (display) raster::plot(sigma, col = viridis::inferno(25), main = "Edge enhanced conductivity with bridge")

  if (pkg == "terra") sigma <- terra::rast(sigma)

  dots = list(...)

  if (!is.null(dots$rawlayer))
  {
   dtm2 = raster::aggregate(dtm, fact = 2, fun = mean)
   slope2 = raster::aggregate(slope, fact = 2, fun = mean)
   roughness2 =  raster::aggregate(roughness, fact = 2, fun = mean)
   chm2 =  raster::aggregate(chm, fact = 2, fun = mean)
   lp2 =  raster::aggregate(lp, fact = 2, fun = mean)
   sobl2 =  raster::aggregate(sobl, fact = 2, fun = mean)
   u = raster::stack(slope2, roughness2, chm2, sigma_i, lp2, d, sobl2)
   names(u) = c("slope", 'rought', "chm", "intensity", "low", "density", "sobel")
   return(u)
  }

  return(sigma)
}

#' @export
rasterize_conductivity2.LAScluster = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  x <- lidR::readLAS(las)
  if (lidR::is.empty(x)) return(NULL)

  sigma <- rasterize_conductivity2(x, dtm, water, param, ...)
  sigma <- lidR:::raster_crop(sigma, lidR::st_bbox(las))
  return(sigma)
}

#' @export
rasterize_conductivity2.LAScatalog = function(las, dtm = NULL, water = NULL, param = alsroads_default_parameters2, ...)
{
  # Enforce some options
  if (lidR::opt_select(las) == "*") lidR::opt_select(las) <- "xyzciap"

  # Compute the alignment options including the case when res is a raster/stars/terra
  alignment <- lidR:::raster_alignment(1)

  if (lidR::opt_chunk_size(las) > 0 && lidR::opt_chunk_size(las) < 2*alignment$res)
    stop("The chunk size is too small. Process aborted.", call. = FALSE)

  # Processing
  options <- list(need_buffer = TRUE, drop_null = TRUE, raster_alignment = alignment, automerge = TRUE)
  output  <- lidR::catalog_apply(las, rasterize_conductivity2, dtm = dtm, water = water, param = param, ..., .options = options)
  return(output)
}

anisotropic_diffusion_filter = function(x, interation = 50, lambda = 0.2, k = 10)
{
  M = raster::as.matrix(x)
  M2 = anisotropic_diffusion(M*255, interation, lambda, k)
  M2 = M2/255
  y = x
  y[] = M2
  y
}

edge_enhancement = function(x, interation = 50, lambda = 0.05, k = 20, pass = 2)
{
  #smx = quantile(x[], na.rm = T, probs = 0.99)
  for (i in 1:pass)
    x = anisotropic_diffusion_filter(x, interation, lambda, k)

  x = raster::stretch(x, minv = 0, maxv = 1, maxq = 0.995)
  x
}

intensity_range_by_flightline <- function(las, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)

  if (length(ids) == 1)
  {
    ans = rasterize_intensityrange(las, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(las, PointSourceID == i)
    q = quantile(psi$Intensity, probs = 0.95)
    psi@data[Intensity>q, Intensity := q]
    ans[[as.character(i)]] <- rasterize_intensityrange(psi, res)
  }

  ans =terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

density_gnd_by_flightline <- function(las, res, drop_angles = 3)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)
  gnd = filter_ground(las)
  angle_max = max(abs(las$ScanAngleRank))
  gnd = filter_poi(gnd, abs(ScanAngleRank) < angle_max-drop_angles)

  if (length(ids) == 1)
  {
    ans = rasterize_density(gnd, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(gnd, PointSourceID == i)
    d <- rasterize_density(psi, res)
    d[d == 0] = NA
    q = quantile(d[], probs = 0.95, na.rm = TRUE)
    d[d>q] = q
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

density_lp_by_flightline <- function(las, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(las$PointSourceID)
  z1  <- 0.5
  z2  <- 3
  tmp <- lidR::filter_poi(las, Z > z1, Z < z2)

  if (length(ids) == 1)
  {
    ans = rasterize_density(gnd, res)
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids)
  {
    psi <- lidR::filter_poi(tmp, PointSourceID == i)
    d <- rasterize_density(psi, res)
    d[d == 0] = NA
    q = quantile(d[], probs = 0.95, na.rm = TRUE)
    d[d>q] = q
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans = terra::stretch(ans, maxv = 1)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = max, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

low_points_by_flightline <- function(nlas, res, th = 2)
{
  .N <- PointSourceID <- NULL

  z1 <- 0.5
  z2 <- 3
  res <- terra::rast(res)
  ids <- unique(nlas$PointSourceID)

  if (length(ids) == 1)
  {
    tmp <- lidR::filter_poi(nlas, Z > z1, Z < z2)
    ans = lidR:::rasterize_fast(tmp, res, 0, "count", pkg = "terra")
    return(raster::raster(ans > th))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(nlas, PointSourceID == i, Z > z1, Z < z2)
    d <- lidR:::rasterize_fast(psi, res, 0, "count", pkg = "terra")
    d[is.na(d)] = 0
    ans[[as.character(i)]] = d <= 3
  }

  ans = terra::rast(ans)
  #plot(ans, col = heat.colors(50))
  ans <- terra::app(ans, fun = min, na.rm = TRUE)
  #plot(ans, col = heat.colors(50))
  return(raster::raster(ans))
}

chm_by_flightline <- function(nlas, res)
{
  .N <- PointSourceID <- NULL

  res <- terra::rast(res)
  ids <- unique(nlas$PointSourceID)

  if (length(ids) == 1)
  {
    ans = rasterize_canopy(nlas, res, p2r())
    return(raster::raster(ans))
  }

  ans <- vector("list", 0)
  for (i in ids) {
    psi <- lidR::filter_poi(nlas, PointSourceID == i)
    d <- rasterize_canopy(psi, res, lidR::p2r())
    ans[[as.character(i)]] = d
  }

  ans = terra::rast(ans)
  ans = terra::mask(ans, terra::vect(water), inverse = T)
  #plot(ans, col = height.colors(50))
  ans = terra::stretch(ans, minv = 0)
  #plot(ans, col = height.colors(50))
  ans <- terra::app(ans, fun = min, na.rm = TRUE)
  #plot(ans, col = height.colors(50))
  ans =  ans - min(ans[], na.rm = TRUE)
  #plot(ans, col = height.colors(50))

  return(raster::raster(ans))
}


rasterize_intensityrange <- function(las, res)
{
  # Detect outliers of intensity and change their value. This is not perfect but avoid many troubles
  Intensity <- NULL
  outliers <- as.integer(stats::quantile(las$Intensity, probs = 0.98))
  las@data[Intensity > outliers, Intensity := outliers]

  # Switch Z and Intensity trick to use fast lidR internal function
  Z <- las[["Z"]]
  las@data[["Z"]] <-  las@data[["Intensity"]]
  imax <- lidR:::rasterize_fast(las, res, 0, "max")
  imin <- lidR:::rasterize_fast(las, res, 0, "min")
  irange <- imax - imin
  return(irange)
}


