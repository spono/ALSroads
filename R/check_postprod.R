#' Check amplitude of differences between corrected and uncorrected roads
#'
#' Check amplitude of differences between corrected and uncorrected roads
#'
#' @param roads  multiple lines (\code{sf} format). Corrected roads.
#' @param ref  multiple lines (\code{sf} format). Original non-corrected roads.
#' @param field  character. Unique identifier field in both road datasets.
#'
#' @return data.frame with metrics about each roads.
#' @export
#' @examples
#' library(sf)
#'
#' road <- system.file("extdata", "j5gr_centerline_971487.gpkg", package="ALSroads")
#' road_cor <- st_read(road, "corrected", quiet = TRUE)
#' road_ori <- st_read(road, "original", quiet = TRUE)
#'
#' df <- check_road_differences(road_cor, road_ori, "objectid")
#' df
check_road_differences <- function(roads, ref, field)
{
  IDs1 <- sort(unique(roads[[field]]))
  IDs2 <- sort(unique(ref[[field]]))
  if (length(IDs1) != length(roads[[field]])) stop("Values in unique identifier field are not unique for 'roads'.", call. = FALSE)
  if (length(IDs2) != length(ref[[field]])) stop("Values in unique identifier field are not unique for 'ref'.", call. = FALSE)
  if (!all(IDs1 == IDs2)) stop("Values in unique identifier field are not the same in both road datasets.", call. = FALSE)

  # Arrange both road datasets to make sure
  # that line indices will match between them
  roads <- dplyr::arrange(roads, field)
  ref <- dplyr::arrange(ref, field)

  # Compute metrics
  ratios_area_perimeter <- mapply(diff_area_perimeter, sf::st_geometry(roads), sf::st_geometry(ref))
  quantiles_along_road <- mapply(diff_along_road, sf::st_geometry(roads), sf::st_geometry(ref), SIMPLIFY = FALSE) |> do.call(what = rbind)

  # Format results
  df_results <- data.frame(
    field = ref[[field]],
    area_over_perimeter = ratios_area_perimeter) |>
    cbind(quantiles_along_road)

  names(df_results)[1] <- field

  return(df_results)
}


#' Compute a difference index between two roads by considering them as a polygon
#'
#' Connect the two roads by their ends to construct a polygon from which
#' a ratio of its area over its perimeter will be calculated.
#' The larger the value, the larger the differences between the two roads. It
#' must be noted that in some edge cases, a low value doesn't mean a low difference.
#'
#' @param road_cor  line (\code{sf} or \code{sfc} format). Corrected road.
#' @param road_ori  line (\code{sf} or \code{sfc} format). Original non-corrected road.
#' @param graph  boolean. Whether or not to display graphics.
#'
#' @return numeric. Ratio of the area over the perimeter of the constructed polygon.
#' @noRd
#' @examples
#' library(sf)
#'
#' road <- system.file("extdata", "j5gr_centerline_971487.gpkg", package="ALSroads")
#' road_cor <- st_read(road, "corrected", quiet = TRUE)
#' road_ori <- st_read(road, "original", quiet = TRUE)
#'
#' ratio <- ALSroads:::diff_area_perimeter(road_cor, road_ori, graph = TRUE)
diff_area_perimeter <- function(road_cor, road_ori, graph = FALSE)
{
  road_ori <- sf::st_geometry(road_ori)
  road_cor <- sf::st_geometry(road_cor)

  if (length(road_ori) > 1 | length(road_cor) > 1) stop("'road_ori' and 'road_cor' must contain only one feature.", call. = FALSE)

  # Extract vertices
  coords_ori <- sf::st_coordinates(road_ori)[,-3]
  coords_cor <- sf::st_coordinates(road_cor)[,-3]


  # Adjust coordinates in order to create a valid sequence
  # of vertices for a polygon
  dist_to_start <- stats::dist(rbind(coords_ori[1,], coords_cor[1,]))[1]
  dist_to_end <- stats::dist(rbind(coords_ori[1,], coords_cor[nrow(coords_cor),]))[1]
  closest_vertex <- which.min(c(dist_to_start, dist_to_end))

  if (closest_vertex == 1) coords_cor <- apply(coords_cor, 2, rev)


  # Construct polygon
  poly <- rbind(coords_ori, coords_cor, coords_ori[1,]) |>
    list() |>
    sf::st_polygon() |>
    sf::st_sfc(crs = sf::st_crs(road_ori))


  # Compute area/perimeter ratio of polygon
  ratio <- as.numeric(sf::st_area(poly) / lwgeom::st_perimeter(poly))


  # Make graphical representation of the differences
  if (graph)
  {
    coords <- rbind(coords_ori, coords_cor)
    limits <- list(x = range(coords[,1]), y = range(coords[,2]))
    bigger <- as.numeric(diff(limits$x) < diff(limits$y)) + 1
    offset <- mean(limits[[bigger]]) - limits[[bigger]][1]
    limits$x <- mean(limits$x) + c(-offset, offset)
    limits$y <- mean(limits$y) + c(-offset, offset)

    plot(poly, col = "orange", axes = TRUE, xlim=limits$x, ylim = limits$y)
    plot(road_ori, lwd = 2, col = "red", add = TRUE)
    plot(road_cor, lwd = 2, col = "darkgreen", add = TRUE)
    graphics::title(sprintf("Ratio area/perimeter: %.1f m\u00b2/m", ratio))
  }
  return(ratio)
}


#' Compute distances quantiles based on how far apart the two roads are
#'
#' Sample points at regular interval each road and measure distances between
#' each pair of points. The larger the quantiles are, the larger the differences
#' between the two roads.
#'
#' @param road_cor  line (\code{sf} or \code{sfc} format). Corrected road.
#' @param road_ori  line (\code{sf} or \code{sfc} format). Original non-corrected road.
#' @param step  numeric (distance unit). Interval on \code{road_ori} at which sample differences.
#' @param graph  boolean. Whether or not to display graphics.
#'
#' @return named numeric. Quantiles P50, P90, P100 along with the maximal difference found at both ends.
#' @noRd
#' @examples
#' library(sf)
#'
#' road <- system.file("extdata", "j5gr_centerline_971487.gpkg", package="ALSroads")
#' road_ori <- st_read(road, "original", quiet = TRUE)
#' road_cor <- st_read(road, "corrected", quiet = TRUE)
#'
#' plot(st_geometry(road_ori))
#' plot(st_geometry(road_cor), add = T)
#'
#' p <- ALSroads:::diff_along_road(road_cor, road_ori, graph = TRUE)
diff_along_road <- function(road_cor, road_ori, step = 10, graph = FALSE)
{
  road_ori <- sf::st_geometry(road_ori)
  road_cor <- sf::st_geometry(road_cor)

  if (length(road_ori) > 1 | length(road_cor) > 1) stop("'road_ori' and 'road_cor' must contain only one feature.", call. = FALSE)

  st_full_line_sample <- function(line, n_steps)
  {
    coords <- sf::st_coordinates(line)[,-3]
    pts_extrems <- c(1, nrow(coords)) |>
      lapply(function(x) { sf::st_point(coords[x,]) }) |>
      sf::st_sfc(crs = sf::st_crs(road_ori))

    pts_middle <- sf::st_cast(sf::st_line_sample(line, n_steps), "POINT")
    pts <- c(pts_extrems[1], pts_middle, pts_extrems[2])
  }


  # Sample points along the two roads
  n_steps <- ceiling(sf::st_length(road_ori) / step) |> as.numeric()
  ls_pts <- lapply(c(road_ori, road_cor), st_full_line_sample, n_steps)
  names(ls_pts) <- c("road_ori", "road_cor")

  # Calculate distance between each pair of points
  dist_step <- sf::st_distance(ls_pts[["road_ori"]], ls_pts[["road_cor"]]) |>
    diag() |>
    as.numeric()

  # Compute P50, P90 and P100
  # Large values, especially at P50, might indicate that the corrected road
  # took a completely new (and wrong) path
  p <- stats::quantile(dist_step, probs = c(0.5, 0.90, 1))
  names(p) <- c("P50", "P90", "P100")

  # Compute the maximal difference distance between the pairs at the first and last vertex.
  # A value near the one of alsroads_default_parameters$extraction$road_buffer might
  # indicate a problem right at the start/end of the least-cost path extraction process
  dist_end_max <- c(end_max = max(dist_step[c(1, length(dist_step))]))


  # Display graphical representation of the differences
  if (graph)
  {
    coords_ori <- sf::st_coordinates(ls_pts[["road_ori"]])
    coords_cor <- sf::st_coordinates(ls_pts[["road_cor"]])

    lines <- 1:nrow(coords_ori) |>
      lapply(function(x) { sf::st_linestring(rbind(coords_ori[x,], coords_cor[x,])) }) |>
      sf::st_sfc(crs = sf::st_crs(road_ori))

    df_dist <- data.frame(idx = 1:length(dist_step), distance = dist_step)

    coords <- rbind(coords_ori, coords_cor)
    limits <- list(x = range(coords[,1]), y = range(coords[,2]))
    bigger <- as.numeric(diff(limits$x) < diff(limits$y)) + 1
    offset <- mean(limits[[bigger]]) - limits[[bigger]][1]
    limits$x <- mean(limits$x) + c(-offset, offset)
    limits$y <- mean(limits$y) + c(-offset, offset)

    plot(road_ori, lwd = 2, col = "red", axes = TRUE, xlim=limits$x, ylim = limits$y)
    plot(road_cor, lwd = 2, col = "darkgreen", add = TRUE)
    plot(lines, col = "purple", add = TRUE)
    graphics::title("Pairwise points along roads")

    plot(df_dist, pch = 20)
    graphics::lines(stats::predict(stats::loess(distance~idx, df_dist)), col = "blue", lwd = 2)
    graphics::abline(h = p)
    graphics::title("Pairwise distances along roads", sprintf("P50: %.1f m | P90: %.1f m | P100: %.1f m", p[1], p[2], p[3]))
  }

  return(c(p, dist_end_max))
}
