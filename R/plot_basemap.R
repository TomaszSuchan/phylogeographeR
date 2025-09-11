#' Plot Base Map
#'
#' @description
#' Creates a base map for plotting geographic data. Accepts either direct coordinates
#' or a data frame with lon/lat columns.
#'
#' @param xmin numeric. Minimum longitude (western boundary) in WGS84
#' @param xmax numeric. Maximum longitude (eastern boundary) in WGS84
#' @param ymin numeric. Minimum latitude (southern boundary) in WGS84
#' @param ymax numeric. Maximum latitude (northern boundary) in WGS84
#' @param data data.frame. Alternative to xmin/xmax/ymin/ymax. Must contain 'lon' and 'lat' columns
#' @param expand numeric. Proportion to expand boundaries (e.g., 0.1 for 10%). Default is 0.
#' @param crs numeric. EPSG code for map projection (default: 4326 for WGS84)
#' @param basemap NULL, SpatRaster, or sf object. Custom basemap. If NULL, uses Natural Earth world boundaries
#' @param land_colour character. Hex code for land color (default: "#d9d9d9")
#' @param sea_colour character. Hex code for sea color (default: "#deebf7")
#' @param add_scale logical. Add scale bar (default: FALSE)
#' @param add_arrow logical. Add north arrow (default: FALSE)
#' @param title character. Plot title (default: "")
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import sf
#' @importFrom ggspatial annotation_scale annotation_north_arrow north_arrow_orienteering layer_spatial
#' @importFrom rnaturalearthdata countries50
#'
#' @examples
#' \dontrun{
#' # Direct coordinates
#' plot_basemap(xmin = -10, xmax = 10, ymin = 40, ymax = 60)
#' 
#' # From data
#' sites <- data.frame(lon = c(-5, 0, 5), lat = c(50, 52, 54))
#' plot_basemap(data = sites, expand = 0.1)
#' 
#' # With raster basemap
#' library(terra)
#' earth <- terra::rast("path/to/raster.tif")
#' plot_basemap(xmin = -10, xmax = 10, ymin = 40, ymax = 60, basemap = earth)
#' 
#' # With decorations
#' plot_basemap(xmin = -10, xmax = 10, ymin = 40, ymax = 60,
#'              add_scale = TRUE, add_arrow = TRUE, title = "Europe")
#' 
#' # Different projection
#' plot_basemap(xmin = -10, xmax = 10, ymin = 40, ymax = 60, crs = 3035)
#'
#' @export
plot_basemap <- function(
    xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
    data = NULL,
    expand = 0,
    crs = 4326,
    basemap = NULL,
    land_colour = "#d9d9d9",
    sea_colour = "#deebf7",
    add_scale = FALSE,
    add_arrow = FALSE,
    title = ""
) {
    # Get boundaries from either direct input or data
    if (!is.null(data)) {
        if (!all(c("lon", "lat") %in% names(data))) {
            stop("data must contain 'lon' and 'lat' columns")
        }
        xmin <- min(data$lon, na.rm = TRUE)
        xmax <- max(data$lon, na.rm = TRUE)
        ymin <- min(data$lat, na.rm = TRUE)
        ymax <- max(data$lat, na.rm = TRUE)
    }
    
    # Validate boundaries
    if (is.null(xmin) || is.null(xmax) || is.null(ymin) || is.null(ymax)) {
        stop("Must provide either xmin/xmax/ymin/ymax or data")
    }
    
    # Apply expansion
    if (expand > 0) {
        x_range <- xmax - xmin
        y_range <- ymax - ymin
        xmin <- xmin - (x_range * expand)
        xmax <- xmax + (x_range * expand)
        ymin <- ymin - (y_range * expand)
        ymax <- ymax + (y_range * expand)
    }
    
    # Setup CRS
    crs <- sf::st_crs(crs)
    
    # Check basemap type and process accordingly
# Check basemap type and process accordingly
is_raster <- FALSE
if (!is.null(basemap)) {
    basemap_class <- class(basemap)[1]
    is_raster <- basemap_class == "SpatRaster"
}

# Get basemap
if (is.null(basemap)) {
    # Use Natural Earth data
    world_data <- rnaturalearthdata::countries50
    map_data <- if ("geometry" %in% names(world_data)) world_data["geometry"] else world_data[, ncol(world_data)]
    map_data <- sf::st_transform(map_data, crs)
    
} else if (is_raster) {
    # Handle raster basemap safely
    if (!requireNamespace("terra", quietly = TRUE)) {
        stop("terra package is required for raster basemaps. Install with: install.packages('terra')")
    }
    if (!requireNamespace("ggspatial", quietly = TRUE)) {
        stop("ggspatial package is required for plotting rasters. Install with: install.packages('ggspatial')")
    }
    
    # Create target bbox polygon
    corners <- data.frame(
        lon = c(xmin, xmax, xmax, xmin, xmin),
        lat = c(ymin, ymin, ymax, ymax, ymin)
    )
    bbox_target <- sf::st_as_sf(corners, coords = c("lon", "lat"), crs = crs) |>
                   dplyr::summarise(geometry = sf::st_union(geometry))
    
    # Transform bbox back to raster CRS
    bbox_orig <- sf::st_transform(bbox_target, sf::st_crs(basemap))
    
    # Convert sf bbox to SpatExtent
    ext_orig <- terra::ext(terra::vect(bbox_orig))
    
    # Crop raster in original CRS
    basemap_cropped <- terra::crop(basemap, ext_orig)
    
    # Then project cropped raster
    basemap_projected <- if (crs$epsg != sf::st_crs(basemap)$epsg) {
        terra::project(basemap_cropped, paste0("EPSG:", crs$epsg))
    } else {
        basemap_cropped
    }
    
    map_data <- basemap_projected
} else {
    # Assume it's an sf object
    map_data <- sf::st_transform(basemap, crs)
}
    
    # Create the plot limits in target CRS
    # Transform corner points from WGS84 to target CRS
    corners <- data.frame(
        lon = c(xmin, xmax, xmax, xmin),
        lat = c(ymin, ymin, ymax, ymax)
    )
    corners_sf <- sf::st_as_sf(corners, coords = c("lon", "lat"), crs = 4326)
    corners_transformed <- sf::st_transform(corners_sf, crs)
    coords_transformed <- sf::st_coordinates(corners_transformed)
    
    xlim <- range(coords_transformed[, "X"])
    ylim <- range(coords_transformed[, "Y"])
    
    # Build plot
    plt <- ggplot2::ggplot()
    
    # Add basemap layer (raster or vector)
    if (is_raster) {
        # Add raster layer
        plt <- plt + ggspatial::layer_spatial(map_data)
    } else {
        # Add vector layer
        plt <- plt + ggplot2::geom_sf(
            data = map_data,
            fill = land_colour,
            colour = "black",
            linewidth = 0.1
        )
    }
    
    # Add coordinate system and styling
    plt <- plt +
        ggplot2::coord_sf(
            xlim = xlim,
            ylim = ylim,
            crs = crs,
            datum = crs,
            expand = FALSE
        ) +
        ggplot2::labs(title = title, x = "Longitude", y = "Latitude") +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = sea_colour),
            panel.border = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.3),
            panel.grid = ggplot2::element_line(colour = "white", linewidth = 0.1),
            axis.text = ggplot2::element_text(colour = "black", size = 8),
            axis.title = ggplot2::element_text(colour = "black", size = 10),
            plot.title = ggplot2::element_text(face = "bold", size = 12)
        )
    
    # Add scale bar
    if (add_scale) {
        plt <- plt + ggspatial::annotation_scale(
            location = "bl",
            bar_cols = c("black", "white"),
            line_width = 0.3,
            height = grid::unit(0.15, "cm"),
            width_hint = 0.15,
            text_cex = 0.6,
            unit_category = "metric",
            style = "ticks"
        )
    }
    
    # Add north arrow
    if (add_arrow) {
        plt <- plt + ggspatial::annotation_north_arrow(
            location = "tr",
            which_north = "true",
            height = grid::unit(0.5, "cm"),
            width = grid::unit(0.5, "cm"),
            pad_x = grid::unit(0.25, "cm"),
            pad_y = grid::unit(0.25, "cm"),
            style = ggspatial::north_arrow_orienteering(
                text_size = 6,
                line_width = 0.5
            )
        )
    }
    
    return(plt)
}
