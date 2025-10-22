#' Plot PCA Results
#'
#' Creates a scatter plot of principal component analysis results with variance explained.
#' Uses combination of colors and shapes for categorical variables. Colors are reused 
#' with different shapes to create unique combinations.
#'
#' @param individuals Character vector of individual/sample names matching rows in eigenvecs and popdata
#' @param eigenvecs Matrix or data frame of eigenvectors (samples x PCs), where rows correspond to individuals
#' @param eigenvals Numeric vector of raw eigenvalues
#' @param popdata Data frame containing population/sample metadata.
#'   First column should contain individual names matching the individuals parameter.
#' @param color_by Integer specifying column index in popdata to color points by. Default is 2.
#' @param pc1 Integer specifying which PC to plot on x-axis. Default is 1
#' @param pc2 Integer specifying which PC to plot on y-axis. Default is 2
#' @param colors Character vector of colors to use for categorical variables. 
#'   Default is Set1 palette (9 colors). Colors will be recycled with different shapes.
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Read popdata
#' popdata <- read.table("popdata.txt", header = FALSE)
#' individuals <- popdata$V1
#' 
#' # Plot with default colors
#' plot_pca(individuals, eigenvecs, eigenvals, popdata, color_by = 2)
#' 
#' # Plot with custom colors
#' my_colors <- c("red", "blue", "green", "yellow", "purple", "orange")
#' plot_pca(individuals, eigenvecs, eigenvals, popdata, color_by = 2, colors = my_colors)
#' }
#'
#' @export
plot_pca <- function(individuals, eigenvecs, eigenvals, popdata, color_by = 2, pc1 = 1, pc2 = 2,
                     colors = RColorBrewer::brewer.pal(9, "Set1")) {
  
  # Convert eigenvectors to dataframe
  pca_df <- as.data.frame(eigenvecs)
  
  # Rename columns to match PC numbers
  colnames(pca_df) <- paste0("PC", 1:ncol(pca_df))
  
  # Add individuals as rownames (for reference)
  rownames(pca_df) <- individuals
  
  # Validate color_by index
  if(color_by < 1 || color_by > ncol(popdata)) {
    stop(sprintf("color_by must be between 1 and %d (number of columns in popdata)", ncol(popdata)))
  }
  
  # Combine PCA data with population data
  combined_df <- merge(pca_df, popdata, by.x = "row.names", by.y = "Individual", all.x = TRUE)
  
  # Calculate variance explained from raw eigenvalues
  total_var <- sum(eigenvals)
  var_explained <- eigenvals / total_var * 100
  
  # Create axis labels with variance explained
  xlab <- sprintf("PC%d (%.1f%%)", pc1, var_explained[pc1])
  ylab <- sprintf("PC%d (%.1f%%)", pc2, var_explained[pc2])
  
  # Get column name for legend
  fill_col_name <- colnames(popdata)[color_by]
  
  # Check if numeric or categorical
  if(is.numeric(combined_df[[fill_col_name]])) {
    # Numeric: use gradient
    p <- ggplot2::ggplot(combined_df, ggplot2::aes(x = .data[[paste0("PC", pc1)]], 
                                                    y = .data[[paste0("PC", pc2)]],
                                                    fill = .data[[fill_col_name]])) +
      ggplot2::geom_point(size = 3, alpha = 0.7, pch = 21, color = "black") +
      ggplot2::labs(x = xlab, y = ylab, fill = fill_col_name) +
      ggplot2::scale_fill_gradient(low = "midnightblue", high = "lightblue")
  } else {
    # Categorical: cycle colors with different shapes
    unique_groups <- unique(combined_df[[fill_col_name]])
    n_groups <- length(unique_groups)
    
    # Number of colors provided
    n_colors <- length(colors)
    
    # Available shapes (5 filled shapes: circle, square, diamond, triangle up, triangle down)
    available_shapes <- 21:25
    n_shapes <- length(available_shapes)
    
    # Assign colors and shapes
    # Color cycles every n_colors groups, shape changes every n_colors groups
    color_values <- rep(colors, length.out = n_groups)
    names(color_values) <- unique_groups
    
    shape_values <- available_shapes[((seq_along(unique_groups) - 1) %/% n_colors) %% n_shapes + 1]
    names(shape_values) <- unique_groups
    
    p <- ggplot2::ggplot(combined_df, ggplot2::aes(x = .data[[paste0("PC", pc1)]], 
                                                    y = .data[[paste0("PC", pc2)]],
                                                    fill = .data[[fill_col_name]],
                                                    shape = .data[[fill_col_name]])) +
      ggplot2::geom_point(size = 3, alpha = 0.7, color = "black") +
      ggplot2::scale_shape_manual(values = shape_values) +
      ggplot2::scale_fill_manual(values = color_values) +
      ggplot2::labs(x = xlab, y = ylab, fill = fill_col_name, shape = fill_col_name)
  }
  
  return(p)
}
