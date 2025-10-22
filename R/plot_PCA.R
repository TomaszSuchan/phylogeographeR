#' Plot PCA Results
#'
#' Creates a scatter plot of principal component analysis results with variance explained.
#'
#' @param individuals Character vector of individual/sample names matching rows in eigenvecs
#' @param eigenvecs Matrix or data frame of eigenvectors (samples x PCs), where rows correspond to individuals
#' @param eigenvals Numeric vector of raw eigenvalues
#' @param popdata Data frame containing population/sample metadata. First column should contain 
#'   individual names matching the individuals parameter. No header required.
#' @param color_by Integer specifying column index in popdata to color points by. Default is 2.
#'   Note: Column index refers to the popdata columns (1 = individuals, 2 = first metadata column, etc.)
#' @param pc1 Integer specifying which PC to plot on x-axis. Default is 1
#' @param pc2 Integer specifying which PC to plot on y-axis. Default is 2
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Example data
#' individuals <- c("Rgl_BAN-01", "Rgl_BAN-02", "Rgl_BAR-01", "Rgl_BAR-02")
#' eigenvecs <- matrix(rnorm(20), ncol = 5)  # 4 samples, 5 PCs
#' eigenvals <- c(50, 30, 10, 7, 3)  # raw eigenvalues
#' 
#' # Population data (no header)
#' popdata <- data.frame(
#'   ind = individuals,
#'   pop = c("Rgl_BAN", "Rgl_BAN", "Rgl_BAR", "Rgl_BAR"),
#'   lat = c(49.20, 49.20, 49.20, 49.20),
#'   lon = c(19.72, 19.72, 20.20, 20.20),
#'   type = c("M", "M", "G", "G")
#' )
#' 
#' # Color by population (column 2)
#' plot_pca(individuals, eigenvecs, eigenvals, popdata, color_by = 2)
#' 
#' # Color by type (column 5)
#' plot_pca(individuals, eigenvecs, eigenvals, popdata, color_by = 5)
#' 
#' # Color by PC2
#' # First need to get PC column index: ncol(popdata) + 2
#' plot_pca(individuals, eigenvecs, eigenvals, popdata, 
#'          color_by = ncol(popdata) + 2)
#' }
#'
#' @export
plot_pca <- function(individuals, eigenvecs, eigenvals, popdata, color_by = 2, pc1 = 1, pc2 = 2) {
  
  # Convert eigenvectors to dataframe
  pca_df <- as.data.frame(eigenvecs)
  
  # Rename columns to match PC numbers
  colnames(pca_df) <- paste0("PC", 1:ncol(pca_df))
  
  # Add individuals as first column
  pca_df <- cbind(individual = individuals, pca_df)
  
  # Ensure popdata has column names
  if(is.null(colnames(popdata)) || colnames(popdata)[1] == "V1") {
    # Assign generic column names if none exist
    colnames(popdata) <- paste0("V", 1:ncol(popdata))
  }
  
  # Merge with population data based on individual names
  # Assuming first column of popdata contains individual names
  merged_df <- merge(pca_df, popdata, by.x = "individual", by.y = colnames(popdata)[1], all.x = TRUE)
  
  # Validate color_by index
  if(color_by < 1 || color_by > ncol(popdata)) {
    stop(sprintf("color_by must be between 1 and %d (total columns in popdata)", ncol(popdata)))
  }
  
  # Calculate variance explained from raw eigenvalues
  total_var <- sum(eigenvals)
  var_explained <- eigenvals / total_var * 100
  
  # Create axis labels with variance explained
  xlab <- sprintf("PC%d (%.1f%%)", pc1, var_explained[pc1])
  ylab <- sprintf("PC%d (%.1f%%)", pc2, var_explained[pc2])
  
  # Get column name for legend from popdata
  fill_label <- colnames(popdata)[color_by]
  
  # Get the color column from merged data
  # It's at position: 1 (individual) + ncol(pca_df)-1 (PCs) + color_by
  color_col_index <- ncol(pca_df) + color_by - 1
  
  # Create base plot
  p <- ggplot2::ggplot(merged_df, ggplot2::aes(x = .data[[paste0("PC", pc1)]], 
                                                y = .data[[paste0("PC", pc2)]],
                                                fill = .data[[color_col_index]])) +
    ggplot2::geom_point(size = 3, alpha = 0.7, pch = 21, color = "black") +
    ggplot2::labs(x = xlab,
                  y = ylab,
                  fill = fill_label)
  
  # Add appropriate scale based on whether color_by column is numeric
  if(is.numeric(merged_df[[color_col_index]])) {
    p <- p + ggplot2::scale_fill_gradient(low = "midnightblue", high = "lightblue")
  } else {
    p <- p + ggplot2::scale_fill_brewer(palette = "Set1")
  }
  
  return(p)
}