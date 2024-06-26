#' PCoA plot function
#' 
#' Plot PCoA from AFLP data (individuals with missing data are removed!), ggplot-based so can be styled with additional ggplot arguments, eg. AFLP_plot_PCoA + theme(axis.text = element_text(size = 10))
#' @import ggplot2
#' @import vegan
#' @param structure_input_path AFLP matrix where rows reperesent individuals and columns the markers, no header, first column are individual names in the format Population-Individual (eg. POP1-1), if any additional column should be ommited (eg. with population number) it should be passed using the 'remove' argument. Assumes that the first row should be removed and duplicated rows as well (that's how the structure input is usually coded for AFLP).
#' @param population_data_path path for the population data file (optional) where the first column is the population, second y-coordinate and third x-coordinate. The rest of the columns can have strata information (defined by the strata argument). Can have a header as long it does not contain any of the population names!
#' @param strata coulumn in the population data file used for data stratification
#' @param axes vector containing the PCoA axes to be plotted
#' @param remove vector of colums to be removed, defalults to usual population column in structure files. Set as remove=NULL if no column should be removed.
#' @export
#' @examples
#' structure_input_path <- "test_data/structure/Aal_carp_structure-inputD.txt"
#' population_data_path <- "test_data/structure/popdata.txt"
#' AFLP_plot_PCoA(structure_input_path, population_data_path, strata=5, remove=NULL, axes=c(2,3))

AFLP_plot_PCoA <- function(structure_input_path, population_data_path=NULL, strata = 1, axes=c(1,2), remove=2) {
  
  aflp_matrix <- load_AFLP(structure_input_path, remove)
  
  # Extract population information from individual names
  colnames(aflp_matrix)[1] <- "Individual"
  split_names <- strsplit(as.character(aflp_matrix$Individual), "-")
  Population <- sapply(split_names, function(x) x[1])
  aflp_matrix <- cbind(Population, aflp_matrix)

  # PCoA
  PCoA <- capscale(aflp_matrix[, -c(1, 2)] ~ 1, distance = "jaccard", na.action = na.omit)

  # Color population groupings if population_data_path is provided:
  if(!is.null(population_data_path)){
  population_data <- read.table(population_data_path)
  populations <- merge(as.data.frame(aflp_matrix$Population),
                       population_data,
                       by.x="aflp_matrix$Population",
                       by.y="V1")
  scores <- cbind(populations[,strata], as.data.frame(scores(PCoA, display = "sites", choices=axes)))
  colnames(scores)[1] <- "Group"
  nmds_plot_new <- ggplot(scores, aes_string(x = names(scores)[2], y = names(scores)[3], fill = names(scores)[1])) +
    coord_fixed() +
    geom_point(size = 2, pch = 21, stroke = 0.3, alpha = 0.8, colour = "black")
  } else {
  nmds_plot_new <- ggplot(as.data.frame(scores(PCoA, display = "sites")), aes(x = MDS1, y = MDS2)) +
    coord_fixed() +
    geom_point(size = 2, pch = 21, stroke = 0.3, alpha = 0.5, colour = "black", fill = "black")
  }
  return(nmds_plot_new)
}
