#' NJ tree plot function
#' 
#' Plot Jaccard distance based NJ tree from AFLP data (individuals with missing data are removed!), ggplot-based so can be styled with additional ggplot arguments, eg. AFLP_plot_NJtree + theme(axis.text = element_text(size = 10))
#' @import ggplot2
#' @import ggtree
#' @import ape 
#' @import vegan
#' @param structure_input_path AFLP matrix where rows reperesent individuals and columns the markers, no header, first column are individual names in the format Population-Individual (eg. POP1-1), if any additional column should be ommited (eg. with population number) it should be passed using the 'remove' argument. Assumes that the first row should be removed and duplicated rows as well (that's how the structure input is usually coded for AFLP).
#' @param population_data_path path for the population data file (optional) where the first column is the population, second is  and the second the strata. Can have a header as long it does not contain any of the population names!
#' @param strata coulumn in the population data file used for data stratification
#' @param remove vector of colums to be removed, defalults to usual population column in structure files. Set as remove=NULL if no column should be removed.
#' @param labels should tip labels be printed or not
#' @export
#' @examples
#' structure_input_path <- "test_data/structure/Aal_carp_structure-inputD.txt"
#' population_data_path <- "test_data/structure/popdata.txt"
#' AFLP_plot_NJtree(structure_input_path, population_data_path, strata=5, remove=NULL)


AFLP_plot_NJtree <- function(structure_input_path, population_data_path=NULL, strata=4, remove = NULL, labels = TRUE) {
  aflp_matrix <- load_AFLP(structure_input_path, remove)
  rownames(aflp_matrix) <- aflp_matrix[,1]
  aflp_matrix <- as.matrix(aflp_matrix[, -1])

  nj_tree <- function(x) nj(vegdist(x, method="jaccard"))

  tr <- nj_tree(aflp_matrix)


  #tr$node.label <- boot.phylo(tr, aflp_matrix, nj_tree, B=1000)

  #add popdata

  # Extract population information from individual names
  
  split_names <- strsplit(as.character(tr$tip.label), "-")
  pop <- sapply(split_names, function(x) x[1])

  nodes_population_data <- as.data.frame(cbind(tr$tip.label, pop))
  colnames(nodes_population_data) <- c("ind", "pop")
  
  # Color population groupings if population_data_path is provided:
  if(!is.null(population_data_path)){
  population_data <- read.table(population_data_path)
  populations <- merge(nodes_population_data,
                       population_data[,c(1, strata)],
                       by.x="pop",
                       by.y="V1")


  groups <- split(populations[, 2], populations[, 3])
  
  tr <- groupOTU(tr, groups, "Group")

  tree_plot <- ggtree(tr, aes(color=group), layout="daylight")

  } else {
  tree_plot <- ggtree(tr, layout="daylight")
  }

  if(labels == TRUE){
    tree_plot <- tree_plot + geom_tiplab(color="black")
  }

  return(tree_plot)
}
