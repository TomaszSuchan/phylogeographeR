#' PCoA plot function
#' 
#' Plot PCoA from AFLP data (individuals with missing data are removed!), ggplot-based so can be styled with additional ggplot arguments, eg. AFLP_plot_PCoA + theme(axis.text = element_text(size = 10))
#' @import ggplot2
#' @import vegan
#' @param aflp_matrix_path AFLP matrix where rows reperesent individuals and columns the markers, no header, first column are individual names in the format Population-Individual (eg. POP1-1), if any additional column should be ommited (eg. with population number) it should be passed using the 'remove' argument. Assumes that the first row should be removed and duplicated rows as well (that's how the structure input is usually coded for AFLP).
#' @param population_data_path path for the population data file where the first column is the population and the second the strata. Can have a header as long it does not contain any of the population names!
#' @param strata coulumn in the population data file used for data stratification
#' @param choices vector containing the PCoA axes to be plotted
#' @param remove vector of colums to be removed, defalults to usual population column in structure files. Set as remove=NULL if no column should be removed.
#' @export
#' @examples
#' aflp_matrix_path <- "struct-new/original_input/Cal_Cal_Structure.txt"
#' population_data_path <- "population_data.txt"
#' AFLP_plot_PCoA(aflp_matrix_path, population_data_path, strata=4, remove=NULL, choices=c(2,3))

AFLP_plot_PCoA <- function(aflp_matrix_path, population_data_path=NULL, strata = 1, choices=c(1,2), remove=2) {
   
  aflp_matrix <- read.table(aflp_matrix_path, head=FALSE, sep="\t", row.names=NULL, skip=1)
  
  # hack to remove last column full of NA
  aflp_matrix <- aflp_matrix[,colSums(is.na(aflp_matrix))<nrow(aflp_matrix)]

  # Remove duplicated rows (as in Structure input files)
  aflp_matrix <- aflp_matrix[!duplicated(aflp_matrix), ]

  colnames(aflp_matrix)[1] <- "Individual"

  # Extract population information from individual names
  split_names <- strsplit(as.character(aflp_matrix$Individual), "-")
  Population <- sapply(split_names, function(x) x[1])
  aflp_matrix <- cbind(Population, aflp_matrix)



  # Remove rows (samples) containing NA in AFLP data
  aflp_matrix[aflp_matrix==-9] <- NA
  aflp_matrix <- na.omit(aflp_matrix)

  # Remove columns containing the same values
  cols1 <- colSums(aflp_matrix[, -c(1, 2, remove)])[colSums(aflp_matrix[, -c(1, 2, remove)]) == dim(aflp_matrix)[1]]
  cols1 <- names(cols1)
  cols0 <- colSums(aflp_matrix[, -c(1, 2, remove)])[colSums(aflp_matrix[, -c(1, 2, remove)]) == 0]
  cols0 <- names(cols0)
  print("Removing columns with monomorphic markers:")
  print(c(cols1, cols0))
  if(length(cols1)>0){
    aflp_matrix <- aflp_matrix[,!(names(aflp_matrix) %in% cols1)]
  }
  if(length(cols0)>0){
    aflp_matrix <- aflp_matrix[,!(names(aflp_matrix) %in% cols0)]
  }

  # PCoA
  PCoA <- capscale(aflp_matrix[, -c(1, 2, remove)] ~ 1, distance = "jaccard", na.action = na.omit)

  # Color population groupings if population_data_path is provided:
  if(!is.null(population_data_path)){
  population_data <- read.table(population_data_path)
  populations <- merge(as.data.frame(aflp_matrix$Population),
                       population_data,
                       by.x="aflp_matrix$Population",
                       by.y="V1")
  scores <- cbind(populations[,strata], as.data.frame(scores(PCoA, display = "sites", choices=choices)))
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
