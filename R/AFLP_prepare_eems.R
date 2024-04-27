#' Prepare EEMS input
#' 
#' Prepare input files for eems. Distances are calculated using jaccard distance, individuals with the same distances are removed (closely related or clonal individuals) 
#' @import vegan
#' @param aflp_matrix_path AFLP matrix where rows reperesent individuals and columns the markers, no header, first column are individual names in the format Population-Individual (eg. POP1-1), if any additional column should be ommited (eg. with population number) it should be passed using the 'remove' argument.
#' @param population_data_path path for the population data file where the first column is the population, second the strata, third y and fourth z coordinate. Can have a header as long it does not contain any of the population names!
#' @param species name of the output files
#' @param remove vector of colums to be removed
#' @export
#' @examples
#' 
#' aflp_matrix_path <- "struct-new/original_input/Cal_Cal_Structure.txt"
#' population_data_path <- "population_data.txt"
#' AFLP_plot_structure_map(aflp_matrix_path, population_data_path, c(2,3))AFLP_plot_PCoA(aflp_matrix_path=, population_data_path=NULL, choices=c(1,2))

AFLP_prepare_eems <- function(species, aflp_matrix_path, population_data_path, remove=2) {
  
  aflp_matrix <- read.table(aflp_matrix_path, head=T, sep="\t")
  colnames(aflp_matrix)[1] <- "Individual"

  # Rename the markers
  colnames(aflp_matrix) <- gsub(".", "-", colnames(aflp_matrix), fixed=T)

  # Extract population information from individual names
  split_names <- strsplit(as.character(aflp_matrix$Individual), "-")
  Population <- sapply(split_names, function(x) x[1])
  aflp_matrix <- cbind(Population, aflp_matrix)

  # Remove rows (samples) containing NA in AFLP data
  aflp_matrix[aflp_matrix==-9] <- NA
  aflp_matrix <- na.omit(aflp_matrix)

  # Remove columns containing the same values
  cols1 <- colSums(aflp_matrix[, -c(1, remove)])[colSums(aflp_matrix[, -c(1, remove)]) == dim(aflp_matrix)[1]]
  cols1 <- names(cols1)
  cols0 <- colSums(aflp_matrix[, -c(1, remove)])[colSums(aflp_matrix[, -c(1, remove)]) == 0]
  cols0 <- names(cols0)
  if(length(cols1)>0){
    aflp_matrix <- aflp_matrix[,!(names(aflp_matrix) %in% cols1)]
  }
  if(length(cols0)>0){
    aflp_matrix <- aflp_matrix[,!(names(aflp_matrix) %in% cols0)]
  }

  # Remove clones
  print("Clonal individuals:")
  print(aflp_matrix[duplicated(aflp_matrix), 1])
  aflp_matrix <- aflp_matrix[!duplicated(aflp_matrix), ]
  
  # distances
  dist_matrix <- as.matrix(vegdist(aflp_matrix[, -c(1, remove)], method="jaccard", na.rm=TRUE))
  duplicates <- duplicated(dist_matrix)
  print("Clonal individuals:")
  print(aflp_matrix[duplicates,1]) 
  write.table(dist_matrix[!duplicates, !duplicates],
              file=paste(species, ".diffs", sep=""),
              quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  
  # popdata
  population_data <- read.table(population_data_path, header=TRUE)
  popdata <- merge(aflp_matrix[,c(1,remove)], population_data, by.x="Population", by.y="IBD_GRID_CELL")

  write.table(popdata[!duplicates, c(5,4)],
              file=paste(species, ".coord", sep=""),
              quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  
  # stats

  print("Number of individuals:")
  print(nrow(dist_matrix[!duplicated(dist_matrix), ]))
  print("Number of markers:")
  print(ncol(aflp_matrix[, -c(1, remove)]))
}