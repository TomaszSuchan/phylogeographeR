#' Prepare EEMS input
#' 
#' Prepare input files for eems. Distances are calculated using jaccard distance, individuals with the same distances are removed (closely related or clonal individuals) 
#' @import vegan
#' @param output basename of the output files
#' @param structure_input_path AFLP matrix where rows reperesent individuals and columns the markers, no header, first column are individual names in the format Population-Individual (eg. POP1-1), if any additional column should be ommited (eg. with population number) it should be passed using the 'remove' argument. Assumes that the first row should be removed and duplicated rows as well (that's how the structure input is usually coded for AFLP).
#' @param population_data_path path for the population data file where the first column is the population, second y-coordinate and third x-coordinate. The rest of the columns are ignored. Can have a header as long it does not contain any of the population names!
#' @param remove vector of colums to be removed, defalults to usual population column in structure files. Set as remove=NULL if no column should be removed.
#' @export
#' @examples
#'
#' output <- "Aal" 
#' structure_input_path <- "test_data/structure/Aal_carp_structure-inputD.txt"
#' population_data_path <- "test_data/structure/popdata.txt"
#' AFLP_prepare_eems(output, structure_input_path, population_data_path, remove=2)

AFLP_prepare_eems <- function(output, structure_input_path, population_data_path, remove=2) {
  
  aflp_matrix <- load_AFLP(structure_input_path, remove)

  # Remove clones
  print("Clonal individuals:")
  print(aflp_matrix[duplicated(aflp_matrix[, -1]), 1])
  aflp_matrix <- aflp_matrix[!duplicated(aflp_matrix[, -1]), ]

  # Extract population information from individual names
  colnames(aflp_matrix)[1] <- "Individual"
  split_names <- strsplit(as.character(aflp_matrix$Individual), "-")
  Population <- sapply(split_names, function(x) x[1])
  aflp_matrix <- cbind(Population, aflp_matrix)

  # distances
  dist_matrix <- as.matrix(vegdist(aflp_matrix[, -c(1, 2)], method="jaccard", na.rm=TRUE))
  duplicates <- duplicated(dist_matrix)
  print("Individuals with the same distances (removed):")
  print(aflp_matrix[duplicates,1]) 
  write.table(dist_matrix[!duplicates, !duplicates],
              file=paste(output, ".diffs", sep=""),
              quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  
  # popdata
  population_data <- read.table(population_data_path)
  popdata <- merge(aflp_matrix[,c(1,2)], population_data, by.x="Population", by.y="V1")

  write.table(popdata[!duplicates, c(4,3)],
              file=paste(output, ".coord", sep=""),
              quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
  
  # stats

  print("Number of individuals:")
  print(nrow(dist_matrix[!duplicated(dist_matrix), ]))
  print("Number of markers:")
  print(ncol(aflp_matrix[, -c(1, 2)]))
}