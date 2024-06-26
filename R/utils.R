#' Load AFLP matrix, remove column if specified, remove duplicated rows (as in Structure input files), remove samples containing NA, remove monomorphic markers
#' @noRd

load_AFLP <- function(structure_input_path, remove_col=NULL) {

  aflp_matrix <- read.table(structure_input_path, head=FALSE, sep="\t", row.names=NULL, skip=1)
  
  # hack to remove last column full of NA
  aflp_matrix <- aflp_matrix[,colSums(is.na(aflp_matrix))<nrow(aflp_matrix)]
  
  # remove coulumn if specified
  if(!is.null(remove)){
    print(paste("Removing column", remove_col))
    aflp_matrix <- aflp_matrix[,-c(remove)]
    }

  # Remove duplicated rows (as in Structure input files)
  aflp_matrix <- aflp_matrix[!duplicated(aflp_matrix), ]

 

  # Remove rows (samples) containing NA in AFLP data
  aflp_matrix[aflp_matrix==-9] <- NA
  aflp_matrix <- na.omit(aflp_matrix)

  # Remove columns containing the same values
  cols1 <- colSums(aflp_matrix[, -1])[colSums(aflp_matrix[, -1]) == dim(aflp_matrix)[1]]
  cols1 <- names(cols1)
  cols0 <- colSums(aflp_matrix[, -1])[colSums(aflp_matrix[, -1]) == 0]
  cols0 <- names(cols0)
  print("Removing columns with monomorphic markers:")
  print(c(cols1, cols0))
  if(length(cols1)>0){
    aflp_matrix <- aflp_matrix[,!(names(aflp_matrix) %in% cols1)]
  }
  if(length(cols0)>0){
    aflp_matrix <- aflp_matrix[,!(names(aflp_matrix) %in% cols0)]
  }
  
  return(aflp_matrix)
}