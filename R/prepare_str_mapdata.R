#' Prepare AFLP data for mapmixture
#' 
#' Prepare data with pophelper for plotting georgraphical map with structure piecharts
#' @name prepare_str_mapdata
#' @import pophelper
#' @param structure_output_path path for the directory with structure output files
#' @param structure_input_path path for the structure input file or plink .fam file
#' @param population_data_path path for the population data file where the first column is the population, second y-coordinate and third x-coordinate. The rest of the columns are ignored. Can have a header as long it does not contain any of the population names!
#' @param filetype type of the input data ('auto', 'structure','tess2','baps','faststructure' or 'clumpp')
#' @param remove vector of colums to be removed
#' @export
#' @examples
#' # For faststructure data:
#' population_data_path <- "test_data/fastStructure/Pta_population_data.txt"
#' structure_output_path <- "test_data/fastStructure/meanQ"
#' structure_input_path <- "test_data/fastStructure/Pta_plink.fam"
#' strmapdata <- prepare_str_mapdata(structure_output_path, structure_input_path, population_data_path, filetype="basic")
#' boundary <- c(xmin=18, xmax=21, ymin=48.5, ymax=49.5)
#' earth <- terra::rast("~/GIS/NE1_HR_LC_SR_W_DR/NE1_HR_LC_SR_W_DR.tif")
#' gg_color_hue <- function(n) {
#'   hues = seq(15, 375, length = n + 1)
#'   hcl(h = hues, l = 65, c = 100)[1:n]
#' }
#' k <- 3
#' mapmixture(strmapdata$str[[as.character(k)]], strmapdata$pop, pie_size = 0.1, basemap = earth, boundary = boundary, cluster_cols = gg_color_hue(k), arrow_position = "tr")
#' 
#' # For structure data:
#' population_data_path <- "test_data/structure/popdata.txt"
#' structure_output_path <- "test_data/structure/Aal_admix_out"
#' structure_input_path <- "test_data/structure/Aal_carp_structure-inputD.txt"
#' strmapdata <- prepare_str_mapdata(structure_output_path, structure_input_path, population_data_path, filetype="structure")
#' boundary <- c(xmin=16, xmax=27.5, ymin=43.5, ymax=50)
#' earth <- terra::rast("~/GIS/NE1_HR_LC_SR_W_DR/NE1_HR_LC_SR_W_DR.tif")
#' gg_color_hue <- function(n) {
#'   hues = seq(15, 375, length = n + 1)
#'   hcl(h = hues, l = 65, c = 100)[1:n]
#' }
#' k <- 3
#' mapmixture(strmapdata$str[[as.character(k)]], strmapdata$pop, pie_size = 0.5, basemap = earth, boundary = boundary, cluster_cols = gg_color_hue(k), arrow_position = "tr")


prepare_str_mapdata <- function(structure_output_path, structure_input_path, filetype="structure", delimiter="_") {
  # Rename "faststructure" to "basic" for pophelper compatibility
  if(filetype=="faststructure"){
    filetype <- "basic"
  }

  # Structure
  if(filetype=="basic"){
    sfiles <- list.files(structure_output_path, pattern="meanQ", full.names = T)
    slist <- pophelper::readQ(sfiles, filetype=filetype, indlabfromfile = F)
  } else{
    sfiles <- list.files(structure_output_path, full.names = T)
    slist <- pophelper::readQ(sfiles, filetype=filetype, indlabfromfile = F)
  }
 
  #clumpp-like:
  aligned_slist <- pophelper::alignK(slist)
  split_slist <- pophelper::splitQ(aligned_slist)

  if(filetype=="structure"){
    reduceApplyListOfArrays<- function(x){
         y<-apply(array(unlist(x), c(dim(x[[1]]), dim(x[[2]]), length(x))), c(1,2), mean)
         colnames(y)<-colnames(x[[1]])
         rownames(y)<-rownames(x[[1]])
      return(y)
    }
    merged_slist <- lapply(split_slist, reduceApplyListOfArrays)
  } else {
    merged_slist <- split_slist 
  }
  
  
  ind <- read.table(structure_input_path, row.names=NULL)[1]
  ind <- ind[!duplicated(ind), ]
  split_names <- strsplit(as.character(ind), delimiter)

  # check that all individuals have the same number of delimiters
  vals <- sapply(split_names, length)
  if(!all(vals == vals[1])) stop("Inconsistent number of delimiters in individual names. Check the delimiter argument.")

  # handle different cases based on number of delimiters
  if (length(split_names[[1]])==1) {
    stop("Delimiter not found in individual names. Please check the delimiter argument.")
  } else if (length(split_names[[1]])==2) {
    site <- sapply(split_names, function(x) x[1])
  } else if (length(split_names[[1]])>2) {
    site <- sapply(split_names, function(x) paste(x[-length(x)], collapse=delimiter))
  }

  addPopdata <- function(x){
    x <- cbind(site, ind, as.data.frame(x))
    return(x)
  }

  merged_slist_popdata <- lapply(merged_slist, addPopdata)

  return(merged_slist_popdata)
}
