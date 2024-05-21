#' Prepare AFLP data for mapmixture
#' 
#' Prepare data for plotting georgraphical map with structure piecharts from AFLP data using the mapmixture package, ggplot-based
#' @name prepare_str_mapdata
#' @import ggplot2
#' @import pophelper
#' @import mapmixture
#' @param structure_output_path path for the directory with structure output files
#' @param population_data_path path for the population data file where the first column is the population and the second the strata. Can have a header as long it does not contain any of the population names!
#' @param choices vector containing the PCoA axes to be plotted
#' @param remove vector of colums to be removed
#' @export
#' @examples
#' population_data_path <- "population_data.txt"
#' structure_output_path <- "struct-new/output/Pmi"
#' strmapdata <- prepare_str_mapdata(structure_output_path, population_data_path)
#' boundary <- c(xmin=5, xmax=27.5, ymin=43.5, ymax=50)
#' earth <- terra::rast("~/GIS/NE1_HR_LC_SR_W_DR/NE1_HR_LC_SR_W_DR.tif")
#' gg_color_hue <- function(n) {
#'   hues = seq(15, 375, length = n + 1)
#'   hcl(h = hues, l = 65, c = 100)[1:n]
#' }
#' k <- 2
#' mapmixture(strmapdata$str[[k]], strmapdata$pop, pie_size = 0.5, basemap = earth, boundary = boundary, cluster_cols = gg_color_hue(k), arrow_position = "tr")


prepare_str_mapdata <- function(structure_output_path, population_data_path) {
  # Structure
  sfiles <- list.files(structure_output_path, full.names = T)
  slist <- readQ(sfiles, filetype="structure", indlabfromfile = T)
 
  #clumpp-like:
  aligned_slist <- alignK(slist)
  split_slist <- splitQ(aligned_slist)
  reduceApplyListOfArrays<- function(x){
         y<-apply(array(unlist(x), c(dim(x[[1]]), dim(x[[2]]), length(x))),
 
                c(1,2), mean)
       colnames(y)<-colnames(x[[1]])
       rownames(y)<-rownames(x[[1]])
    return(y)
  }
  merged_slist <- lapply(split_slist, reduceApplyListOfArrays)
  addPopdata <- function(x){
    Ind <- rownames(x)
    split_names <- strsplit(as.character(Ind), "-")
    Site <- sapply(split_names, function(x) x[1])
    x <- cbind(Site, Ind, as.data.frame(x))
    return(x)
  }
  merged_slist_popdata <- lapply(merged_slist, addPopdata)
  
  Ind <- rownames(readQ(sfiles[1], filetype="structure", indlabfromfile = T)[[1]])
  # Extract population information from individual names
  split_names <- strsplit(as.character(Ind), "-")
  Site <- sapply(split_names, function(x) x[1])
  # Crate a dataframe with population information
  population_data <- read.table(population_data_path, header=TRUE)
  populations <- population_data[population_data[,1] %in% Site, c(1,3,4)]
  colnames(populations) <- c("Site", "Lat", "Lon")

  return(list(str = merged_slist_popdata, pop = populations))
}
