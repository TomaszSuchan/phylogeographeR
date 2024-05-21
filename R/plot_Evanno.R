#' Evanno plot function
#' 
#' Plots Evanno delta K plot from structure output data, ggplot-based so can be styled with additional ggplot arguments, eg. plot_Evanno + theme(axis.text = element_text(size = 10))
#' @import ggplot2
#' @import pophelper
#' @param structure_output_path path for the directory with structure output files
#' @export
#' @examples
#' structure_output_path <- "struct-new/output/Pmi"
#' plot_Evanno(structure_output_path)

plot_Evanno <- function(structure_output_path) {
  # Structure
  sfiles <- list.files(structure_output_path, full.names = T)
  slist <- readQ(sfiles, filetype="structure", indlabfromfile = T)
  # runs summary:
  stats_str <- summariseQ(tabulateQ(slist))
  evanno <- evannoMethodStructure(data=stats_str, exportplot=F,returnplot=F,returndata=T)
  # plot
  evanno_plot <- ggplot(data=evanno[2:(nrow(evanno)-1),], aes(x = k, y = deltaK)) +
                 geom_line() +
                 geom_point() +
                 scale_x_continuous(breaks=evanno$k)
  
  return(evanno_plot)
}
