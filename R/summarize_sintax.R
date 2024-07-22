#' Combine OTU table with sintax classification
#' 
#' Outputs a table with sintax classification plus sequence abundances from the OTU table 
#' @param sintax file with the output of sintax classification
#' @param otu_table a tex table with samples as columns and OTUs as rows
#' @export
#' @examples
#' sintax <- "test_data/metabarcoding/ITS2-otus.sintax-UNITE.txt"
#' otu_table <- "test_data/metabarcoding/ITS2-otutab.txt"
#' summarize_sintax(sintax, otu_table)

summarize_sintax <- function(sintax, otu_table) {
  otutab <- read.table(otu_table, header=TRUE, sep="\t", fill = TRUE, comment.char = "")
  colnames(otutab)[1] <- "OTU"
  classification <- read.table(sintax, header=F, sep="\t", fill = TRUE)
  classification$V3 <- NULL
  
  otu <- data.frame(do.call('rbind', strsplit(as.character(classification$V1),';',fixed=TRUE)))[1]
  
  classification <- tidyr::separate_wider_delim(classification, cols='V2', names=c('d', 'p', 'c', 'o', 'f', 'g', 's'), delim= ",", too_few = 'align_start', too_many = 'drop')
  classification <- tidyr::separate_wider_delim(classification, cols='V1', names=c('OTU', 'size'), delim= ";size=", too_few = 'align_start')

  classification <- tidyr::separate_wider_delim(classification, cols='d', names=c('d', 'd_prob'), delim= "(", too_few = 'align_start') 
  classification <- tidyr::separate_wider_delim(classification, cols='p', names=c('p', 'p_prob'), delim= "(", too_few = 'align_start') 
  classification <- tidyr::separate_wider_delim(classification, cols='c', names=c('c', 'c_prob'), delim= "(", too_few = 'align_start') 
  classification <- tidyr::separate_wider_delim(classification, cols='o', names=c('o', 'o_prob'), delim= "(", too_few = 'align_start') 
  classification <- tidyr::separate_wider_delim(classification, cols='f', names=c('f', 'f_prob'), delim= "(", too_few = 'align_start') 
  classification <- tidyr::separate_wider_delim(classification, cols='g', names=c('g', 'g_prob'), delim= "(", too_few = 'align_start') 
  classification <- tidyr::separate_wider_delim(classification, cols='s', names=c('s', 's_prob'), delim= "(", too_few = 'align_start') 

  classification$d <- gsub("d:", "", classification$d)
  classification$p <- gsub("p:", "", classification$p)
  classification$c <- gsub("c:", "", classification$c)
  classification$o <- gsub("o:", "", classification$o)
  classification$f <- gsub("f:", "", classification$f)
  classification$g <- gsub("g:", "", classification$g)
  classification$s <- gsub("s:", "", classification$s)
  
  classification$d_prob <- gsub(")", "", classification$d_prob)
  classification$p_prob <- gsub(")", "", classification$p_prob)
  classification$c_prob <- gsub(")", "", classification$c_prob)
  classification$o_prob <- gsub(")", "", classification$o_prob)
  classification$f_prob <- gsub(")", "", classification$f_prob)
  classification$g_prob <- gsub(")", "", classification$g_prob)
  classification$s_prob <- gsub(")", "", classification$s_prob)

  classification$size <- as.numeric(classification$size)
  classification$d_prob <- as.numeric(classification$d_prob)
  classification$p_prob <- as.numeric(classification$p_prob)
  classification$c_prob <- as.numeric(classification$c_prob)
  classification$o_prob <- as.numeric(classification$o_prob)
  classification$f_prob <- as.numeric(classification$f_prob)
  classification$g_prob <- as.numeric(classification$g_prob)
  classification$s_prob <- as.numeric(classification$s_prob)
  
  aggregated <- merge(classification, otutab, by="OTU")

  return(aggregated)
}