#' Gene expression data lymphoblastoid cell lines of all 210 unrelated HapMap 
#' individuals from four populations
#'
#' A 210 times 47293 matrix with indviduals along the rows and expression data along the columns
#'
#' The first column corresponds to column number 24266 (with name GI_6005726-S) in the orignal data.
#' Column names give the name of the variables, row names the "name" of the individuals. 
#'
#'
#' @docType data
#' @keywords datasets
#' @name SangerData
#' @usage data(SangerData)
#' @format A numerical matrix  with 210 rows and 47293 variables
#' @source Dataset downloaded from 
#' \url{https://ftp.sanger.ac.uk/pub/genevar/}
#'
#' Reference:
#'
#' Stranger, BE et al (2007): Relative impact of nucleotide and copy number variation on gene expression phenotypes
#' Science, 2007•science.org
#'
"SangerData"