#' @title  Use monocle to generate trajectory
#'
#' @description Use monocle2 DDRTree algorithm to generate trajectory and perform pseudotime assignments.
#'
#' @param mono a monocle2 CellDataSet object
#' @param varGenes a vector of feature ids used for ordering cells. Must be a subset of all features detected.
#' @param reverse Whether to reverse the beginning and end points of the learned biological process. Helpful for cases where the trajectory returned is oriented opposite to
#' actual biology.
#'
#' @import monocle
#' @import ggplot2
#' @importFrom  BiocGenerics estimateSizeFactors
#' @importFrom  BiocGenerics estimateDispersions
#'
#' @return an updated CellDataSet object, in which phenoData contains values for State and Pseudotime for each cell
#' @export
#'
RunMonocle2 <- function(mono, varGenes, reverse = F){
  ###can be taken from sc@assays$RNA@var.features
  mono <- setOrderingFilter(mono, ordering_genes = varGenes)
  mono <- estimateSizeFactors(mono)
  mono <- estimateDispersions(mono)
  mono <- reduceDimension(mono, reduction_method = 'DDRTree',
                          auto_param_selection = F)
  mono <- orderCells(mono, reverse = reverse)
  mono
}
