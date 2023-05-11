#' @title Find trajectory-unstable genes
#' @description Find trajectory-unstable genes that show multiple changes in variation but which may not be neccesarily associated with pseudotime progression.
#' @param egret An Egret object
#' @param sigPval p-value cutoff for expression varaible feature detection. 
#' @param sigDiff Threshold for maximum expression vairance between any two bins (log2 fold change)
#' @import plyr
#' @return  t1
#' @export
#'
FindUnstableGenes <- function(egret, sigPval = 0.01, sigDiff = 0.25){
  bindegs <- egret@BinDEGs
  bindegs <- bindegs[bindegs$PVal < sigPval, ]
  bindegs <- bindegs[bindegs$StepChange > sigDiff, ]

  t1 <- as.data.frame(table(bindegs$Gene))
  colnames(t1)[1] <- 'Gene'
  colnames(t1)[2] <- 'DiffFreq'
  addinf <- egret@GeneInfo
  addinf$Gene <- row.names(addinf)
  addinf <- addinf[addinf$Gene %in% t1$Gene, ]
  t1 <- join(t1, addinf, by = 'Gene')
  t1
}
