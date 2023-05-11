#' @title Find trajectory-independent background genes
#' @description Find independent background genes that do not significantly change in expression relative to trajectory progression.
#' Identification of these background genes may be helpful for identifiying context-specific reference genes to use as controls for
#' downstream validation work such as qPCR.
#' @param egret An Egret object
#' @param sigDet Cutoff for detection rate of a feature across the dataset being considered. For stably expressed reference genes, 
#' higher detection rates are advisable for picking out genes for more reliable validation. Default = 0.3 (30% of all cells)
#' @param sigExp Cutoff for average detection rate of a feature across the dataset being considered. Genes with higher average detection are likely to be more easily detected. 
#' Default cutoff =2 (log2count)
#' @param sigVar Variance cutoff for expression variance of a feature across an entire dataset. Calculated as standard deviation of log2Count.
#' @param sigPval p-value cutoff for expression varaible feature detection. 
#' Features with p-values < threshold are excluded as possible variants regardless of actual difference.
#' @param sigDiff High threshold for maximum expression vairance between any two bins (log2 fold change) that would rule out a feature as being potentially unstable.
#'
#' @return t1
#' @export
#'
FindBackgroundGenes <- function(egret,
                                sigDet = 0.3,
                                sigExp = 2,
                                sigVar = 10,
                                sigPval = 0.05,
                                sigDiff = 0.3){

  bindegs <- egret@BinDEGs
  bindegs <- bindegs[bindegs$PVal < sigPval, ]
  bindegs <- bindegs[bindegs$StepChange < sigDiff, ]

  t1 <- egret@GeneInfo
  sexpgene <- setdiff(row.names(t1), bindegs$Gene)
  t1 <- t1[sexpgene, ]
  t1 <- t1[t1$DetectionRate > sigDet, ]
  t1 <- t1[t1$Mean > sigExp, ]
  t1 <- t1[t1$SD < sigVar, ]
  t1
}
