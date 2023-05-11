#' @title Find trajectory-correlated genes
#' @description Find genes with expression changes tightly correlated with a given trajectory.
#' These include genes with relatively stable changes (increase or decrease) that reflect the dominant trajectory changes. 
#' @param egret An Egret object
#' @param sigDet Cutoff for detection rate of a feature across the dataset being considered. For tarjectory-associated genes, 
#' higher detection rates are advisable for picking out genes for more reliable validation. Default = 0.2 (20% of all cells)
#' @param sigExp Cutoff for average detection rate of a feature across the dataset being considered. Genes with higher average detection are likely to be more easily detected. 
#' Default cutoff =0.1 (log2count)
#' @param sigVar Variance cutoff for expression variance of a feature across an entire dataset. Genes with too low variance are likely to have low-magnitude expression changes #' which may be difficult to validate. Calculated as standard deviation of log2Count.
#' @import plyr
#' @return  cor.out
#' @export
#'
FindCorrelatedGenes <- function(egret,
                                sigDet = 0.2,
                                sigExp = 0.1,
                                sigVar = 1){
  ginfo <- egret@GeneInfo
  ginfo <- ginfo[ginfo$DetectionRate > sigDet, ]
  ginfo <- ginfo[ginfo$Mean > sigExp, ]
  ginfo <- ginfo[ginfo$SD > sigVar, ]

  t1 <- egret@AvgExp
  t1 <- t1[row.names(ginfo), ]

  cor.out <- lapply(1:nrow(t1), function(x)
    stats::cor(t(t1[x, ]), as.numeric(names(t1))))
  cor.out <- as.data.frame(cor.out)
  cor.out <- as.data.frame(t(cor.out))
  row.names(cor.out) <- row.names(t1)
  colnames(cor.out)[1] <- 'Pearson_Correlation'
  ginfo$Gene <- row.names(ginfo)
  cor.out$Gene <- row.names(cor.out)
  cor.out <- join(cor.out, ginfo, by = 'Gene')

  cor.out
}
