#' @title number of DEGs between each bin
#' @description number of DEGs between each bin calculated using pairwise comparsion across each two given pseudotime bins.
#' @param egret An Egret object
#' @param sigPval p-value cutoff used to judge DEG signficance
#' @param sigDiff Expression fold change (log2FC) cutoff used to determine DEG significance
#'
#' @return a plot
#' @export
#'
Soar_DEGs <- function(egret, sigPval = 0.01 ,sigDiff = 0.25){
  Pseudostep <- DEGs <- CellNumber <- NULL
  bindegs <- egret@BinDEGs
  bindegs <- bindegs[bindegs$PVal < sigPval, ]
  bindegs <- bindegs[bindegs$StepChange > sigDiff, ]

  posdeg <- bindegs[bindegs$StepChange > 0, ]
  negdeg <- bindegs[bindegs$StepChange < 0, ]

  posdeg <- as.data.frame(table(posdeg$PseudoStep))
  negdeg <- as.data.frame(table(negdeg$PseudoStep))

  egret@H_Summary$DEGs_Pos <- unlist(posdeg$Freq)
  egret@H_Summary$DEGs_Neg <- unlist(negdeg$Freq)

  egret@H_Summary$DEGs <- egret@H_Summary$DEGs_Pos + egret@H_Summary$DEGs_Neg

  ggplot(egret@H_Summary, aes(Pseudostep, DEGs,
                              fill = CellNumber, group=1))+
    geom_point(size =5, shape =21, alpha = 0.8, stroke = 1.5)+
    geom_line(size = 2, alpha = 0.6)+
    theme_classic()+theme(text = element_text(size =18, face = 'bold'))
}
