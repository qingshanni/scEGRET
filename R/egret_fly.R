#' @title Set flight. Calculates gene expression variation with respect to reference trajectory
#' @description Fly calculates GEX changes with respect to reference trajectory order in terms of three seperate metrics, encompassing number of differentially expressed genes #' between each pseudotime step, mean expression profile similarity between each pseduotime step, and expression enthalpy between each pseudotime step. Calculated values are
#' stored in associated slots of the Egret object returned as output.
#' @param egret An egret object
#' @param steps Number of pseudotime steps to bin cells into (integer). Larger values will lead to more finely delinated pseudotime steps, but may also introduce  excessive noise caused by technical variation in sequencing in small numbers of cells. As an approximate guide, the Rice rule (2*cube root of cell count) for histogram binning may be helpful; 1000 cells corresponds to 20 bins. Default value = 10.
#' @param minCellCount Minimum number of cells required to define a discrete pseudotime bin. As cells in the same bin are merged together to generate an average expression profile, it is recommended for each bin to include several cells to limit technical noise. Sensible values may range from 3-100s based on the detection sensitivty of the underlying single-cell library construction method and known dropout rate. Default value = 3.
#' @param min.detect Minimum detection rate cutoff for a feature to be considered in these analyses. Introducing too many rarely detected genes may introduce technical noise. Default value of 0.1 filters out features expressed in less than 10% of all cells.
#' @param min.pval.cutoff Minimium cutoff p-value for differential expression analysis. Stricter or looser thresholds set the number of DEGs stored. Default t-test p-value = 0.05. Recommended to set a looser threshold for initial calculation, and modify downstream for specific visualizations.
#'
#' @import matrixTests
#' @import matrixStats
#' @import reshape2
#' @import utils
#' @return an egret object
#' @export
#'
egret.fly <- function(egret,
                      steps = 10,
                      minCellCount = 3,
                      min.detect = 0.1,
                      min.pval.cutoff = 0.05){
  Pseudotime <- rowSds <- NULL
  egret@Metadata$Pseudotime <- egret@Pseudotime
  egret@Metadata <- egret@Metadata  %>% mutate(BinInterval = cut(Pseudotime, breaks= steps))
  ###bins cells according to pseudotime step
  egret@Metadata$TimeBin <- as.factor(egret@Metadata$BinInterval)
  levels(egret@Metadata$TimeBin) <- 1:length(levels(egret@Metadata$TimeBin))
  egret@Metadata <- droplevels.data.frame(egret@Metadata)

  cellbins <- lapply(levels(egret@Metadata$TimeBin), function(x)
    row.names(egret@Metadata[egret@Metadata$TimeBin == x, ]))
  names(cellbins) <- levels(egret@Metadata$TimeBin)
  cellbincount <- lapply(1:length(cellbins), function(x)
    length(cellbins[[x]]))
  cellbins <- cellbins[cellbincount > minCellCount]
  egret@CellBins <- cellbins

  ###filter genes by detection rate
  detrate <- rowSums(egret@Counts != 0)
  detrate <- as.data.frame(detrate)
  colnames(detrate)[1] <- "Nonzero_Count"
  detrate$DetectionRate <- detrate$Nonzero_Count / ncol(egret@Counts)
  detrate <- detrate[detrate$DetectionRate > min.detect, ]

  egret@FilteredCounts <- egret@Counts[row.names(detrate), ]
  egret@GeneInfo <- detrate
  egret.sd <- rowSds(as.matrix(egret@FilteredCounts))
  egret.sd <- log(1+egret.sd)
  egret@GeneInfo$SD <- egret.sd
  egret.mean <- rowMeans(as.matrix(egret@FilteredCounts))
  egret.mean <- log(1+egret.mean)
  egret@GeneInfo$Mean <- egret.mean

  ###average gene expression in each bin
  avgexp <- lapply(egret@CellBins, function(x) rowMeans(egret@FilteredCounts[x]))
  avgexp <- as.data.frame(avgexp)
  avgexp <-log(1+avgexp)
  names(avgexp) <- names(cellbins)
  egret@AvgExp <- avgexp


  ###calculate step enthalpy
  stepEXP <- lapply(2:ncol(avgexp), function(x) avgexp[[x]] - avgexp[[x-1]])
  stepEXP <- as.data.frame(stepEXP)
  egret@stepChange <- stepEXP
  names(egret@stepChange) <- 1:(ncol(avgexp)-1)
  stepEXP <- abs(stepEXP)

  egret.enth <- lapply(1:nrow(stepEXP), function(x)
    stepEXP[x,]/egret.mean[x])
  egret.enth <- as.data.frame(do.call(rbind, egret.enth))
  egret.enth <- egret.enth[!is.na(egret.enth[1]), ]
  sum.enth <- colSums(egret.enth)
  sum.enth <- as.data.frame(sum.enth)
  row.names(sum.enth) <- 1:nrow(sum.enth)
  colnames(sum.enth)[1] <- 'StepEnthalpy'

  cellcount <- as.data.frame(table(egret@Metadata$TimeBin))
  cellcount <- cellcount[cellcount$Freq > minCellCount, ]
  cellcount <- cellcount$Freq
  adjcount <- lapply(2:length(cellcount),
                     function(x) cellcount[[x]] - cellcount[[x-1]])
  sum.enth$CellNumber <- unlist(adjcount)
  sum.enth$Pseudostep <- as.numeric(row.names(sum.enth))
  sum.enth$Relative_Step_Change <- sum.enth$StepEnthalpy / nrow(egret@FilteredCounts)
  egret@H_Summary <- sum.enth

  ###calculate DEGs count per bin
  #using stepEXP
  wilcoxtable <- lapply(2:length(cellbins), function(x)
    row_wilcoxon_twosample(egret@FilteredCounts[cellbins[[(x-1)]]],
                           egret@FilteredCounts[cellbins[[x]]],
                           correct = T))
  wilcoxtable <- lapply(1:(length(cellbins)-1), function(x)
    wilcoxtable[[x]]$pvalue)
  egret@WilcoxPval <- as.data.frame(wilcoxtable)
  names(egret@WilcoxPval) <- 1:(ncol(avgexp)-1)

  ##assemble a df for significant bin-wise GEX diff
  ex1 <- egret@stepChange
  pv1 <- egret@WilcoxPval
  ex1$Gene <- row.names(egret@GeneInfo)
  ex1 <- melt(ex1, id.vars = 'Gene')
  pv1$Gene <- row.names(egret@GeneInfo)
  pv1 <- melt(pv1, id.vars = 'Gene')
  colnames(ex1)[2] <- 'PseudoStep'
  colnames(ex1)[3] <- 'StepChange'
  ex1$PVal <- pv1$value
  ex1 <- ex1[ex1$PVal < min.pval.cutoff, ]
  egret@BinDEGs <- as.data.frame(ex1)

  #######
  ###additional M by M calculations
  ###to generate pyramid plots for each possible
  ###combination of 2 bins
  #######
  ###calculate DEGs count per bin
  cellbins <- egret@CellBins
  posscomb <- combn(names(cellbins), 2)
  posscomb <- as.data.frame(posscomb)

  egret@AllCombns <- posscomb

  wilcoxtable2 <- lapply(1:ncol(posscomb), function(x)
    row_wilcoxon_twosample(egret@FilteredCounts[cellbins[[posscomb[[x]][1]]]],
                           egret@FilteredCounts[cellbins[[posscomb[[x]][2]]]],
                           correct = T))
  wilcoxtable2 <- lapply(1:ncol(posscomb), function(x)
    wilcoxtable2[[x]]$pvalue)
  wilcoxtable2 <- as.data.frame(wilcoxtable2)
  names(wilcoxtable2) <- names(posscomb)
  egret@WilcoxPvalfull <- as.data.frame(wilcoxtable2)

  ##
  ###calculate step enthalpy
  avgexp <- egret@AvgExp
  stepEXP2 <- lapply(1:ncol(posscomb), function(x) avgexp[,(posscomb[x][1,])] -
                       avgexp[,(posscomb[x][2,])])
  stepEXP2 <- as.data.frame(stepEXP2)
  names(stepEXP2) <- names(posscomb)

  egret@stepChangefull <- stepEXP2
  stepEXP2 <- abs(stepEXP2)

  egret.enth2 <- lapply(1:nrow(stepEXP2), function(x)
    stepEXP2[x,]/egret.mean[x])
  egret.enth2 <- as.data.frame(do.call(rbind, egret.enth2))
  egret.enth2 <- egret.enth2[!is.na(egret.enth2[1]), ]
  sum.enth2 <- colSums(egret.enth2)
  sum.enth2 <- as.data.frame(sum.enth2)
  row.names(sum.enth2) <- 1:nrow(sum.enth2)
  colnames(sum.enth2)[1] <- 'StepEnthalpy'
  egret@H_Summaryfull <- as.data.frame(sum.enth2)

  ##assemble a df for significant bin-wise GEX diff
  ex1 <- egret@stepChangefull
  pv1 <- egret@WilcoxPvalfull
  ex1$Gene <- row.names(egret@GeneInfo)
  ex1 <- melt(ex1, id.vars = 'Gene')
  pv1$Gene <- row.names(egret@GeneInfo)
  pv1 <- melt(pv1, id.vars = 'Gene')
  colnames(ex1)[2] <- 'PseudoStep'
  colnames(ex1)[3] <- 'StepChange'
  ex1$PVal <- pv1$value
  ex1 <- ex1[ex1$PVal < min.pval.cutoff, ]
  egret@BinDEGsfull <- as.data.frame(ex1)

  ###
  #object out
  ###
  egret
}
