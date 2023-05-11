#' @title View Expression of given genes over discrete pseudotime bins
#' @description View expression of given genes over discrete pseudotime bins. Useful for manually picked genes (e.g. known markers or functional complexes).
#' potential visualization for trajectory-stable reference and trajectory-unstable noise genes as well.
#' @param egret An Egret object
#' @param pick.features A character vector of feature(s) contained in the gene expression data.
#' @param smooth.span Value controlling the smoothness of the loess span for curve fitting. 
#' @param jit.width Value controlling jitter width. Increase jitter to help prevent overplotting at certain points.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return null
#' @export
#'
Soar_GEXtimebin <- function(egret, pick.features, smooth.span = 0.75, jit.width = 0.5){
  TimeBin <- logExpression <- Gene <- NULL
  pick.df <- egret@FilteredCounts[pick.features, ]
  pick.df <- as.data.frame(t(pick.df))
  pick.df$TimeBin <- egret@Metadata$TimeBin
  pick.df <- melt(pick.df, id.vars = 'TimeBin')
  pick.df$TimeBin <- as.numeric(as.character(pick.df$TimeBin))
  colnames(pick.df)[2] <- 'Gene'
  colnames(pick.df)[3] <- 'Expression'
  pick.df$logExpression <- log(pick.df$Expression+1, 2)

  ggplot(pick.df, aes(TimeBin, logExpression, fill = Gene, color= Gene))+
    geom_jitter(color = 'black', alpha = 0.6, shape = 21, width = jit.width)+
    geom_smooth(method = 'loess', alpha = 0.3, span = smooth.span)+
    theme_classic()+
    theme(text = element_text(size =18, face = 'bold'))
}
