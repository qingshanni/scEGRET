#' @title View Expression of given genes over pseudotime
#' @description View expression of given genes over pseudotime. Useful for manually picked genes (e.g. known markers or functional complexes).
#' potential visualization for trajectory-stable reference and trajectory-unstable noise genes as well.
#' @param egret An Egret object
#' @param pick.features A character vector of feature(s) contained in the gene expression data.
#' @param smooth.span Value controlling the smoothness of the loess span for curve fitting. 
#'
#' @return a plot
#' @export
#'
Soar_GEXpseudotime <- function(egret, pick.features, smooth.span = 0.75){
  Pseudotime <- logExpression <- Gene <- NULL
  pick.df <- egret@FilteredCounts[pick.features, ]
  pick.df <- as.data.frame(t(pick.df))
  pick.df$Pseudotime <- as.numeric(egret@Pseudotime)
  pick.df <- melt(pick.df, id.vars = 'Pseudotime')
  pick.df$Pseudotime <- as.numeric(as.character(pick.df$Pseudotime))
  colnames(pick.df)[2] <- 'Gene'
  colnames(pick.df)[3] <- 'Expression'
  pick.df$logExpression <- log(pick.df$Expression+1, 2)

  ggplot(pick.df, aes(Pseudotime, logExpression, fill = Gene, color= Gene)) +
    geom_point(color = 'black', alpha = 0.6, shape = 21)+
    geom_smooth(method = 'loess', alpha = 0.3, span = smooth.span)+
    theme_classic()+
    theme(text = element_text(size =18, face = 'bold'))
}
