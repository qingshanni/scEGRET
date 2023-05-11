#' @title Pyramid plot showing gene expression correlation matrix
#' @description Downstream visualizations of average gene expression correlation across each pseudotime bin rendered as a pyramid of all possible combinations. 
#' Users can modify colors and color range for optimal visualization. Correlation calculated using Pearson's R and expressed as R value. 
#' @param egret An Egret object
#' @param l.low Minimum value for bins (ideally this value is not greater than the minimum value in any bin).
#' @param l.mid Medium value for bins. Set to define color blending of figure legend.
#' @param l.high Maximum value for bins.
#' @param c.low Low color for bins
#' @param c.mid Middle color for bins
#' @param c.high High color for bins
#'
#' @import plyr
#' @importFrom stats cor
#' @return null
#' @export
#'
Soar_Mountain_GEXCor <- function(egret, l.low = 0.92, l.mid = 0.97, l.high = 1,
                                 c.low = 'blue', c.mid = 'white', c.high ='red'
){
  PseudoStep <- Var1 <- value <- NULL
  c1 <- cor(egret@AvgExp)

  get_tri <- function(matrix){
    matrix[upper.tri(matrix)]<- NA
    matrix <- melt(matrix, na.rm = TRUE)
    return(matrix)
  }
  lc1 <- get_tri(c1)
  colnames(lc1)[2] <- 'PseudoStep'

  g <- ggplot(data = lc1, aes(PseudoStep, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = c.low, high = c.high, mid = c.mid,
                         midpoint = l.mid, limit = c(l.low,l.high), space = "Lab",
                         name="Expression\nPearson\nCorrelation") +
    theme_classic()+
    scale_y_reverse()+
    scale_x_reverse()+
    theme(axis.text.x = element_text(angle = 225, vjust = 1,
                                     size = 16, hjust = 1),
          axis.title.x = element_text(angle = 180, size = 18),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank())+
    coord_fixed()

  p.legend <- get_legend(g)
  g+theme(legend.position = 'none')

  grid.newpage()                                 # Create new plot page
  print(g+theme(legend.position = 'none'),          # Draw rotated plot
        vp = viewport(width = 0.8,
                      height = 0.8,
                      angle = 135))
  vp = viewport(x=0.2, y=0.85, width=0, height=0)
  pushViewport(vp)
  grid.draw(p.legend)
}
