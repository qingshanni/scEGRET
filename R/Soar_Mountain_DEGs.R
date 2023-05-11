#' @title Downstream visualizations of DEGs in each bin
#' @description Downstream visualizations of DEGs across each pseudotime bin rendered as a pyramid of all possible combinations. Users can modify colors and color range for
#' optimal visualization of the range of DEGs between bins. 
#' @param egret An egret object
#' @param sigPval p-value cutoff for DEG analysis
#' @param sigDiff logFC cutoff for DEG analysis
#' @param l.low Minimum value for bins (ideally this value is not greater than the minimum value in any bin).
#' @param l.mid Medium value for bins. Set to define color blending of figure legend.
#' @param l.high Maximum value for bins.
#' @param c.low Low color for bins
#' @param c.mid Middle color for bins
#' @param c.high High color for bins
#'
#' @import reshape2
#' @import cowplot
#' @import grid
#' @return null
#' @export
#'
Soar_Mountain_DEGs <- function(egret, sigPval = 0.01 ,sigDiff = 0.25,
                               l.low = 0, l.mid = 250, l.high = 1500,
                               c.low = 'blue', c.mid = 'white', c.high ='red')
{
  PseudoStep<- Var1 <- value <- NULL
  cobm <- egret@AllCombns
  cobm <- as.data.frame(t(cobm))
  degfull <- egret@BinDEGsfull

  degfull <- degfull[degfull$PVal < sigPval, ]
  degfull <- degfull[abs(degfull$StepChange) > sigDiff, ]

  t1 <- as.data.frame(table(degfull$PseudoStep))
  t1$S1 <- cobm$V1
  t1$S2 <- cobm$V2

  #write a blank matrix to fill zeros
  nguide <- union(levels(as.factor(t1$S1)), levels(as.factor(t1$S2)))
  sup1 <- as.data.frame(matrix(0, length(nguide), 3))
  colnames(sup1)[1] <- 'Freq'
  colnames(sup1)[2] <- 'S1'
  colnames(sup1)[3] <- 'S2'
  sup1$S1 <- nguide
  sup1$S2 <- nguide
  t1 <- t1[2:4]
  t1 <- rbind(t1, sup1)

  t1 <- dcast(t1, S1~S2, value.var = 'Freq')
  t1$S1 <- as.numeric(t1$S1)
  t1 <- t1[order(t1$S1, decreasing = F), ]
  row.names(t1) <- t1$S1
  t1 <- t1[2:ncol(t1)]
  t1 <- as.data.frame(t(t1))
  t1$S2 <- row.names(t1)
  t1$S2 <- as.numeric(t1$S2)
  t1 <- t1[order(t1$S2, decreasing = F), ]
  t1 <- t1[1:(ncol(t1)-1)]
  t1 <- melt(as.matrix(t1), na.rm = TRUE)
  colnames(t1)[2] <- 'PseudoStep'
  t1$Var1 <- as.numeric(t1$Var1)
  t1$PseudoStep <- as.numeric(t1$PseudoStep)

  g <- ggplot(data = t1, aes(PseudoStep, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = c.low, high = c.high, mid = c.mid,
                         midpoint = l.mid, limit = c(l.low,l.high), space = "Lab",
                         name="Number of\nDEGs") +
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
