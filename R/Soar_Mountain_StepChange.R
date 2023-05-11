#' @title  Pyramid plot showing gene expression enthalpy change across each bin
#' @description Downstream visualizations of the overall gene expression difference across each pseudotime bin rendered as a pyramid of all possible combinations. 
#' Users can modify colors and color range for optimal visualization. Correlation calculated as the sum of all variation in gene expression. May be useful to conceptualize 
#' as being similar to thermodynamic reaction enthalpy. 
#' @param egret An Egret object
#' @param l.low Minimum value for bins (ideally this value is not greater than the minimum value in any bin).
#' @param l.mid Medium value for bins. Set to define color blending of figure legend.
#' @param l.high Maximum value for bins.
#' @param c.low Low color for bins
#' @param c.mid Middle color for bins
#' @param c.high High color for bins
#'
#' @return null
#' @export
#'
Soar_Mountain_StepChange <- function(egret, l.low = 0, l.mid = 0.25, l.high = 0.5,
                                     c.low = 'blue', c.mid = 'white', c.high ='red')
{
  PseudoStep <- S2 <- NULL
  cobm <- egret@AllCombns
  cobm <- as.data.frame(t(cobm))
  stepc <-  egret@stepChangefull
  egret.mean <- rowMeans(as.matrix(egret@FilteredCounts))
  egret.mean <- log(1+egret.mean)
  egret.mean <- as.vector(egret.mean)
  stepc <- as.data.frame(t(abs(stepc)))

  stepc <- lapply(1:ncol(stepc), function(x)
    stepc[,x]/egret.mean[x])
  stepc <- as.data.frame(stepc)
  stepc <- rowSums(stepc)

  t1 <- as.data.frame(stepc)
  t1$stepc <- t1$stepc/nrow(egret@FilteredCounts)
  t1$S1 <- cobm$V1
  t1$S2 <- cobm$V2

  #write a blank matrix to fill zeros
  nguide <- union(levels(as.factor(t1$S1)), levels(as.factor(t1$S2)))
  sup1 <- as.data.frame(matrix(0, length(nguide), 3))
  colnames(sup1)[1] <- 'stepc'
  colnames(sup1)[2] <- 'S1'
  colnames(sup1)[3] <- 'S2'
  sup1$S1 <- nguide
  sup1$S2 <- nguide
  t1 <- rbind(t1, sup1)


  colnames(t1)[2] <- 'PseudoStep'
  t1$PseudoStep <- as.numeric(t1$PseudoStep)
  t1$S2 <- as.numeric(t1$S2)

  g <- ggplot(data = t1, aes(PseudoStep, S2, fill = stepc))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = c.low, high = c.high, mid = c.mid,
                         midpoint = l.mid, limit = c(l.low,l.high), space = "Lab",
                         name="Relative\nStep\nChange") +
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
