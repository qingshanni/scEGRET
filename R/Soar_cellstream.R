#' @title Plot the number of cells in each time bin
#' @description Stream plot showing number of cells belonging to
#' each group within each time bin and the relative composition changes between bins. Useful as an approximate guide of changes, but more precise calculations of compositional
#' changes should be performed using more rigorous statistical methods and visualized in linearized plots to avoid possible distortions. 
#' @param egret an EGRET object
#' @param metadata A factor classifier of cell metadata (e.g. cell type, sample, condition, etc.)
#'
#' @import ggstream
#' @return a plot
#' @export
#'
Soar_cellstream <- function(egret, metadata = 'CellType'){
  TimeBin <- CellCount <- Var2 <- NULL
  t1 <- egret@Metadata
  t1 <- table(t1$TimeBin, t1[,metadata])
  t1 <- as.data.frame(t1)
  colnames(t1)[1] <- 'TimeBin'
  colnames(t1)[3] <- 'CellCount'
  t1$TimeBin <- as.numeric(t1$TimeBin)

  ggplot(t1, aes(TimeBin, CellCount, fill = Var2)) +
    geom_stream(color = 'black', alpha = 0.75)+labs(fill = metadata)+
    theme_classic()+
    theme(text = element_text(size =18, face = 'bold'))
}
