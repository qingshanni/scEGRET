#' @title Make a CellDataSet object from a Seurat object
#'
#' @description Make a CellDataSet object used in monocle from a Seurat object to proceed to trajectory analysis
#'
#' @param sc a Seurat object
#' @param lowerDetectionLimit the minimum expression level that consistitutes true expression
#'
#' @return Returns a CellDataSet object used in monocle
#' @export
#'
seurat2Monocle <- function(sc,lowerDetectionLimit = 0.01){
  fdata           <-  data.frame(gene_short_name = rownames(sc@assays$RNA@counts))
  rownames(fdata) <- rownames(sc@assays$RNA@counts)
  fd <- new('AnnotatedDataFrame', data = fdata)
  pd <- new('AnnotatedDataFrame', data = sc@meta.data)
  newCellDataSet(sc@assays$RNA@counts,phenoData = pd, featureData = fd,
                 lowerDetectionLimit = lowerDetectionLimit)
}
