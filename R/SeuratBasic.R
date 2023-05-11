#' Make a Seurat object from a data frame.
#' @description Make a Seurat object and process through dimension reduction steps from a data frame. Includes a very basic data filtering and dimension selection option.
#' Useful for a cursory analysis. It is highly recommended that single-cell is not cleaned through this function, but rather through careful consideration using dedicated
#' single-cell data processing packages (e.g. Seurat, SingleCellExperiment).
#'
#' @param sc a data frame with cells as columns and features as rows
#' @param nfeatures Number of features to select as top variable features
#' @param dims dimensions to use as input features
#' @param Feature_low Include cells where at least this many features are detected.
#' @param Feature_high Exclude cells where more than this many features are detected.
#'
#' @return Returns a Seurat object containing a UMAP representation
#'
#' @import Seurat
#'
#' @export
#'
SeuratBasic <- function(sc,nfeatures = 3000,dims = 1:30,
                        Feature_low = 500,
                        Feature_high = 10000){
  nFeature_RNA <- NULL
  sc <- CreateSeuratObject(sc)
  sc <- subset(sc, subset = nFeature_RNA > Feature_low)
  sc <- subset(sc, subset = nFeature_RNA < Feature_high)
  sc <- NormalizeData(sc)
  sc <- FindVariableFeatures(sc, nfeatures = nfeatures)
  sc <- ScaleData(sc)
  sc <- RunPCA(sc)
  sc <- RunUMAP(sc, dims = dims)
  sc
}
