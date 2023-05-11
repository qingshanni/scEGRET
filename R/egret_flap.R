#' @title Creates an EGRET object
#' @description Creates an EGRET object； takes as input an expression counts table (one feature per row, one cell per column),
#' metadata associated with each cell (e.g. batch, sample origin, treatment, etc.),
#' and a vector of pseudotime values (single value per cell).
#' @param df an expression data frame with with cells as columns and features as rows
#' @param metadata Additional cell-level metadata to add to the egret object. Used for visualizations, and may include both descriptive and numerical variables.
#' Should be a data.frame where the rows are cell names and the columns are
#' additional metadata fields. Row names in the metadata need to match all
#' column names of the counts data.frame; NA values are allowed.
#' @param pseudotimes a vector providing pseudotime values for each cell (calculated externally).
#'
#'
#' @import methods
#' @return  egret.obj
#' @export
#'
egret.flap <- function(df, metadata, pseudotimes){
  setClass("Egret", slots=list(Counts = 'data.frame',
                               FilteredCounts = 'data.frame',
                               GeneInfo = 'data.frame',
                               Metadata = 'data.frame',
                               Pseudotime = 'numeric',
                               CellBins = 'list',
                               H_Summary = 'data.frame',
                               stepChange = 'data.frame',
                               WilcoxPval = 'data.frame',
                               AvgExp  = 'data.frame',
                               BinDEGs = 'data.frame',
                               AllCombns = 'data.frame',
                               WilcoxPvalfull= 'data.frame',
                               stepChangefull= 'data.frame',
                               H_Summaryfull= 'data.frame',
                               BinDEGsfull = 'data.frame'
  ), where=.GlobalEnv)
  egret.obj <- new("Egret", Counts = df,
                   Metadata = metadata,
                   Pseudotime = pseudotimes)
  egret.obj
}
