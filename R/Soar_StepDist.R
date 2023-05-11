#' @title Downstream visualizations: Expression step change (enthalpy) across each bin
#' @description Downstream visualizations step change across each bin, showing the nearest pseudotime bins.
#' Useful as a more direct visualization compared to the pyramid plots.
#' @param egret  An egret object
#' @param Ethreshold a threshold value depicted as a line in the plot. Useful for pointing out potential points of discontinuity.
#' Outliers can be selected acocrding to user choice. As one possibility, values over 3 median absolute deviations from the median can be useful as a guide.
#' @return null
#' @export
#'
Soar_StepDist <- function(egret, Ethreshold = 0.5){
  Pseudostep <- Relative_Step_Change <- CellNumber <- NULL
  ggplot(egret@H_Summary, aes(Pseudostep, Relative_Step_Change,
                              fill = CellNumber, group=1))+
    geom_point(size =5, shape =21, alpha = 0.8, stroke = 1.5)+
    geom_line(size = 2, alpha = 0.6)+
    geom_hline(yintercept = Ethreshold, color = 'firebrick', alpha = 0.8)+
    theme_classic()+theme(text = element_text(size =18, face = 'bold'))
}
