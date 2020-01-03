#' Threshold
#'
#' This function thresholds a connectome to only include those edges which meet the input parameters
#'
#' @param connectome A connectomic edgelist
#' @param pct.thresh Percent of cluster expressing marker
#' @param p.thresh P-value cutoff, defaults to 0.05
#' @param scale.thresh Scale value cutoff, defaults to NULL
#' @export

Threshold <- function(connectome,pct.thresh = 0.05,p.thresh = 0.05,scale.thresh = NULL){
  temp <- subset(connectome,percent.source > pct.thresh & percent.target > pct.thresh & lig.p < p.thresh & rec.p < p.thresh)
  if(!is.null(x = scale.thresh)){
    temp <- subset(connectome,percent.source > pct.thresh & percent.target > pct.thresh & lig.p < p.thresh & rec.p < p.thresh & ligand.scale > scale.thresh & recept.scale > scale.thresh)
  }
  return(temp)
}
