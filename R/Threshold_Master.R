#' Threshold_Master
#'
#' This function thresholds a connectome to only include those edges which meet the input parameters
#'
#' @param connectome A connectomic edgelist
#' @param pct.thresh Percent of cluster expressing marker
#' @param p.thresh P-value cutoff, defaults to 0.05
#' @param scale.thresh Scale value cutoff, defaults to NULL
#' @export

Threshold_Master <- function(connectome,pct.thresh = 0.05,p.thresh = 0.05,scale.thresh = NULL){
  temp <- subset(connectome,master.rec.pct > pct.thresh & master.lig.pct > pct.thresh & lig.p.fisher < p.thresh & rec.p.fisher < p.thresh)
  if(!is.null(x = scale.thresh)){
    temp <- subset(connectome,master.rec.pct > pct.thresh & master.lig.pct > pct.thresh & lig.p.fisher < p.thresh & rec.p.fisher < p.thresh & master.lig.wt > scale.thresh & master.rec.wt > scale.thresh)
  }
  return(temp)
}
