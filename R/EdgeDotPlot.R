#' EdgeDotPlot
#'
#' Plotting function to make a DotPlot of edges within a specified network. This function first finds all edges meeting the desired thresholding criteria, and then plots complete information regarding all mechanisms and celltype vectors implicated.
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param weight.attribute Column to use to define edgeweights for network analysis. 'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'. If 'weight_sc', function will automatically filter at min.z = 0 to remove negative source/sink values.
#' @param min.z Minimum z-score for ligand and receptor.
#' @param ... Arguments passed to FilterConnectome
#' @export

EdgeDotPlot <- function(connectome,
                       weight.attribute = 'weight_sc',
                       min.z = NULL,...){
  library(ggplot2)
  
  #Filter
  if (weight.attribute == 'weight_sc' & is.null(min.z)){
    temp <- FilterConnectome(connectome, remove.na = T,min.z = 0,...)
  }else{
    temp <- FilterConnectome(connectome,remove.na = T,min.z = min.z,...)
  }
  
  # Identify vectors and mechanisms remaining
  vectors <- unique(temp$vector)
  mechanisms <- unique(temp$pair)
  
  # Subset to these
  temp2 <- subset(connectome,pair %in% mechanisms & vector %in% vectors)
  
  #Plot
  p1 <- ggplot(data = temp2,aes(x = vector, y = pair, size = weight_sc,alpha = weight_sc)) +geom_point(color = 'blue',stroke = 0)+
    theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  return(p1)
}
