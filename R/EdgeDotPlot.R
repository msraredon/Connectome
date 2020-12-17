#' EdgeDotPlot
#'
#' Plotting function to make a DotPlot of edges within a specified network. This function first finds all edges meeting the desired thresholding criteria, and then plots complete information regarding all mechanisms and celltype vectors implicated.
#'
#' @param connectome A connectomic object
#' @param ... Arguments passed to FilterConnectome
#' @export

EdgeDotPlot <- function(connectome,...){

  library(ggplot2)

  #Filter

  temp <- FilterConnectome(connectome,remove.na = T,...)


  # Identify vectors and mechanisms remaining
  vectors <- unique(temp$vector)
  mechanisms <- unique(temp$pair)

  # Subset to these
  temp2 <- subset(connectome,pair %in% mechanisms & vector %in% vectors)

  #Plot

  my_palette <- colorRampPalette(c("blue", "yellow", "red"), alpha=TRUE)(n=399)


  p1 <- ggplot(data = temp2,aes(x = vector, y = pair,
                                #alpha = log(weight_norm+1)
                                )) +
    geom_point(aes(size=log(weight_norm+1),color=log(weight_sc+1),stroke = 0))+
   # theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_gradientn('log(weight_scale + 1)', colors=my_palette)+theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  return(p1)
}
