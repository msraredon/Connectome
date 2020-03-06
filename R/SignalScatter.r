#' SignalScatter
#'
#' Scatter plot of signaling vectors involving features of interest. Returns a ggplot object, colored by mechanism, labeled by cell-cell vector. ##
#'
#' @param connectome A connectomic edgelist
#' @param features Gene of interest. All edges containing these features will be plotted.
#' @param label.threshold Threshold for labeling of plot, applied to both x- and y- axes. Defaults to 1.

#' @export

SignalScatter <- function(connectome,features,label.threshold = 1){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)

  if (!is.null(features)){
    # Subset
    if (length(features)==1){
      cell_cell <- subset(connectome, ligand == features | receptor == features)
    }else{
      cell_cell <- subset(connectome, ligand %in% features | receptor %in% features)
    }
    # Plot
    ggplot(cell_cell,aes(recept.scale, ligand.scale, size = weight_sc, color = pair)) +
      geom_point() +
      theme_light()+
      guides(colour = guide_legend(override.aes = list(size=6)))+
      ggtitle(paste(c(features),'Signaling'))+
      theme(legend.position="right")+
      ggrepel::geom_text_repel(data=subset(cell_cell, recept.scale > label.threshold & ligand.scale > label.threshold),aes(label=vector))
  }
}
