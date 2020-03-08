#' SignalScatter
#'
#' Scatter plot of signaling vectors involving features of interest. Returns a ggplot object, colored by mechanism, labeled by cell-cell vector.
#'
#' @param connectome A connectomic edgelist
#' @param features Features of interest
#' @param ... Arguments passed to FilterConnectome
#' @param label.threshold Threshold for labeling of plot. Sum of receptor score and ligand score must be greater than this number.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc or weight_norm'. Defaults to 'weight_sc'

#' @export

SignalScatter <- function(connectome,
                          features,
                          label.threshold = 1,
                          weight.attribute = 'weight_sc',...){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)

    # Subset
    signal <- FilterConnectome(connectome,remove.na = T,features = features,...)

    # Plot
    if (weight.attribute == 'weight_sc'){
      p1 <- ggplot(signal,aes(recept.scale, ligand.scale, size = weight_sc, color = pair))+
      ggrepel::geom_text_repel(data=subset(signal, (recept.scale + ligand.scale) > label.threshold),aes(label=vector))
    }
    if (weight.attribute == 'weight_norm'){
      p1 <- ggplot(signal,aes(recept.expression, ligand.expression, size = weight_norm, color = pair))+
      ggrepel::geom_text_repel(data=subset(signal, (recept.expression + ligand.expression) > label.threshold),aes(label=vector))
    }

    p1 +
      geom_point() +
      theme_light()+
      guides(colour = guide_legend(override.aes = list(size=6)))+
      theme(legend.position="right")
}
