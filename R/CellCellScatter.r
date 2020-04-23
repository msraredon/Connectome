#' CellCellScatter
#'
#' Scatter plot of a single cell-cell vector, involving features of interest. Returns a ggplot object, colored by mechanism. Points are labeled if above label.threshold on both x- and y-axes.
#'
#' @param connectome A connectomic edgelist
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param ... Arguments passed to FilterConnectome
#' @param label.threshold Threshold for labeling of plot, applied to both x- and y- axes. Defaults to 1.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc or weight_norm'. Defaults to 'weight_sc'

#' @export

CellCellScatter <- function(connectome,
                            sources.include,
                            targets.include,
                            label.threshold = 1,
                            weight.attribute = 'weight_sc',...){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)
  require(ggrepel)

  # Subset
  cell_cell <- FilterConnectome(connectome,remove.na = T,sources.include = sources.include,targets.include = targets.include,...)

  # Plot
  if (weight.attribute == 'weight_sc'){
    p1 <- ggplot(cell_cell,aes(recept.scale,ligand.scale,size = weight_sc,color = pair))+
        ggrepel::geom_text_repel(data=subset(cell_cell, (recept.scale^2 + ligand.scale^2) > label.threshold),aes(label=pair))
  }
  if (weight.attribute == 'weight_norm'){
    p1 <- ggplot(cell_cell,aes(recept.expression,ligand.expression,size = weight_norm,color = pair))+
        ggrepel::geom_text_repel(data=subset(cell_cell, (recept.expression^2 + ligand.expression^2) > label.threshold),aes(label=pair))
  }

  p1 <-  p1 +
    geom_point() +
    theme_light()+
    guides(colour = guide_legend(override.aes = list(size=6)))+
    ggtitle(paste(sources.include,'to',targets.include,'Signaling'))+
    theme(legend.position="right")

  return(p1)
}
