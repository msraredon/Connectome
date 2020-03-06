#' CellCellScatter
#'
#' Scatter plot of a single cell-cell vector, involving features of interest. Returns a ggplot object, colored by mechanism. Points are labeled if above label.threshold on both x- and y-axes.
#'
#' @param connectome A connectomic edgelist
#' @param cell.1 The source cell
#' @param cell.2 The receiving cell
#' @param label.threshold Threshold for labeling of plot, applied to both x- and y- axes. Defaults to 1.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc','weight_norm', or 'weight_raw'. Defaults to 'weight_sc'

#' @export

CellCellScatter <- function(connectome,
                            cell.1,
                            cell.2,
                            label.threshold = 1,
                            weight.attribute = 'weight_sc'){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)
  require(ggrepel)
  # Subset
  cell_cell <- subset(connectome, source == cell.1 & target == cell.2)
  # Plot
  if (weight.attribute == 'weight_sc'){
    p1 <- ggplot(cell_cell,aes(recept.scale,ligand.scale,size = weight_sc,color = pair))+
        ggrepel::geom_text_repel(data=subset(cell_cell, recept.scale > label.threshold & ligand.scale > label.threshold),aes(label=pair))
  }
  if (weight.attribute == 'weight_norm'){
    p1 <- ggplot(cell_cell,aes(recept.expression,ligand.expression,size = weight_norm,color = pair))+
        ggrepel::geom_text_repel(data=subset(cell_cell, recept.expression > label.threshold & ligand.expression > label.threshold),aes(label=pair))
  }
  if (weight.attribute == 'weight_raw'){
    p1 <- ggplot(cell_cell,aes(recept.raw,ligand.raw,size = weight_raw,color = pair))+
        ggrepel::geom_text_repel(data=subset(cell_cell, recept.raw > label.threshold & ligand.raw > label.threshold),aes(label=pair))
  }
  p1 <-  p1 +
    geom_point() +
    theme_light()+
    guides(colour = guide_legend(override.aes = list(size=6)))+
    ggtitle(paste(cell.1,' to ',cell.2,'Signaling'))+
    theme(legend.position="right")

  return(p1)
}
