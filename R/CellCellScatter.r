#' CellCellScatter
#'
#' This function takes a connectome and two cells and plots the edges based on ligand versus receptor weights
#'
#' @param connectome A connectomic edgelist
#' @param cell.1 The source cell
#' @param cell.2 The receiving cell
#' @param lab.thresh Threshold for labeling of plot, applied to both x- and y- axes. Defaults to 1.
#' @export

CellCellScatter <- function(connectome,cell.1,cell.2,lab.thresh = 1){
  # Subset
  cell_cell <- subset(connectome, source == cell.1 & target == cell.2)
  # Plot
  ggplot(cell_cell,aes(recept.scale,ligand.scale,size = weight_sc,color = pair)) +
    geom_point() +
    theme_light()+
    guides(colour = guide_legend(override.aes = list(size=6)))+
    ggtitle(paste(cell.1,' to ',cell.2,'Signaling'))+
    theme(legend.position="right")+
    geom_text_repel(data=subset(cell_cell, recept.scale > lab.thresh | ligand.scale > lab.thresh),aes(label=pair))
}
