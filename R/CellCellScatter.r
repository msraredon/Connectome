#' CellCellScatter
#'
#' This function takes a connectome and two cells and plots the interactions based on ligand versus receptor weights
#'
#' @param connectome A connectomic edgelist
#' @param cell1 The source cell
#' @param cell2 The receiving cell
#' @param weight The column name for the edgeweight of interest (determines dot sizes). Defaults to 'weight_sc'
#' @param rec.wt The column name for receptor expression (x axis). Defaults to 'recept.scale'
#' @param lig.wt The column name for ligand expression (y axis). Defaults to 'ligand.scale'
#' @export
CellCellScatter <- function(connectome,cell1,cell2,weight = 'weight_sc',rec.wt = 'recept.scale',lig.wt = 'ligand.scale'){
  require(ggrepel)
  cell_cell <- subset(connectome, source == cell1 & target == cell2)
  ggplot(cell_cell,aes_string(rec.wt,lig.wt,size = weight,color = 'pair')) +
    geom_jitter() +
    theme_light()+
    guides(colour = guide_legend(override.aes = list(size=6)))+
    ggtitle(paste(cell1,' to ',cell2,'Signaling'))+
    theme(legend.position="right")+
    geom_text(data=subset(cell_cell, rec.wt > 1 | lig.wt > 1),aes(label=pair))
}
