#' ModeScatter
#'
#' This function takes a connectome and plots the interactions within a single signaling modality. x-axis is receptor expression and y-axis is ligand expression. Color is exact mechanism, label is cell-cell vector.
#'
#' @param connectome A connectomic edgelist
#' @param MOI The mode of interest (i.e., 'Wnt')
#' @param cell2 The receiving cell
#' @param weight The column name for the edgeweight of interest (determines dot sizes). Defaults to weight_sc
#' @param rec.wt The column name for receptor expression (x axis). Defaults to recept.scale
#' @param lig.wt The column name for ligand expression (y axis). Defaults to ligand.scale
#' @export
ModeScatter <- function(connectome,MOI = 'Wnt',weight = 'weight_sc',rec.wt = 'recept.scale',lig.wt = 'ligand.scale'){
  require(ggrepel)
  cell_cell <- subset(connectome, mode == MOI)
  cell_cell$cellpair <- paste(cell_cell$source,cell_cell$target,sep='-')
  ggplot(cell_cell,aes_string(rec.wt,lig.wt,size = weight,color = 'pair')) +
    geom_jitter() +
    theme_light()+
    guides(colour = guide_legend(override.aes = list(size=6)))+
    ggtitle(paste(MOI,'Signaling'))+
    theme(legend.position="right")+
    geom_text_repel(aes(label=cellpair),hjust=0, vjust=0,size = 6)
}
