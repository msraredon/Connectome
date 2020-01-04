#' SignalScatter
#'
#' Scatter plot of signaling vectors involving features of interest.
#'
#' @param connectome A connectomic edgelist
#' @param feature Gene of interest. All edges containing this feature will be plotted.
#' @export

SignalScatter <- function(connectome,feature){
  if (!is.null(feature)){
    # Subset
    cell_cell <- subset(connectome, ligand == feature | receptor == feature)
    # Plot
    ggplot(cell_cell,aes(recept.scale, ligand.scale, size = weight_sc, color = pair)) +
      geom_point() +
      theme_light()+
      guides(colour = guide_legend(override.aes = list(size=6)))+
      ggtitle(paste(feature,'Signaling'))+
      theme(legend.position="right")+
      geom_text_repel(aes(label=vector),hjust=0, vjust=0)
  }
}
