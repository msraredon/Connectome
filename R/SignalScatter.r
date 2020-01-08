#' SignalScatter
#'
#' Scatter plot of signaling vectors involving features of interest.
#'
#' @param connectome A connectomic edgelist
#' @param features Gene of interest. All edges containing these features will be plotted.
#' @export

SignalScatter <- function(connectome,features){
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
      geom_text_repel(aes(label=vector),hjust=0, vjust=0)
  }
}
