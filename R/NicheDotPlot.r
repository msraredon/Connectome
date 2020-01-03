# This requires the weight column to be called 'weight_sc'.  Should be revised to be more flexible
#' NicheDotPlot
#'
#' This function takes a connectomic edgelist and creates a dot plot of a given niche, organized by source cell and mode. The y-axis is binned by source cell type,
#' and the x-axis is the sum of the weight of all edges for a given ligand made by each cell.
#' @param connectome A connectomic edgelist
#' @param NOI A single node of interest. This function will plot an illustration of the niche of this node.
#' @export

NicheDotPlot <- function(connectome,NOI,return = F){
    require(igraph)
    master <- connectome
    master_sub <- subset(master, target %in% NOI)
    cells <- unique(master_sub$source)
    ligands <- unique(master_sub$ligand)
    df <- data.frame()
    for (i in 1:length(cells)){
      temp <- subset(master_sub,source == cells[[i]])
      #net <- graph_from_data_frame(temp, directed = T)
      #hub <- hub_score(net,weights = temp$weight_sc,scale = T)$vector
      #auth <- authority_score(net,weights = temp$weight_sc)$vector
      for (j in 1:length(ligands)){
        temp2 <- subset(temp,ligand == ligands[[j]])
        wt.source <- sum(temp2$weight_sc)
        #temp2 <- subset(temp,target == cells[[j]])
        #wt.sink <- sum(temp2$weight_sc)
        row <- data.frame(ligand = ligands[[j]], cells = cells[[i]], wt.source = wt.source,row.names = NULL)
        df <- rbind(df,row)
      }
    }

    # Plots
    p1 <- ggplot(df,aes(ligand,wt.source,color = reorder(cells)))+
      geom_point(size = df$wt.source*1,alpha = 0.6)+
      coord_flip()

      print(p1)
}
