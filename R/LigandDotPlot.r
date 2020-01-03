# This requires the weight column to be called 'weight_sc'.  Should be revised to be more flexible
#' LigandDotPlot
#'
#' This function takes a connectomic edgelist and creates a source and a sink LIGAND- and cell- organized dot plot. The y-axis is binned by ligand,
#' and the x-axis is the sum of the weight of all edges for each ligand made by each cell. Points are binned by cell type, with the size of the point
#' correlating to the Kleinberg hub score (for source graph) and Kleinberg authority score (for sink graph)
#' @param connectome A connectomic edgelist
#' @param NOI The nodes of interest to include in the network analysis and subsequent plotting
#' @export

LigandDotPlot <- function(connectome,NOI,return = F){
    require(igraph)
    master <- connectome
    master_sub <- subset(master, source %in% NOI)
    master_sub <- subset(master_sub, target %in% NOI)
    lig.list <- unique(master_sub$ligand)
    cells <- unique(union(master_sub$source,master_sub$target))
    df <- data.frame()
    for (i in 1:length(lig.list)){
      temp <- subset(master_sub,ligand == lig.list[[i]])
      net <- graph_from_data_frame(temp, directed = T)
      #hub <- hub_score(net,weights = temp$weight_sc,scale = T)$vector
      #auth <- authority_score(net,weights = temp$weight_sc)$vector
      hub <- strength(net,mode = 'out',weights = temp$weight_sc)
      auth <- strength(net,mode = 'in',weights = temp$weight_sc)
      for (j in 1:length(cells)){
        temp2 <- subset(temp,source == cells[[j]])
        wt.source <- sum(temp2$weight_sc)
        temp2 <- subset(temp,target == cells[[j]])
        wt.sink <- sum(temp2$weight_sc)
        row <- data.frame(ligand = lig.list[[i]], cells = cells[[j]], hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]], wt.source = wt.source, wt.sink = wt.sink,row.names = NULL)
        df <- rbind(df,row)
      }
    }
    # Plots
    p1 <- ggplot(df,aes(ligand,wt.source,color = reorder(cells)))+
      geom_point(size = df$wt.source/5,alpha = 0.6)+
      coord_flip()

    p2 <- ggplot(df,aes(ligand,wt.sink,color = reorder(cells)))+
      geom_point(size = df$wt.sink/5,alpha = 0.6)+
      coord_flip()

      plot_grid(p1,p2)
}
