# This requires the weight column to be called 'weight_sc'.  Should be revised to be more flexible
#' ModalDotPlot
#'
#' This function takes a connectomic edgelist and creates a source and a sink mode- and cell- organized dot plot. The y-axis is binned by mode,
#' and the x-axis is the sum of the weight of all edges for each mode made by each cell. Points are binned by cell type, with the size of the point
#' correlating to the Kleinberg hub score (for source graph) and Kleinberg authority score (for sink graph)
#' @param connectome A connectomic edgelist
#' @param NOI The nodes of interest to include in the network analysis and subsequent plotting
#' @export

ModalDotPlot <- function(connectome,NOI,return = F){
    require(igraph)
    master <- connectome
    master_sub <- subset(master, source %in% NOI)
    master_sub <- subset(master_sub, target %in% NOI)
    modes <- unique(master_sub$mode)
    cells <- unique(union(master_sub$source,master_sub$target))
    df <- data.frame()
    for (i in 1:length(modes)){
      temp <- subset(master_sub,mode == modes[[i]])
      net <- graph_from_data_frame(temp, directed = T)
      hub <- hub_score(net,weights = temp$weight_sc, scale = T)$vector
      auth <- authority_score(net,weights = temp$weight_sc, scale = T)$vector
      for (j in 1:length(cells)){
        temp2 <- subset(temp,source == cells[[j]])
        wt.source <- sum(temp2$weight_sc)
        temp2 <- subset(temp,target == cells[[j]])
        wt.sink <- sum(temp2$weight_sc)
        row <- data.frame(mode = modes[[i]], cells = cells[[j]], hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]], wt.source = wt.source, wt.sink = wt.sink,row.names = NULL)
        df <- rbind(df,row)
      }
    }
    # Plots
    p1 <- ggplot(df,aes(mode,wt.source,color = reorder(cells)))+
      geom_point(size = df$hub.score*10,alpha = 0.6)+
      coord_flip()+
      guides(colour = guide_legend(override.aes = list(size=10)))

    p2 <- ggplot(df,aes(mode,wt.sink,color = reorder(cells)))+
      geom_point(size = df$auth.score*10,alpha = 0.6)+
      coord_flip()+
      guides(colour = guide_legend(override.aes = list(size=10)))

      plot_grid(p1,p2)
}
