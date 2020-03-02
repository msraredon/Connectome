#' ModalDotPlot
#'
#' This function takes a connectomic edgelist and creates a source and a sink mode- and cell- organized dot plot. The y-axis is the discrete variable 'mode',
#' and the x-axis is the sum of the weight (weight_sc column) of all edges for each mode made by each cell. Points are organized by cell type, with the size of the point
#' correlating to the Kleinberg hub score (for source graph) and Kleinberg authority score (for sink graph)
#'
#' @param connectome A connectomic edgelist
#' @param NOI The nodes (cell identities) of interest to include in the network analysis and subsequent plotting. Defaults to all nodes.
#' @export

ModalDotPlot <- function(connectome,NOI = NULL){
    require(igraph)
    require(ggplot2)
    master <- connectome
    if (!is.null(NOI)){
      if (length(NOI) == 1){
        master_sub <- subset(master, source == NOI & target == NOI)
      }else{
        master_sub <- subset(master, source %in% NOI & target %in% NOI)
      }
    }else{
      master_sub <- master
    }

    modes <- unique(master_sub$mode)
    cells <- as.factor(unique(union(master_sub$source,master_sub$target)))
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
      coord_flip()+ theme(legend.position="none") + ggtitle('Outgoing edgeweight')

    p2 <- ggplot(df,aes(mode,wt.sink,color = reorder(cells)))+
      geom_point(size = df$auth.score*10,alpha = 0.6)+
      coord_flip()+ theme(legend.position="none") + ggtitle('Incoming edgeweight')

    legend <- get_legend(
        p1 +
          guides(color = guide_legend(nrow = 2,byrow=TRUE,override.aes = list(size=5))) +
          theme(legend.position = "bottom")
      )

    plot.top <- plot_grid(p1, p2,nrow = 1)
    return(plot_grid(plot.top,legend,ncol = 1,rel_heights = c(1, .1)))
}
