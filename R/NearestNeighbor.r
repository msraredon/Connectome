# requires weight column to be called 'weight_sc' <-- should be revised to be more flexible
#' NearestNeighbor
#'
#' This function takes a connectome and plots a visualization of cell-cell connectivity showing the sum of the weights of all edges between each cell combo.
#'
#' @param connectome A connectomic edgelist
#' @param NOI The nodes of interest to include
#' @export

NearestNeighbor <- function(connectome,NOI,return = F){
    require(igraph)
    master <- connectome
    master_sub <- subset(master, source %in% NOI)
    master_sub <- subset(master_sub, target %in% NOI)
    cells <- unique(union(master_sub$source,master_sub$target))
    nn <- data.frame()
    for (i in 1:length(cells)){
      temp <- subset(master_sub,target == cells[[i]])
      net <- graph_from_data_frame(temp, directed = T)
      hub <- hub_score(net,weights = temp$weight_sc,scale = T)$vector
      auth <- authority_score(net,weights = temp$weight_sc)$vector
      for (j in 1:length(cells)){
        temp2 <- subset(temp,source == cells[[j]])
        wt.source <- sum(temp2$weight_sc)
        temp2 <- subset(temp,target == cells[[j]])
        wt.sink <- sum(temp2$weight_sc)
        row <- data.frame(sink = cells[[i]], source = cells[[j]], hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]], wt.source = wt.source, wt.sink = wt.sink,row.names = NULL)
        nn <- rbind(nn,row)
      }
    }
    # Dot plot
    p <- ggplot(nn,aes(sink,wt.source,color = reorder(source)))+
      geom_point(size = nn$hub.score*10,alpha = 0.6)+
      coord_flip()
    print(p)
    if (return == T){
      return(nn)
    }

}
