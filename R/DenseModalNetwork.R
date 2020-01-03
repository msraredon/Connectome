# This requires the weight column to be called 'lig.wt' and 'rec.wt'.  Should be revised to be more flexible

#' DenseModalNetwork
#'
#' This function takes a connectomic edgelist and creates a visualization of the cell-cell network for a given mode with edges that are the sum of all weights. Must have modal annotations for this to work.
#' @param connectome A connectomic edgelist
#' @param MOI The signaling mode of interest to include in the network analysis and subsequent plotting
#' @export

DenseModalNetwork <- function(connectome,MOI){
  require(igraph)
  require(RColorBrewer)
  require(scales)
  nodes <- sort(union(droplevels(connectome$source),droplevels(connectome$target)))
  edgelist <- subset(connectome,connectome$mode == MOI)
  #Condense
      aff <- Affinity(edgelist,nodes.include = union(unique(edgelist$source),unique(edgelist$target)))
      net <- graph_from_adjacency_matrix(aff, mode = "directed", diag = FALSE,
      add.colnames = NULL, add.rownames = NA,weighted = TRUE)
      lay <- layout_in_circle(net)
      V(net)$label.family <- "Helvetica"
      V(net)$color <- alpha(hue_pal()(length(nodes)),0.8)
      V(net)$size <- 20
      V(net)$label.cex <- 1.5
      V(net)$frame.color <- NA
      E(net)$weight <- get.edge.attribute(net, "weight")
      E(net)$color <- "grey"
      rbPal <- colorRampPalette(c('gray92','black'))
      try(E(net)$color <- rbPal(10)[as.numeric(cut(E(net)$weight,breaks = 10))])
      E(net)$arrow.size <- 1.5
      E(net)$width <- 4
      plot(net, layout=lay, edge.label.family="Helvetica", edge.label.cex=0.4,
        edge.label.color = "black",
        main=paste('Species Conserved',MOI,"Signaling"),
        vertex.label.color = "black",edge.curved = T)
}
