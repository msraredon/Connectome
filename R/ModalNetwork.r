# This requires the weight column to be called weight_sc.  Should be revised to be more flexible
#' ModalNetwork
#'
#' This function takes a connectomic edgelist and creates a visualization of the cell-cell network for a given mode. Must have modal annotations for this to work.
#' @param connectome A connectomic edgelist
#' @param MOI The signaling mode of interest to include in the network analysis and subsequent plotting
#' @export

ModalNetwork <- function(connectome,MOI){
  require(igraph)
  require(RColorBrewer)
  require(scales)
  nodes <- sort(union(droplevels(connectome$source),droplevels(connectome$target)))
  edgelist <- subset(connectome,connectome$mode == MOI)

      net <- graph_from_data_frame(d = edgelist, vertices = nodes, directed = T)
      lay <- layout_in_circle(net)
      V(net)$label.family <- "Helvetica"
      V(net)$color <- hue_pal()(length(nodes))
      V(net)$size <- 20
      V(net)$label.cex <- 0.5
      E(net)$weight <- get.edge.attribute(net, "weight_sc")
      E(net)$color <- "grey"
      rbPal <- colorRampPalette(c('gray92','black'))
      try(E(net)$color <- rbPal(10)[as.numeric(cut(E(net)$weight,breaks = 10))])
      E(net)$arrow.size <- 1.5
      E(net)$width <- 4
      plot(net, layout=lay,
        edge.label=E(net)$pair, edge.label.family="Helvetica", edge.label.cex=0.4,
        edge.label.color = "black",
        main=paste('Species Conserved',MOI,"Signaling"),
        vertex.label.color = "black")
}
