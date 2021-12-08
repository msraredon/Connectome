#' NetworkPlot
#'
#' Creates a network plot of a connectomic object.  Wrapper for igraph functionality.
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param weight.attribute Column to use to define edgeweights for network analysis. 'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'
#' @param title Description of the network being plotted
#' @param ... Arguments passed to FilterConnectome
#' @param cols.use Optional. Colors for plotting nodes.
#' @param min.z Minimum z-score for ligand and receptor.
#' @param mar Default 1. Symmetric margin around plot.

#' @export

NetworkPlot <- function(connectome,
                        weight.attribute = 'weight_sc',
                        title = NULL,
                        cols.use = NULL,
                        include.all.nodes = F,
                        min.z = NULL,
                        mar = 1.5,
                        layout = 'circle',
                        edge.label.cex = 0.4,...){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)
  require(scales)

  # Record total nodes
  nodes <- as.character(sort(unique(union(connectome$source, connectome$target))))

  # Perform filtration
  if (weight.attribute == 'weight_sc' & is.null(min.z)){
    connectome <- FilterConnectome(connectome, remove.na = T,min.z = 0,...)
  }else{
  connectome <- FilterConnectome(connectome,remove.na = T,min.z = min.z,...)
  }

  # Define nodes for plot
  if(!include.all.nodes){
    nodes <- as.character(sort(unique(union(connectome$source, connectome$target))))
  }

  # igraph based plotting
  
    # Set margins
    par(mar = c(mar, mar, mar, mar),
        xpd = NA,
        bg = "transparent")
    # Prep data
    edgelist <- connectome
    net <- igraph::graph_from_data_frame(d = edgelist, vertices = nodes, directed = T)
    lay <- layout_in_circle(net)
    
    if (layout == 'force.directed'){
      lay <- layout_with_fr(net)
    }

    # Set node colors
      if (!is.null(cols.use)){
        V(net)$color <- cols.use
      }else{V(net)$color <- scales::hue_pal()(length(nodes))}

      V(net)$size <- 20
      #V(net)$size <- hub_score(net)$vector*30
      V(net)$frame.color <- NA
      V(net)$label.dist <- 0 #2.5

                radian.rescale <- function(x, start=0, direction=1) {
                  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
                  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
                }

      V(net)$label.degree <- radian.rescale(x=1:length(nodes), direction=-1, start=0)
      V(net)$label.family <- "Helvetica"
      V(net)$label.cex <- 0.6

      E(net)$weight <- get.edge.attribute(net, weight.attribute)
      E(net)$color <- "grey"
      rbPal <- colorRampPalette(c('gray92','black'))
      try(E(net)$color <- rbPal(10)[as.numeric(cut(E(net)$weight,breaks = 10))])
      E(net)$arrow.size <- 0.5
      E(net)$width <- E(net)$weight
    if (!is.null(title)){
      plot(net, layout=lay,
        edge.label=E(net)$pair, edge.label.family="Helvetica", edge.label.cex=edge.label.cex,
        edge.label.color = "black",
        main=title,
        vertex.label.color = "black")
    }else{
      plot(net, layout=lay,
        edge.label=E(net)$pair, edge.label.family="Helvetica", edge.label.cex=edge.label.cex,
        edge.label.color = "black",
        vertex.label.color = "black")
    }
      p1.base <- recordPlot()
      return(p1.base)
}
