#' NetworkPlot
#'
#' Creates a network plot of a connectomic object.  Wrapper for igraph functionality.
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param features Genes of interest. All edges containing these features will be plotted.
#' @param modes.include Modes to be plotted. Defaults to all modes. Can be used to look at a narrow category of signaling.
#' @param weight.attribute Column to use to define edgeweights for network analysis. 'weight_sc','weight_norm', or 'weight_raw'. Defaults to 'weight_sc'
#' @param title Description of the network being plotted
#' @param min.pct Minimum fraction of cells within a given cluster expressing the ligand or receptor. Defaults to 0.10, allows NULL.
#' @export

NetworkPlot <- function(connectome,
                        features = NULL,
                        weight.attribute = 'weight_sc',
                        title = NULL,
                        modes.include = NULL,
                        cols.use = NULL,
                        min.pct = 0.10,...){
  require(igraph)
  require(ggplot2)
  require(cowplot)
  require(dplyr)
  require(scales)
  #Define nodes for plot
    nodes <- as.character(sort(unique(union(connectome$source, connectome$target))))
  # Subset to modes of interest
    if (!is.null(modes.include)){
      if (length(modes.include) == 1){
        connectome <- subset(connectome, mode == modes.include)
      }else{
        connectome <- subset(connectome, mode %in% modes.include)
      }
    }
  # Subset to features of interest
    if (!is.null(features)){
      if (length(features) == 1){
        connectome <- subset(connectome, ligand == features | receptor == features)
      }else{
        connectome <- subset(connectome, ligand %in% features | receptor %in% features)
      }
    }
  # Subset to expression percentage (and/or remove NA values, if applicable)
  if (!is.null(min.pct)){
    connectome <- subset(connectome, percent.source > min.pct & percent.target > min.pct)
  }else{
    connectome <- connectome[!is.na(connectome$percent.source),]
    connectome <- connectome[!is.na(connectome$percent.target),]
  }


  # igraph based plotting
    edgelist <- connectome
    net <- graph_from_data_frame(d = edgelist, vertices = nodes, directed = T)
    lay <- layout_in_circle(net)

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
        edge.label=E(net)$pair, edge.label.family="Helvetica", edge.label.cex=0.4,
        edge.label.color = "black",
        main=paste(title,"Network",sep = ' '),
        vertex.label.color = "black")
    }else{
      plot(net, layout=lay,
        edge.label=E(net)$pair, edge.label.family="Helvetica", edge.label.cex=0.4,
        edge.label.color = "black",
        vertex.label.color = "black")
      }
}
