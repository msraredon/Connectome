#' NetworkPlot
#'
#' Creates a network 'hairball' plot of a filtered object.  Wrapper for igraph functionality.
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param weight.attribute The desired column to use for edgeweights. Defaults to 'weight_sc'
#' @export

NetworkPlot <- function(connectome, weight.attribute = 'weight_sc',title = NULL,...){

    con <- connectome
    nodes <- sort(unique(union(con$source,con$target)))
    edgelist <- con
    net <- graph_from_data_frame(d = edgelist, vertices = nodes, directed = T)
    lay <- layout_in_circle(net)
      V(net)$color <- hue_pal()(length(nodes))
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
      E(net)$width <- 1
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
