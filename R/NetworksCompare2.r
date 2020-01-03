#' NetworksCompare2
#'
#' This function takes a list of connectomes and plots the networks side by side. Performs on a mode of interest (MOI)
#' Output will be an image
#' @param con.list A list of connectomes
#' @param MOI The mode of interest to plot
#' @export

NetworksCompare2 <- function(con.list,MOI,weight.attribute = 'weight_sc',...){
  require(igraph)
  require(RColorBrewer)
  require(scales)
  pdf(file = paste(MOI,"Networks Compare.pdf",sep=" "),width = 11,height = 2)
  par(mfrow=c(1,length(con.list)))
  par(mar = c(0,0,1,0))
  par(cex.main = 1)

  for (j in 1:length(con.list)){
    con <- con.list[[j]]
    nodes <- sort(unique(union(con$source,con$target)))
    edgelist <- subset(con,con$mode == MOI)
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

      plot(net, layout=lay,
        edge.label=E(net)$pair, edge.label.family="Helvetica", edge.label.cex=0.4,
        edge.label.color = "black",
        main=paste(MOI,"Network",names(con.list)[j],sep = ' '),
        vertex.label.color = "black")
    }
    dev.off()
}
