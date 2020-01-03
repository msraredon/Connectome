# This function takes a list of affinity matrices and prints them as network plots side by side
#' NetworksCompare
#'
#' This function takes a list of affinity matrices (made with 'Affinity') and plots the networks side-by-side.
#' Output will be a PDF in your working directory.
#' @param aff.list A list of affinity matrices made with Affinity

#' @export

NetworksCompare <- function(aff.list,output = c("full"),...){
  require(igraph)
  require(RColorBrewer)
  require(scales)
  pdf(file = paste("AffinityNetworksCompare",Sys.Date(),".pdf",sep=""),width = 9,height = 2)
  par(mfrow=c(1,length(aff.list)))
  par(mar = c(0,0,1,0))
  par(cex.main = 1)
  for (j in 1:length(aff.list)){
    aff <- aff.list[[j]]
    nodes <- colnames(aff)
    if ("full" %in% output){
      net <- graph_from_adjacency_matrix(aff, mode = "directed", diag = FALSE,
      add.colnames = NULL, add.rownames = NA,weighted = TRUE)
      V(net)$color <- hue_pal()(length(nodes))
      V(net)$size <- hub_score(net)$vector*30
      require(RColorBrewer)
      rbPal <- colorRampPalette(c('gray','blue2'))
      E(net)$color <- rbPal(100)[as.numeric(cut(E(net)$weight,breaks = 100))]
      E(net)$width <- E(net)$weight/50
      E(net)$arrow.size <- 0#log(E(net)$weight)/8
      lay.use <- layout_in_circle(net)
      V(net)$label.dist <- 0 #2.5

          radian.rescale <- function(x, start=0, direction=1) {
            c.rotate <- function(x) (x + start) %% (2 * pi) * direction
            c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
          }

      V(net)$label.degree <- radian.rescale(x=1:length(nodes), direction=-1, start=0)
      V(net)$label.family <- "Helvetica"
      V(net)$label.cex <- 0.6
      V(net)$frame.color <- NA
      plot(net, layout=lay.use,
        main=paste("Full Network (Hub) - ",names(aff.list)[j]),
        vertex.label.color = "black")
    }
  }
  dev.off()
}
