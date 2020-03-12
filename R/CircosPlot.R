#' CircosPlot
#'
#' Plotting function to make Circos plots using the circlize package, following the vignette by the Saeys Lab at: https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param weight.attribute Column to use to define edgeweights for network analysis. 'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'. If 'weight_sc', function will automatically filter at min.z = 0 to remove negative source/sink values.
#' @param ... Arguments passed to FilterConnectome
#' @param cols.use Optional. Colors for plotting nodes.
#' @param min.z Minimum z-score for ligand and receptor.
#'
#' @export

CircosPlot <- function(connectome,
                      weight.attribute = 'weight_sc',
                      cols.use = NULL,
                      min.z = NULL,...){
  library(tidyverse)
  library(circlize)
  library(dplyr)
  library(scales)

  # Perform filtration
    if (weight.attribute == 'weight_sc' & is.null(min.z)){
      connectome <- FilterConnectome(connectome, remove.na = T,min.z = 0,...)
    }else{
    connectome <- FilterConnectome(connectome,remove.na = T,min.z,...)
    }

  # Pull the dataframe of interest for plotting and format with weight as third column
  df <- data.frame(connectome %>% select(ligand,receptor))
  df$ligand <- make.unique(as.character(df$ligand))
  df$receptor <- make.unique(as.character(df$receptor))
  df$weight <- connectome[,weight.attribute]
  temp <- connectome[,!colnames(connectome) %in% colnames(df)]
  df <- cbind(df,temp)
  #df$id <- 1:nrow(df)

  #Establish ordering
  source.order <- df[order(df$source), ]
  target.order <- df[order(df$target), ]
  source.order$id <- 1:nrow(source.order)
  target.order$id <- 1:nrow(target.order)
  sector.order <- c(as.character(source.order$ligand),as.character(target.order$receptor))

  # Coloring based on cell types (not ligand and/or receptor, which are the sectors in this type of plot)
  if (is.null(cols.use)){
    nodes <- as.character(unique(union(df$source,df$target)))
    cols.use <- hue_pal()(length(nodes))
    names(cols.use) <- nodes
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }

  # Map to get ligand colorings
  map <- base::merge(source.order, cols.use, by.x = "source", by.y = "cell", all = TRUE)
  map <- map[order(map$id), ]
  lig.cols <- as.character(map$cols.use)
  names(lig.cols) <- map$ligand

  # Map to get receptor colorings
  map <- base::merge(target.order, cols.use, by.x = "target", by.y = "cell", all = TRUE)
  map <- map[order(map$id), ]
  rec.cols <- as.character(map$cols.use)
  names(rec.cols) <- map$receptor

  # Make sector colors (grid.cols)
  sectors <- c(source.order$ligand,target.order$receptor)
  sector.cols <- c(as.character(lig.cols),as.character(rec.cols))

  # Plotting
  circos.clear()
  chordDiagram(df[,1:3],
              order = sector.order,
              col = lig.cols,
              grid.col = sector.cols,
              directional = 1,
              direction.type = "arrows",
              link.arr.type = "big.arrow",annotationTrack = "grid",preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)

}
