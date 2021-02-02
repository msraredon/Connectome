#' CircosPlot
#'
#' Plotting function to make Circos plots using the circlize package, following the vignette by the Saeys Lab at: https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md Note that this plotting type is incompatible with edges where the ligand and the receptor are the exact same gene.
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param weight.attribute Column to use to define edgeweights for network analysis. 'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'. If 'weight_sc', function will automatically filter at min.z = 0 to remove negative source/sink values.
#' @param cols.use Optional. Colors for plotting nodes.
#' @param min.z Minimum z-score for ligand and receptor.
#' @param lab.cex Text size for gene names
#' @param balanced.edges Edges in this plot can change thickness along their length. This parameter decides whether to scale edges by a single edgeweight (chosen in weight.attribute) or by the separate cell-specific ligand and receptor values.  Default balanced (TRUE).  If FALSE, the edges will expand or contract to join ligand weight to receptor weight.
#' @param edge.color.by.source Default TRUE - edges will be colored by their source cell type. If false, edges will be colored by receiving cell instead.
#' @param small.gap Default 1. Amount of distance between sectors.  If the number of edges is very large, this will have to be reduced in size.
#' @param big.gap Default 10. Amount of distance between the source cells and the target cells (top and bottom arc of graph).  If the number of edges is very large, this can be reduced in size in addition to 'small.gap'
#' @param title Character string for title of plot.
#' @param ... Arguments passed to FilterConnectome
#' @export

CircosPlot <- function(connectome,
                      weight.attribute = 'weight_sc',
                      cols.use = NULL,
                      min.z = NULL,
                      lab.cex = 1,
                      balanced.edges = T,
                      edge.color.by.source = T,
                      small.gap = 1,
                      big.gap = 10,
                      title = NULL,...){
  library(tidyverse)
  library(circlize)
  library(dplyr)
  library(scales)
  library(ComplexHeatmap)

  # If (weight.attribute != 'weight_norm'){
    if (weight.attribute == 'weight_sc' & is.null(min.z)){
      connectome <- FilterConnectome(connectome, remove.na = T,min.z = 0,...)
    }else{
    connectome <- FilterConnectome(connectome,remove.na = T,min.z = min.z,...)
    }
  #}
  # Pull the dataframe of interest for plotting and format with weight as third column
  connectome$lig.stash <- as.character(connectome$ligand)
  connectome$rec.stash <- as.character(connectome$receptor)
  df <- data.frame(connectome %>% select(ligand,receptor))
  df$ligand <- make.unique(as.character(df$ligand))
  df$receptor <- make.unique(as.character(df$receptor))
  #df$weight <- connectome[,weight.attribute]
  temp <- connectome[,!colnames(connectome) %in% colnames(df)]
  df <- cbind(df,temp)

  # Squash ligands back together to single name if they are duplicates (different edges on same cell type)
  for (i in 1:length(unique(df$lig.stash))){
    temp <- subset(df,lig.stash == unique(df$lig.stash)[i])
    for (j in 1:length(unique(temp$source))){
      temp2 <- subset(temp,source == unique(temp$source)[j])
      dummy <- paste(rep(' ',j-1),collapse = '') # Add number of spaces corresponding to number of unique sources
      df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
    }
    #if(length(unique(temp$source)) == 1){
    #  df[rownames(temp),]$ligand <- as.character(temp$lig.stash)
    #}
  }

  # Squash receptors back together to single name if they are duplicates (different edges on same cell type)
  for (i in 1:length(unique(df$rec.stash))){
    temp <- subset(df,rec.stash == unique(df$rec.stash)[i])
    for (j in 1:length(unique(temp$target))){
      temp2 <- subset(temp,target == unique(temp$target)[j])
      dummy <- paste(rep(' ',j-1),collapse = '') # Add number of spaces corresponding to number of unique targets
      df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
    }
    #if(length(unique(temp$target)) == 1){
    #  df[rownames(temp),]$receptor <- as.character(temp$rec.stash)
    #}
  }

  # Squash ligands back together, by cell type, if they are expressed on multiple cell types
  #temp <- subset(df,ligand != lig.stash) # this is a problem
  #if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$source))){
  #    temp2 <- subset(temp,source == unique(temp$source)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique sources
  #    df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
  #  }
  #}
  # Squash receptors back together, by cell type, if they are expressed on multiple cell types
  #temp <- subset(df,receptor != rec.stash) # this is a problem
  #if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$target))){
  #    temp2 <- subset(temp,target == unique(temp$target)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique targets
  #    df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
  #  }
  #}

  #Establish ordering (order) so that genes are grouped nicely by celltype
  source.order <- df[order(df$source), ]
  target.order <- df[order(df$target), ]
  source.order.un <- unique(source.order[,c('ligand','source')])
  target.order.un <- unique(target.order[,c('receptor','target')])

  source.order$id <- 1:nrow(source.order)
  target.order$id <- 1:nrow(target.order)
  source.order.un$id <- 1:nrow(source.order.un)
  target.order.un$id <- 1:nrow(target.order.un)

  sector.order.un <- c(as.character(source.order.un$ligand),
                       as.character(target.order.un$receptor))

  # Coloring setup
  if (is.null(cols.use)){
    nodes <- as.character(unique(union(df$source,df$target)))
    cols.use <- hue_pal()(length(nodes))
    names(cols.use) <- nodes
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }else{
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }


  # Map to get ligand colorings (edges)
  map <- base::merge(source.order, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.edges <- as.character(map$cols.use)
  names(lig.cols.edges) <- map$ligand

  # Map to get receptor colorings (edges) # this does not work
  map <- base::merge(target.order, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.edges <- as.character(map$cols.use)
  names(rec.cols.edges) <- map$receptor

  # Map to get ligand colorings (sectors)
  map <- base::merge(source.order.un, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.sect <- as.character(map$cols.use)
  names(lig.cols.sect) <- map$ligand

  # Map to get receptor colorings (sectors)
  map <- base::merge(target.order.un, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.sect <- as.character(map$cols.use)
  names(rec.cols.sect) <- map$receptor

  # Make sector colors (grid.cols)
  sectors <- c(source.order.un$ligand,target.order.un$receptor)
  sector.cols <- c(as.character(lig.cols.sect),as.character(rec.cols.sect))

  # Plotting
  # Decide edge order and edge color order
  if (edge.color.by.source == T){
    edge.color <- lig.cols.edges
    df.plot <- source.order
  }else{
    edge.color <- rec.cols.edges
    df.plot <- target.order
  }
  # Decide weight attributes and balanced vs. not
  if (weight.attribute == 'weight_norm'){
    if (balanced.edges == T){
      df.plot <- df.plot[,c('ligand','receptor','weight_norm')]
    }else{
      df.plot <- df.plot[,c('ligand','receptor','ligand.expression','recept.expression')]
    }
  }
  if (weight.attribute == 'weight_sc'){
    if (balanced.edges == T){
      df.plot <- df.plot[,c('ligand','receptor','weight_sc')]
    }else{
      df.plot <- df.plot[,c('ligand','receptor','ligand.scale','recept.scale')]
    }
  }
  #if (weight.attribute == 'score'){
  #  if (balanced.edges == T){
  #    df.plot <- df.plot[,c('ligand','receptor','score')]
  #  }else{
  #    df.plot <- df.plot[,c('ligand','receptor','ligand.norm.lfc','recept.norm.lfc')]
  #  }
  #}

  circos.clear()
  #circos.par(gap.degree = gap.degree)
  chordDiagram(df.plot,
              order = sector.order.un,
              col = edge.color,
              grid.col = sector.cols,
              directional = 1,
              direction.type = "arrows",
              link.arr.type = "big.arrow",
              annotationTrack = "grid",
              preAllocateTracks = 1,
              small.gap = small.gap,
              big.gap = big.gap)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),cex = lab.cex)
    #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  # Make and add legend
  legend <- Legend(at = as.character(unique(union(df$source,df$target))),
                  type = "grid",
                  legend_gp = gpar(fill = as.character(cols.use[as.character(unique(union(df$source,df$target))),]$cols.use)),
                  title_position = "topleft",
                  title = "Cell Type")
  draw(legend, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
  if(!is.null(title)){title(title)}

  p1.base <- recordPlot()
  return(p1.base)
}
