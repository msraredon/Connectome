#' DiffEdgeDotPlot
#'
#' Plotting function to make a DotPlot of differential edges for a network comparison. This function first finds all edges meeting the desired thresholding criteria, and then plots complete information regarding all mechanisms and celltype vectors implicated.
#'
#' @param differential.connectome A differential connectome object
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param features Gene of interest. Output will be limited to edges including these specific genes.
#' @param min.score Default NULL. Will limit output to edges with a differential score greater than this value.
#' @param min.pct Default NULL. Threshold to return clusterwise observations for both ligand and receptor. Only needs to be satisfied in connect.1 OR in connect.2.
#' @param verbose Whether to output feedback to user
#' @param infinity.to.max Default TRUE.  If TRUE, will create a pseudo value to replace values of "Inf"
#' @export

DiffEdgeDotPlot <- function(differential.connectome,
                            features = NULL,
                            sources.include = NULL,
                            targets.include = NULL,
                            min.score = NULL,
                            min.pct = NULL,
                            verbose = T,
                            infinity.to.max = T,...){
  
  library(ggplot2)
  
  # Setup vector column
  differential.connectome$vector <- paste(differential.connectome$source,
                                          differential.connectome$target,
                                          sep = ' - ')
  
  # Set 'Inf' values to a maximum score
  if (infinity.to.max == T){
    if (length(differential.connectome[differential.connectome$ligand.norm.lfc == 'Inf',]$ligand.norm.lfc) > 0){
      differential.connectome[differential.connectome$ligand.norm.lfc == 'Inf',]$ligand.norm.lfc <- max(differential.connectome[differential.connectome$ligand.norm.lfc != 'Inf',]$ligand.norm.lfc)*1.01
    }
    if (length(differential.connectome[differential.connectome$recept.norm.lfc == 'Inf',]$recept.norm.lfc) > 0){
      differential.connectome[differential.connectome$recept.norm.lfc == 'Inf',]$recept.norm.lfc <- max(differential.connectome[differential.connectome$recept.norm.lfc != 'Inf',]$recept.norm.lfc)*1.01
    }
    if (length(differential.connectome[differential.connectome$score == 'Inf',]$score) > 0){
      differential.connectome[differential.connectome$score == 'Inf',]$score <- max(differential.connectome[differential.connectome$score != 'Inf',]$score)*1.01
    }
    
    if (length(differential.connectome[differential.connectome$ligand.norm.lfc == '-Inf',]$ligand.norm.lfc) > 0){
      differential.connectome[differential.connectome$ligand.norm.lfc == '-Inf',]$ligand.norm.lfc <- min(differential.connectome[differential.connectome$ligand.norm.lfc != '-Inf',]$ligand.norm.lfc)*1.01
    }
    if (length(differential.connectome[differential.connectome$recept.norm.lfc == '-Inf',]$recept.norm.lfc) > 0){
      differential.connectome[differential.connectome$recept.norm.lfc == '-Inf',]$recept.norm.lfc <- min(differential.connectome[differential.connectome$recept.norm.lfc != '-Inf',]$recept.norm.lfc)*1.01
    }
    if (length(differential.connectome[differential.connectome$score == '-Inf',]$score) > 0){
      differential.connectome[differential.connectome$score == '-Inf',]$score <- min(differential.connectome[differential.connectome$score != '-Inf',]$score)*1.01
    }
  }
  
  #Filter
  data <- differential.connectome
  pre.filter <- nrow(data)

  # Subset based on min.score
  if (!is.null(min.score)){
    data <- subset(data,score > min.score)
  }
  
  # Subset on nodes (cell types) of interest
  if (!is.null(sources.include)){
    data <- subset(data, source %in% sources.include)
  }
  if (!is.null(targets.include)){
    data <- subset(data, target %in% targets.include)
  }
  # Subset on features of interest
  if (!is.null(features)){
    data <- subset(data,ligand %in% features | receptor %in% features)
  }
  # Subset based on min.pct
  if (!is.null(min.pct)){
    data <- subset(data,pct.source.1 > min.pct | pct.source.2 > min.pct)
    data <- subset(data, pct.target.1 > min.pct | pct.target.2 > min.pct)
  }
  
  # Postfilter
  post.filter <- nrow(data)
  
  if (verbose){
    message(paste("\nPre-filter edges: ",as.character(pre.filter)))
    message(paste("\nPost-filter edges: ",as.character(post.filter)))
    message("\nConnectome filtration completed")
  }
  
  
  # Identify vectors and mechanisms remaining
  vectors <- unique(data$vector)
  mechanisms <- unique(data$pair)
  
  # Subset to these
  data2 <- subset(differential.connectome,pair %in% mechanisms & vector %in% vectors)
  
  #Plot
  
  my_palette <- colorRampPalette(c("blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  
  p1 <- ggplot(data = data2,aes(x = vector, y = pair,
                                #alpha = log(score+1)
  )) +
    geom_point(aes(size=log(score+1),color=log(score+1),stroke = 0))+
    # theme_minimal()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_gradientn('log(score + 1)', colors=my_palette)+theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  return(p1)
}
