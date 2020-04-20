#' DifferentialScoringPlot
#'
#' Currently in beta testing. Creates x3 aligned heatmaps allowing visualization of ligand, receptor, and perturbation scores for a given cell-system of interest.
#'
#' @param differential.connectome A differential connectome, made with DifferentialConnectome. May be filtered as desired prior to plotting.
#' @param sources.include Source nodes of interest. Output will be limited to edges coming from these sources.
#' @param targets.include Target nodes of interest. Output will be limited to edges landing on these targets.
#' @param features Gene of interest. Output will be limited to edges including these specific genes.
#' @param min.score Default NULL. Will limit output to edges with a differential score greater than this value.
#' @param min.pct Default NULL. Threshold to return clusterwise observations for both ligand and receptor. Only needs to be satisfied in connect.1 OR in connect.2.
#' @param verbose Whether to output feedback to user
#' @param infinity.to.max Default FALSE.  If TRUE, will create a pseudo value to replace values of "Inf"

#' @export

DifferentialScoringPlot <- function(differential.connectome,
                                    features = NULL,
                                    sources.include = NULL,
                                    targets.include = NULL,
                                    min.score = NULL,
                                    min.pct = NULL,
                                    verbose = T,
                                    infinity.to.max = T){
  require(ggplot2)
  require(cowplot)
  require(dplyr)

  data <- differential.connectome
  pre.filter <- nrow(data)

  # Setup vector column
  data$vector <- paste(data$source,data$target,sep = ' - ')

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

  # Set 'Inf' values to a maximum score
  if (infinity.to.max == T){
    if (length(data[data$ligand.norm.lfc == 'Inf',]$ligand.norm.lfc) > 0){
    data[data$ligand.norm.lfc == 'Inf',]$ligand.norm.lfc <- max(data[data$ligand.norm.lfc != 'Inf',]$ligand.norm.lfc)*1.01
    }
    if (length(data[data$recept.norm.lfc == 'Inf',]$recept.norm.lfc) > 0){
    data[data$recept.norm.lfc == 'Inf',]$recept.norm.lfc <- max(data[data$recept.norm.lfc != 'Inf',]$recept.norm.lfc)*1.01
    }
    if (length(data[data$score == 'Inf',]$score) > 0){
    data[data$score == 'Inf',]$score <- max(data[data$score != 'Inf',]$score)*1.01
    }

    if (length(data[data$ligand.norm.lfc == '-Inf',]$ligand.norm.lfc) > 0){
    data[data$ligand.norm.lfc == '-Inf',]$ligand.norm.lfc <- min(data[data$ligand.norm.lfc != '-Inf',]$ligand.norm.lfc)*1.01
    }
    if (length(data[data$recept.norm.lfc == '-Inf',]$recept.norm.lfc) > 0){
    data[data$recept.norm.lfc == '-Inf',]$recept.norm.lfc <- min(data[data$recept.norm.lfc != '-Inf',]$recept.norm.lfc)*1.01
    }
    if (length(data[data$score == '-Inf',]$score) > 0){
    data[data$score == '-Inf',]$score <- min(data[data$score != '-Inf',]$score)*1.01
    }
  }  

  p1 <- ggplot(data,aes(x = vector, y = pair)) +
    geom_tile(aes(fill = ligand.norm.lfc )) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low="blue",mid = 'grey',high="red",midpoint = 0) + ggtitle('Ligand Log2 Fold Change')

  p2 <- ggplot(data,aes(x = vector, y = pair)) +
    geom_tile(aes(fill = recept.norm.lfc )) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low="blue",mid = 'grey',high="red",midpoint = 0) + ggtitle('Receptor Log2 Fold Change')

  p3 <- ggplot(data,aes(x = vector, y = pair)) +
    geom_tile(aes(fill = score )) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_gradient2(low="blue",mid = 'grey',high="red",midpoint = 0) + ggtitle('Perturbation Score')

  plot_grid(p1,p2,p3,nrow = 1)
}
