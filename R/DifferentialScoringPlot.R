#' DifferentialScoringPlot
#'
#' Currently in beta testing. Creates x3 aligned heatmaps allowing visualization of ligand, receptor, and perturbation scores for a given cell-system of interest.
#'
#' @param differential.connectome A differential connectome, made with DifferentialConnectome. May be filtered as desired prior to plotting.
#' @param min.score Default NULL. Threshold to prioritize only strongly perturbed edges.
#' @param sources.include
#' @param targets.include

#' @export

DifferentialScoringPlot <- function(differential.connectome,
                                    sources.include = NULL,
                                    targets.include = NULL,
                                    min.score = NULL){

  data <- differential.connectome

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

  p1 <- ggplot(data,aes(x = vector, y = pair)) +
    geom_tile(aes(fill = ligand.norm.lfc )) +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="blue",mid = 'grey',high="red",midpoint = 0) + ggtitle('Ligand Log2 Fold Change')

  p2 <- ggplot(data,aes(x = vector, y = pair)) +
    geom_tile(aes(fill = recept.norm.lfc )) +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="blue",mid = 'grey',high="red",midpoint = 0) + ggtitle('Receptor Log2 Fold Change')

  p3 <- ggplot(data,aes(x = vector, y = pair)) +
    geom_tile(aes(fill = score )) +
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_gradient2(low="blue",mid = 'grey',high="red",midpoint = 0) + ggtitle('Perturbation Score')

  plot_grid(p1,p2,p3,nrow = 1)
}
