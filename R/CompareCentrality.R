#' CompareCentrality
#'
#' @param connectome.list A NAMED list of connectomic edgelists. The output plot will be split by the names of the list.
#' @param features Genes of interest to generate the network for centrality analysis.
#' @param modes.include Modes to include in centrality analysis.
#' @param nodes.include Cell types to include in the network analysis.
#' @export


CompareCentrality <- function (connectome.list,
                                  features = NULL,
                                  modes.include = NULL,
                                  nodes.include = NULL,
                                  weight.attribute = 'weight_sc'){
require(ggplot2)
require(dplyr)
require(igraph)

if (!is.null(features) & !is.null(modes.include)){message("\nFeatures and Mode filtrations have both been requested. Output will reflect the combination of these selection choices.")}

  nn_total <- data.frame()
  for (n in 1:length(connectome.list)) {
    master_func <- connectome.list[[n]]
    master_func$wt <- master_func$weight_sc
    master_sub <- subset(master_func, source %in% NOI & target %in%
                           NOI)
    modes <- unique(master_sub$mode)
    cells <- unique(union(master_sub$source, master_sub$target))
    df <- data.frame()
    for (i in 1:length(modes)) {
      temp <- subset(master_sub, mode == modes[[i]])
      net <- graph_from_data_frame(temp, directed = T)
      hub <- hub_score(net, weights = temp$wt, scale = T)$vector
      auth <- authority_score(net, weights = temp$wt)$vector
      for (j in 1:length(cells)) {
        temp2 <- subset(temp, source == cells[[j]])
        wt.source <- sum(temp2$wt)
        temp2 <- subset(temp, target == cells[[j]])
        wt.sink <- sum(temp2$wt)
        row <- data.frame(mode = modes[[i]], cells = cells[[j]],
                          hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]],
                          wt.source = wt.source, wt.sink = wt.sink, row.names = NULL,
                          species = master_sub$species[1])
        df <- rbind(df, row)
      }
    }
    nn_total <- rbind(nn_total, df)
  }

  p1 <- ggplot(subset(nn_total, mode == mode_plot), aes(species,
                                                        wt.source, color = reorder(cells))) + geom_point(size = subset(nn_total,
                                                                                                                       mode == mode_plot)$hub.score * 10, alpha = 0.6) +
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    scale_y_continuous(trans='log10')+
    geom_text(data=subset(nn_total, mode == mode_plot) %>% group_by(species) %>% top_n(1,hub.score),aes(species,wt.source,label=cells))

  p2 <- ggplot(subset(nn_total, mode == mode_plot), aes(species,
                                                        wt.sink, color = reorder(cells))) + geom_point(size = subset(nn_total,
                                                                                                                     mode == mode_plot)$auth.score * 10, alpha = 0.6) +
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    scale_y_continuous(trans='log10')+
    geom_text(data=subset(nn_total, mode == mode_plot) %>% group_by(species) %>% top_n(1,auth.score),aes(species,wt.sink,label=cells))

  print(plot_grid(p1, p2))
}
