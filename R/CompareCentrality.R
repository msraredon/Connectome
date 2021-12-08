#' CompareCentrality
#'
#' Takes a list of connectomes and will compare sending- and receiving- centrality, side-by-side, for the given networks, or requested network subsets.
#'
#' @param connectome.list A named list of connectomic edgelists. The output plot will be split by the names of the list.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc' or 'weight_norm'. Defaults to 'weight_norm'
#' @param min.z Minimum z-score for ligand and receptor.
#' @param ... Arguments passed to FilterConnectome.  Will be applied to each connectome within the list.
#' @param cols.use Desired colors for cell types, alphabetized. Defaults to standard ggplot colors.
#' @param normalize Default TRUE. Scales each mode to have equivalent x-axes.

#' @export


CompareCentrality <- function (connectome.list,
                              weight.attribute = 'weight_norm',
                              min.z = NULL,
                              cols.use = NULL,
                              normalize = T,...){
require(ggplot2)
require(dplyr)
require(igraph)
require(ggthemes)

if(weight.attribute == 'weight_sc' & is.null(min.z)){
  message("\nWeight attribute is 'weight_sc', recommend also setting min.z = 0 to avoid negative ligand and receptor scores")
}

for(i in 1:length(connectome.list)){
  connectome.list[[i]]$name <- names(connectome.list)[[i]]
}

for(i in 1:length(connectome.list)){
  connectome.list[[i]] <- FilterConnectome(connectome.list[[i]],min.z = min.z,remove.na = T,...)
}


  nn_total <- data.frame()
  for (n in 1:length(connectome.list)) {
    master_func <- connectome.list[[n]]
    master_func$wt <- master_func[,weight.attribute]
    #modes <- as.character(unique(master_func$mode))
    cells <- as.character(unique(union(master_func$source, master_func$target)))

    df <- data.frame()
    #for (i in 1:length(modes)) {
      #temp <- subset(master_func, mode == modes[[i]])
      temp <- master_func
      net <- igraph::graph_from_data_frame(temp, directed = T)
      hub <- hub_score(net, weights = temp$wt, scale = T)$vector
      auth <- authority_score(net, weights = temp$wt,, scale = T)$vector
      total.edgeweight <- sum(temp[,weight.attribute])
      
      for (j in 1:length(cells)) {
        temp2 <- subset(temp, source == cells[[j]])
        wt.source <- sum(temp2$wt)
        temp2 <- subset(temp, target == cells[[j]])
        wt.sink <- sum(temp2$wt)
        if (normalize == T){
          wt.source <- wt.source/total.edgeweight
          wt.sink <- wt.sink/total.edgeweight
        }
        row <- data.frame(#mode = modes[[i]],
                          cells = cells[[j]],
                          hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]],
                          wt.source = wt.source, wt.sink = wt.sink, row.names = NULL,
                          name = master_func$name[1])
        df <- rbind(df, row)
      }
    #}
    nn_total <- rbind(nn_total, df)
  }

  p1 <- ggplot(nn_total, aes(name, wt.source, color = as.factor(cells))) +
    geom_point(size = nn_total$hub.score * 10, alpha = 0.6) +
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    geom_text(data=nn_total %>% group_by(name) %>% top_n(1,hub.score),aes(name,wt.source,label=cells))+
    theme_hc()+
    ylim(0,max(nn_total$wt.source)*1.1)+
    ggtitle('Outgoing Centrality')+
    ylab('Outgoing Edgeweight by Cell Type')+
    xlab('System')+
    theme(plot.title = element_text(hjust = 0.5,face = 'bold'),
          axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
  if (normalize == T){
    p1 <- p1+
      ylab('Outgoing Edgeweight Fraction by Cell Type')
  }
  p2 <- ggplot(nn_total, aes(name,wt.sink, color = as.factor(cells))) +
    geom_point(size = nn_total$auth.score * 10, alpha = 0.6) +
    guides(colour = guide_legend(override.aes = list(size = 10)))+
    geom_text(data=nn_total %>% group_by(name) %>% top_n(1,auth.score),aes(name,wt.sink,label=cells))+
    theme_hc()+
    ylim(0,max(nn_total$wt.sink)*1.1)+
    ggtitle('Incoming Centrality')+
    ylab('Incoming Edgeweight by Cell Type')+
    xlab('System')+
    theme(plot.title = element_text(hjust = 0.5,face = 'bold'),
          axis.text.x = element_text(angle = 45,vjust = 1, hjust=1))
  if (normalize == T){
    p2 <- p2+
      ylab('Incoming Edgeweight Fraction by Cell Type')
  }
  
  # Modify colors if desired
  if (!is.null(cols.use)){
    p1 <- p1 + scale_colour_manual(values = cols.use)
    p2 <- p2 + scale_colour_manual(values = cols.use)
  }
  print(plot_grid(p1, p2))
}

