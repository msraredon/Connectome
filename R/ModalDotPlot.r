#' ModalDotPlot
#'
#' This function takes a connectomic edgelist and creates a source and a sink mode- and cell- organized dot plot. The y-axis is the discrete variable 'mode',
#' and the x-axis is the sum of the weights of all edges for each mode made by each cell. Points are organized by cell type, with the size of the point
#' correlating to the Kleinberg hub score (for source graph) and Kleinberg authority score (for sink graph). Network filtration is performed prior to network centrality calculations.
#'
#' @param connectome A connectomic edgelist
#' @param cols.use Desired colors for cell types, alphabetized. Defaults to standard ggplot colors.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'
#' @param ... Arguments passed to FilterConnectome
#' @param min.z Minimum z-score for ligand and receptor.

#' @export

ModalDotPlot <- function(connectome,
                        cols.use = NULL,
                        weight.attribute = 'weight_sc',
                        min.z = NULL,...){
    require(igraph)
    require(ggplot2)
    require(cowplot)
    require(dplyr)

    if(weight.attribute == 'weight_sc' & is.null(min.z)){
      message("\nWeight attribute is 'weight_sc', recommend also setting min.z = 0 to avoid negative ligand and receptor scores")
    }

    # Store
    master <- connectome

    # Filter as demanded (remove NAs at minimum)
    master_sub <- FilterConnectome(master,min.z = min.z,remove.na = T,...)

    # Set up to plot ModalDotPlot

    modes <- as.character(unique(master_sub$mode))
    cells <- as.character(unique(union(levels(master_sub$source),levels(master_sub$target))))

    df <- data.frame()
    for (i in 1:length(modes)){
      temp <- subset(master_sub,mode == modes[[i]])
      net <- graph_from_data_frame(temp, directed = T)
      hub <- hub_score(net,weights = temp[,weight.attribute], scale = T)$vector
      auth <- authority_score(net,weights = temp[,weight.attribute], scale = T)$vector
      for (j in 1:length(cells)){
        temp2 <- subset(temp,source == cells[[j]])
        wt.source <- sum(temp2[,weight.attribute])
        temp2 <- subset(temp,target == cells[[j]])
        wt.sink <- sum(temp2[,weight.attribute])
        row <- data.frame(mode = modes[[i]], cells = cells[[j]], hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]], wt.source = wt.source, wt.sink = wt.sink,row.names = NULL)
        df <- rbind(df,row)
      }
    }
    # Plots
    p1 <- ggplot(df,aes(mode,wt.source,color = reorder(cells)))+
      geom_point(size = df$hub.score*10,alpha = 0.6)+
      coord_flip()+ theme(legend.position="none") + ggtitle('Outgoing edgeweight')+
      geom_text(data=df %>% group_by(mode) %>% top_n(1,hub.score),aes(mode,wt.source,label=cells))


    p2 <- ggplot(df,aes(mode,wt.sink,color = reorder(cells)))+
      geom_point(size = df$auth.score*10,alpha = 0.6)+
      coord_flip()+ theme(legend.position="none") + ggtitle('Incoming edgeweight')+
      geom_text(data=df %>% group_by(mode) %>% top_n(1,auth.score),aes(mode,wt.sink,label=cells))

    # Modify colors if desired
    if (!is.null(cols.use)){
      p1 <- p1 + scale_colour_manual(values = cols.use)
      p2 <- p2 + scale_colour_manual(values = cols.use)
    }
    # Put legend on bottom
    legend <- get_legend(
        p1 +
          guides(color = guide_legend(nrow = 2,byrow=TRUE,override.aes = list(size=5))) +
          theme(legend.position = "bottom")
      )
    # Assemble plot
    plot.top <- plot_grid(p1, p2,nrow = 1)
    return(plot_grid(plot.top,legend,ncol = 1,rel_heights = c(1, .1)))
}
