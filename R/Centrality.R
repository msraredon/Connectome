#' Centrality
#'
#' This function takes a connectomic edgelist and creates a source and a sink mode- and cell- organized dot plot. The y-axis is the discrete variable 'mode',
#' and the x-axis is the sum of the weights of all edges for each mode made by each cell. Points are organized by cell type, with the size of the point
#' correlating to the Kleinberg hub score (for source graph) and Kleinberg authority score (for sink graph). Network filtration is performed prior to network centrality calculations.
#'
#' @param connectome A connectomic edgelist
#' @param cols.use Desired colors for cell types, alphabetized. Defaults to standard ggplot colors.
#' @param weight.attribute Column to use to define edgeweights for network analysis.  'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'
#' @param min.z Minimum z-score for ligand and receptor.
#' @param normalize Default TRUE. Scales each mode to have equivalent x-axes.
#' @param group.by Default 'mode'. Determines how to subdivide network for centrality analysis.  Accepts 'mode','ligand','receptor', 'gene', or 'mechanism'.
#' @param ... Arguments passed to FilterConnectome


#' @export

Centrality <- function(connectome,
                        cols.use = NULL,
                        weight.attribute = 'weight_sc',
                        min.z = NULL,
                        normalize = T,
                        group.by = 'mode',...){
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
    ligands <- as.character(unique(master_sub$ligand))
    recepts <- as.character(unique(master_sub$receptor))
    genes <- as.character(union(unique(master_sub$ligand),unique(master_sub$receptor)))
    pairs <- as.character(unique(master_sub$pair))

    cells <- as.character(unique(union(master_sub$source,master_sub$target)))

    # Determine which type of network plot:
    if (group.by == 'mode'){
      groups <- modes
    }
    if (group.by == 'ligand'){
      groups <- ligands
    }
    if (group.by == 'receptor'){
      groups <- recepts
    }
    if (group.by == 'gene'){
      groups <- genes
    }
    if (group.by == 'mechanism'){
      groups <- pairs
    }

    # Process data
    df <- data.frame()

    for (i in 1:length(groups)){
      if (group.by == 'mode'){
      temp <- subset(master_sub,mode == groups[[i]])
      }
      if (group.by == 'ligand'){
      temp <- subset(master_sub,ligand == groups[[i]])
      }
      if (group.by == 'receptor'){
      temp <- subset(master_sub,receptor == groups[[i]])
      }
      if (group.by == 'gene'){
      temp <- subset(master_sub,ligand == groups[[i]] | receptor == groups[[i]])
      }
      if (group.by == 'mechanism'){
      temp <- subset(master_sub,pair == groups[[i]])
      }

      # Make network for analysis
      net <- igraph::graph_from_data_frame(temp, directed = T)
      hub <- hub_score(net,weights = temp[,weight.attribute], scale = T)$vector
      auth <- authority_score(net,weights = temp[,weight.attribute], scale = T)$vector
      total.edgeweight <- sum(temp[,weight.attribute])
      # Get cumulative source/sink edgeweights per cell type
      for (j in 1:length(cells)){
        temp2 <- subset(temp,source == cells[[j]])
        wt.source <- sum(temp2[,weight.attribute])
        temp2 <- subset(temp,target == cells[[j]])
        wt.sink <- sum(temp2[,weight.attribute])
        if (normalize == T){
          wt.source <- wt.source/total.edgeweight
          wt.sink <- wt.sink/total.edgeweight
        }
      # Compile info for single network
        row <- data.frame(
          group = groups[[i]],
          cells = cells[[j]],
          hub.score = hub[cells[[j]]],
          auth.score = auth[cells[[j]]],
          wt.source = wt.source,
          wt.sink = wt.sink,
          row.names = NULL)

        df <- rbind(df,row)
      }
    }
    
    #Alphabetize
    df$group <- factor(df$group,levels = sort(as.character(unique(df$group)),decreasing = T))
    
    # Plots
    p1 <- ggplot(df,aes(x=group,y=wt.source,color = as.factor(cells)))+
      geom_point(size = df$hub.score*10,alpha = 0.6)+
      coord_flip()+
      theme(legend.position="none") + ggtitle('Outgoing Centrality')+
      geom_text(data=df %>% group_by(group) %>% top_n(1,hub.score),aes(group,wt.source,label=cells))+
      ylab('Outgoing Edgeweight by Cell Type')
      if (normalize == T){
        p1 <- p1+
        ylab('Outgoing Edgeweight Fraction by Cell Type')+
          ylim(0,1)
      }

    p2 <- ggplot(df,aes(group,wt.sink,color = as.factor(cells)))+
      geom_point(size = df$auth.score*10,alpha = 0.6)+
      coord_flip()+
      theme(legend.position="none") + ggtitle('Incoming Centrality')+
      geom_text(data=df %>% group_by(group) %>% top_n(1,auth.score),aes(group,wt.sink,label=cells))+
      ylab('Incoming Edgeweight by Cell Type')
    if (normalize == T){
      p2 <- p2+
      ylab('Incoming Edgeweight Fraction by Cell Type')+
        ylim(0,1)
    }

    # Define vertical axis label
    if (group.by == 'mode'){
      p1 <- p1 + xlab('Network (by Family)')
      p2 <- p2 + xlab('Network (by Family)')
    }
    if (group.by == 'ligand'){
      p1 <- p1 + xlab('Network (by Ligand)')
      p2 <- p2 + xlab('Network (by Ligand)')
    }
    if (group.by == 'receptor'){
      p1 <- p1 + xlab('Network (by Receptor)')
      p2 <- p2 + xlab('Network (by Receptor)')
    }
    if (group.by == 'gene'){
      p1 <- p1 + xlab('Network (by Gene)')
      p2 <- p2 + xlab('Network (by Gene)')
    }
    if (group.by == 'mechanism'){
      p1 <- p1 + xlab('Network (by Mechanism)')
      p2 <- p2 + xlab('Network (by Mechanism)')
    }

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
