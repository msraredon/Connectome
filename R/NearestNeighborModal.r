#' Nearest Neighbor
#'
#' This function makes a mode-specific dotplot showing nearest neighbors within a networks
#' Can compare a LIST of connectomes
#' requires weight column to be called 'weight_sc' <-- should be revised to be more flexible

#' @param connectome.list A list of connectomic edgelists
#' @param NOI The nodes of interest to include in the network analysis and subsequent plotting
#' @param mode_plot The mode of interest to plot
#' @export

NearestNeighborModal <- function(connectome.list,NOI,mode_plot = 'Growth factors'){

    nn_total <- data.frame()
    for (n in 1:length(connectome.list)){
        master_func  <- connectome.list[[n]]
        master_func$wt <- master_func$weight_sc ####LOOK CLOSELY HERE!!!
        master_sub <- subset(master_func, source %in% NOI & target %in% NOI)

        modes <- unique(master_sub$mode)
        cells <- unique(union(master_sub$source,master_sub$target))

        df <- data.frame()
        for (i in 1:length(modes)){
          temp <- subset(master_sub,mode == modes[[i]])
          net <- graph_from_data_frame(temp, directed = T)
          hub <- hub_score(net,weights = temp$wt,scale = T)$vector
          auth <- authority_score(net,weights = temp$wt)$vector
          for (j in 1:length(cells)){
            temp2 <- subset(temp,source == cells[[j]])
            wt.source <- sum(temp2$wt)
            temp2 <- subset(temp,target == cells[[j]])
            wt.sink <- sum(temp2$wt)
            row <- data.frame(mode = modes[[i]], cells = cells[[j]], hub.score = hub[cells[[j]]], auth.score = auth[cells[[j]]], wt.source = wt.source, wt.sink = wt.sink,row.names = NULL,species = master_sub$species[1])
            df <- rbind(df,row)
          }
        }
        nn_total <- rbind(nn_total,df)
    }

    # Plots
    p1 <- ggplot(subset(nn_total,mode == mode_plot),aes(species,wt.source,color = reorder(cells)))+
      geom_point(size = subset(nn_total,mode == mode_plot)$hub.score*10,alpha = 0.6)+
      guides(colour = guide_legend(override.aes = list(size=10)))
    p2 <- ggplot(subset(nn_total,mode == mode_plot),aes(species,wt.sink,color = reorder(cells)))+
      geom_point(size = subset(nn_total,mode == mode_plot)$auth.score*10,alpha = 0.6)+
      guides(colour = guide_legend(override.aes = list(size=10)))

    plot_grid(p1,p2)
}
