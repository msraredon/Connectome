
# This function processes connectome data to act as input for http://wodaklab.org/hivegraph/graph
# The input connectome must have already been run through the 'Annotate' function to work (must have modal categorizations)

Hive_Organization <- function(connectome){
  master_sub <- connectome
  # Set up initial axes
  cells <- sort(unique(union(droplevels(master_sub$source),droplevels(master_sub$target))))
  modes <- sort(unique(droplevels(master_sub$mode)))
  # Color by mode
  require(randomcoloR)
  cols <- distinctColorPalette(length(modes),altCol = T, runTsne = F)
  cols <- col2rgb(cols)
  #OR
  pal <- colorRampPalette(brewer.pal(9, "Set1"))
  cols <- pal(length(modes))
  cols <- col2rgb(cols)
  colnames(cols) <- modes
  nodes1 <- data.frame()
    for (i in 1:length(cells)){
      for (j in 1:length(modes)){
            temp <- data.frame(cell = cells[i],mode = modes[j])
            nodes1 <- rbind(nodes1,temp)
      }
    }
  nodes1$id <- paste(nodes1$cell,nodes1$mode)
  nodes2 <- data.frame(id = modes)
  nodes3 <- data.frame(id = cells)
  # Get 'heights' and 'values' of each node
  for (m in 1:length(cells)){
    NOI = cells[m]
    n <- master_sub[master_sub$target == NOI,]
    n$wt <- n$lig.wt + n$rec.wt
    # For Axis 1 nodes, set heights and values:
    nodes1$height <- NA
    for (i in 1:nrow(nodes1)){
      temp <- subset(n,source == as.character(nodes1$cell[i]) & mode == nodes1$mode[i])
      nodes1$height[i] <- sum(temp$wt)
    }
    nodes1$height_tr <- nodes1$height
    nodes1[nodes1$height_tr == 0,]$height_tr <- 0.5 #Set a minimum node height
    nodes1 <- within(nodes1,value <- cumsum(height_tr))
    nodes1$value[1] <- 15
    for (i in 2:nrow(nodes1)){
      nodes1$value[i] <- nodes1$value[i-1] + nodes1$height_tr[i-1] + 5
    }
    nodes1$r <- NA
    nodes1$g <- NA
    nodes1$b <- NA
    for (i in 1:length(modes)){
      nodes1[nodes1$mode == modes[i],]$r <- cols[1,modes[i]]
      nodes1[nodes1$mode == modes[i],]$g <- cols[2,modes[i]]
      nodes1[nodes1$mode == modes[i],]$b <- cols[3,modes[i]]
    }
    n1 <- data.frame(id = nodes1$id, axis = "a", label = nodes1$id, value = nodes1$value, height = nodes1$height_tr, r = nodes1$r, g = nodes1$g, b = nodes1$b, a = 200)

    # For Axis 2 nodes,set heights and values:
    nodes2$height <- tapply(n$wt, n$mode, FUN=sum)[nodes2$id]
    nodes2$height_tr <- nodes2$height
    nodes2[is.na(nodes2$height_tr),]$height_tr <- 1 #Set a minimum node height
    nodes2 <- within(nodes2,value <- cumsum(height_tr))
    nodes2$value[1] <- 5
    for (i in 2:nrow(nodes2)){
      nodes2$value[i] <- nodes2$value[i-1] + nodes2$height_tr[i-1] + 2
    }
    nodes2$r <- NA
    nodes2$g <- NA
    nodes2$b <- NA
    for (i in 1:length(modes)){
      nodes2[nodes2$id == modes[i],]$r <- cols[1,modes[i]]
      nodes2[nodes2$id == modes[i],]$g <- cols[2,modes[i]]
      nodes2[nodes2$id == modes[i],]$b <- cols[3,modes[i]]
    }
    n2 <- data.frame(id = nodes2$id, axis = "b", label = nodes2$id, value = nodes2$value, height = nodes2$height_tr, r = nodes2$r, g = nodes2$g, b = nodes2$b, a = 200)

    # For Axis 3 nodes,set heights and values:
    nodes3$height <- 5
    nodes3 <- within(nodes3,value <- cumsum(height))
    nodes3$value[1] <- 5
    for (i in 2:nrow(nodes3)){
      nodes3$value[i] <- nodes3$value[i-1] + nodes3$height[i-1] + 1
    }
    n3 <- data.frame(id = nodes3$id, axis = "c", label = nodes3$id, value = nodes3$value, height = nodes3$height, r = 100, g = 0, b = 100, a = 200)

    # Bind altogether to match web format
    hive_n <- rbind(n1,n2,n3)

    #Define links:
    temp <- nodes1[nodes1$height != 0,]
    hive_l <- unique(data.frame("start node label" = temp$id, "end node label" = temp$mode))
    temp <- nodes2[nodes2$height != 0,]
    hive_l <- rbind(hive_l,unique(data.frame("start node label" = temp$id, "end node label" = NOI)))
    # Write to tab delimited files
    write.table(hive_n,file = paste(NOI,"hive_n.txt"),sep = "\t",row.names = F)
    write.table(hive_l,file = paste(NOI,"hive_l.txt"),sep = "\t",row.names = F)

  }
}
