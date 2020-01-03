#' Conserved_Edges
#'
#' This function takes multiple connectomic edgelists as input and returns a single master connectome identifying conserved edges and average properties
#'
#' @param con.list A list of connectomes to calculate conservation across
#' @export

Conserved_Edges <- function(con.list){
  require(BisRNA)
  require(data.table)
  edges.use <- c()
  master <- data.frame()
  # make edgetype if it doesn't exist and bind all connectomes together
  for (i in 1:length(con.list)) {
    con.list[[i]]$edgetype <- paste(con.list[[i]]$source,con.list[[i]]$pair,con.list[[i]]$target, sep = "-")
    con.list[[i]]$cellpair <- paste(con.list[[i]]$source,con.list[[i]]$target, sep = "-")
    master <- rbind(master,con.list[[i]])
  }
  master <- as.data.table(master)
  master[ , master.lig.wt      := mean(ligand.scale),   by = edgetype]
  master[ , master.rec.wt      := mean(recept.scale),   by = edgetype]
  master[ , master.score      := mean(weight_sc),   by = edgetype]
  master[ , master.lig.pct      := mean(percent.source),   by = edgetype]
  master[ , master.rec.pct      := mean(percent.target),   by = edgetype]
  # Correct 0 p-values
  master[master$lig.p == 0,]$lig.p <- 1e-305
  master[master$rec.p == 0,]$rec.p <- 1e-305
  # Combine p-values
  master[ , lig.p.fisher:= fisher.method(lig.p), by = edgetype]
  master[ , rec.p.fisher:= fisher.method(rec.p), by = edgetype]
  # KS p-value computation
  interactions_rank_test <- function(ranks, NMech = 100) {
   ks.test(x=ranks, y = stepfun(x = c(1:NMech), y = c(0:NMech)/NMech), alternative = "greater")$p.value
  }
  # Identify cell-cell pairs and species present
  allpairs <- unique(master$cellpair)
  allspecies <- unique(master$species)
  # Initialize rankings
  master$interaction.rank <- as.numeric(0)
  # For each cell pair:
  for (j in 1:length(allpairs)){
    # Create normalized rankings of each interaction, in each species:
    for (k in 1:length(allspecies)){
    master[master$cellpair == allpairs[j] & master$species == allspecies[k],][order(-weight_sc),]$interaction.rank <-
      as.numeric(rownames(master[master$cellpair == allpairs[j] & master$species == allspecies[k],][order(-weight_sc),]))/nrow(master[master$cellpair == allpairs[j] & master$species == allspecies[k],])
    }
  }
  # Normalize to 100 (NMech)
  master$interaction.rank <- master$interaction.rank*100
  # Calculate KS p-value across species
  master$p.KS <- as.numeric(NA)
  master[ , p.KS  := interactions_rank_test(interaction.rank),  by = edgetype]
  return(master)
}
