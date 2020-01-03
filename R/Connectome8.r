#' Connectome Version 8
#'
#' Creates a connectomic edgelist from a Seurat 3.0 object.  Allows orphan ligand and receptors and fully non-expressed pairs. Useful for comparative connectomics.
#'
#' @param object A Seurat object
#' @param species The species of the object that is being processed.  Allows 'human','mouse','rat', or 'pig'
#' @param use.supported_known Determines which set of ligand receptor pairs to use. Defaults to TRUE
#' @param use.supported If TRUE, overrides the above and includes all literature supported ligand receptor pairs.
#' @param use.all If TRUE, overrides the above and includes all possible ligand receptor pairs from the FANTOM5 database
#' @export
Connectome8 <- function(object,species,use.supported_known = F,use.supported = F,use.all = T,...){
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(plotrix)

# Load ground-truth database (FANTOM5, species-converted if necessary)
    if (species == 'human'){
        fantom5 <- ncomms8866
    }
    if (species == 'mouse'){
      fantom5 <- ncomms8866_mouse
    }
    if (species == 'rat'){
      fantom5 <- ncomms8866_rat
    }
    if (species == 'pig'){
      fantom5 <- ncomms8866_pig
    }

# Import ligand-receptor pairs and metadata
all.ligands <- fantom5[,2]
all.recepts <- fantom5[,4]
ncomms.supported <- fantom5[fantom5$Pair.Evidence == "literature supported",]
ncomms.supported.known <- ncomms.supported[ncomms.supported$Pair.Source == "known",]
# Determine which list of pairs to perform analysis on:
if (use.supported_known){
ligands <- ncomms.supported.known[,2] #Determines the ligand list to use
recepts <- ncomms.supported.known[,4] #Determines the receptor list to use
}
if (use.supported){
ligands <- ncomms.supported[,2] #Determines the ligand list to use
recepts <- ncomms.supported[,4] #Determines the receptor list to use
}
if (use.all){
  ligands <- all.ligands
  recepts <- all.recepts
}

# Identify ligands and receptors expressed in the object
genes.use.ligands <- intersect(ligands,rownames(object@assays$RNA))
genes.use.recepts <- intersect(recepts,rownames(object@assays$RNA))
genes.use = union(genes.use.ligands,genes.use.recepts)

# Create averages and other relevant cluster-wise metrics, of only these GOI
cluster.avgs <- AverageExpression(object,features = genes.use, assays = 'RNA')$RNA
cluster.avgs.raw <- AverageExpression(object,features = genes.use,slot = 'counts', assays = 'RNA')$RNA
cluster.avgs.scale <- AverageExpression(object,features = genes.use,slot = "scale.data", assays = 'RNA')$RNA

cluster.avgs.SE <- AverageExpressionSE_v2(object,features = genes.use,assays = 'RNA')$RNA
cluster.avgs.raw.SE <- AverageExpressionSE_v2(object,features = genes.use,slot = 'counts',assays = 'RNA')$RNA
cluster.avgs.scale.SE <- AverageExpressionSE_v2(object,features = genes.use,slot = "scale.data",assays = 'RNA')$RNA
cluster.pcts <- PercentExpression_v2(object,features = genes.use,slot = 'counts')$RNA #percent cells in cluster expressing greater than zero

# Generate full connectome
sources <- colnames(cluster.avgs)
targets <- colnames(cluster.avgs)
pb <- txtProgressBar(min = 0, max = length(sources), initial = 0,style = 3)

connectome <- data.frame()
  for (i in 1:length(sources)){
    temp <- data.frame()
    for (j in 1:length(targets)){
      vector <- data.frame(source = sources[i],
                                target = targets[j],
                                ligand = ligands,
                                receptor = recepts,
                                pair = paste(ligands,recepts,sep = ' - '),
                                ligand.expression = cluster.avgs[ligands,][,sources[i]],
                                recept.expression = cluster.avgs[recepts,][,targets[j]],
                                ligand.scale = cluster.avgs.scale[ligands,][,sources[i]],
                                recept.scale = cluster.avgs.scale[recepts,][,targets[j]],
                                ligand.raw = cluster.avgs.raw[ligands,][,sources[i]],
                                recept.raw = cluster.avgs.raw[recepts,][,targets[j]],
                                ligand.exp.SE = cluster.avgs.SE[ligands,][,sources[i]],
                                recept.exp.SE = cluster.avgs.SE[recepts,][,targets[j]],
                                ligand.scale.SE = cluster.avgs.scale.SE[ligands,][,sources[i]],
                                recept.scale.SE = cluster.avgs.scale.SE[recepts,][,targets[j]],
                                ligand.raw.SE = cluster.avgs.raw.SE[ligands,][,sources[i]],
                                recept.raw.SE = cluster.avgs.raw.SE[recepts,][,targets[j]],
                                percent.source = cluster.pcts[ligands,][,sources[i]],
                                percent.target = cluster.pcts[recepts,][,targets[j]]
                                )
      temp <- rbind(temp,vector)
    }
    connectome <- rbind(connectome,temp)
    Sys.sleep(0.5)
    setTxtProgressBar(pb,i)
  }



# Add weights and additional bulk columns
connectome$weight <- connectome$ligand.expression + connectome$recept.expression
connectome$weight_sc <- connectome$ligand.scale + connectome$recept.scale
connectome$weight_raw <- connectome$ligand.raw + connectome$recept.raw
connectome$vector <- paste(connectome$source,connectome$target,sep = ' - ')
connectome$edge <- paste(connectome$source,connectome$ligand,connectome$receptor,connectome$target,sep = ' - ')
connectome$source.ligand <- paste(connectome$source,connectome$ligand,sep = ' - ')
connectome$receptor.target <- paste(connectome$receptor,connectome$target,sep = ' - ')
return(connectome)
}
