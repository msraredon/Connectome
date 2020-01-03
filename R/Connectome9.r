#' Connectome Version 9
#'
#' Creates a connectomic edgelist from a Seurat 3.0 object.  Allows orphan ligand and receptors and fully non-expressed pairs. Useful for comparative connectomics.
#'
#' @param object A Seurat object
#' @param species The species of the object that is being processed.  Requires input, and allows 'human','mouse','rat', or 'pig'
#' @param include.putative Default TRUE. Includes ligand-receptor pairs deemed putative in FANTOM5 database.
#' @param include.all Default FALSE. If TRUE, includes gene pairs labeled EXCLUDED in FANTOM5 database.  See ncomms8866 .rda file for qualifications for exclusion.
#' @param p.values Default FALSE. Runs a Wilcoxon Rank test to calculate adjusted p-value for ligand and receptor expression within the input object
#' @param max.cells.per.ident Default NULL. If a value is input, input object will be downsampled to requested number of cells per identity, to speed run time.
#' @param return.thresh Default 0.01. All ligands or receptors with a calculated p-value larger than this number will be reported as 1.
#' @export

Connectome9 <- function(object,species,include.putative = T,include.excluded = F,p.values = F,max.cells.per.ident = NULL,return.thresh = 0.01,...){
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(plotrix)

# Downsample input object
if (max.cells.per.ident){
  object <- SubsetData(object = object,max.cells.per.ident = max.cells.per.ident)
}

# Load ground-truth database (FANTOM5, species-converted if necessary)
    if (species == 'human'){
        fantom5 <- ncomms8866_human
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
lit.put <- fantom5[fantom5$Pair.Evidence %in% c("literature supported","putative"),]
lit <- fantom5[fantom5$Pair.Evidence %in% c("literature supported"),]
# Determine which list of pairs to perform analysis on:
if (include.putative){
  ligands <- lit.put[,2] #Determines the ligand list to use
  recepts <- lit.put[,4] #Determines the receptor list to use
  modes <- lit.put[,'mode']
}
else{
  ligands <- lit[,2] #Determines the ligand list to use
  recepts <- lit[,4] #Determines the receptor list to use
  modes <- lit[,'mode']
}
if (include.excluded){
  ligands <- fantom5[,2]
  recepts <- fantom5[,4]
  modes <- fantom5[,'mode']
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

# Include Wilcoxon Rank P-values?
if (p.values){
  message(paste("Calculating p-values using Wilcoxon Rank"))
  cluster.p.values <- FindAllMarkers(object,assay = 'RNA',features = genes.use, test.use = 'wilcox',
    logfc.threshold = 0,min.pct = 0,return.thresh = return.thresh)
  cluster.p.values$cell.gene <- paste(cluster.p.values$cluster,cluster.p.values$gene,sep = ' - ')
  cluster.p.values$gene.cell <- paste(cluster.p.values$gene,cluster.p.values$cluster,sep = ' - ')
}else{}

# Generate full connectome
sources <- colnames(cluster.avgs)
targets <- colnames(cluster.avgs)
message(paste("Generating Connectome"))
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
                                mode = modes,
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

# Add p-values if required
  if (p.values){
    connectome$lig.p <- 1
    connectome$rec.p <- 1
    message(paste("Mapping ligand p-values"))
    pb <- txtProgressBar(min = 0, max = length(cluster.p.values$cell.gene), initial = 0,style = 3)
    for (i in 1:length(cluster.p.values$cell.gene)){
      if (nrow(connectome[connectome$source.ligand == cluster.p.values$cell.gene[i],])>0 & cluster.p.values[i,]$p_val_adj != 1){
        connectome[connectome$source.ligand == cluster.p.values$cell.gene[i],]$lig.p <- cluster.p.values[i,]$p_val_adj
      }
      Sys.sleep(0.5)
      setTxtProgressBar(pb,i)
    }
    message(paste("Mapping receptor p-values"))
    pb <- txtProgressBar(min = 0, max = length(cluster.p.values$gene.cell), initial = 0,style = 3)
    for (i in 1:length(cluster.p.values$gene.cell)){
      if (nrow(connectome[connectome$receptor.target == cluster.p.values$gene.cell[i],])>0 & cluster.p.values[i,]$p_val_adj != 1){
        connectome[connectome$receptor.target == cluster.p.values$gene.cell[i],]$rec.p <- cluster.p.values[i,]$p_val_adj
      }
      Sys.sleep(0.5)
      setTxtProgressBar(pb,i)
    }
  }
message(paste("Connectome complete"))
return(connectome)
}
