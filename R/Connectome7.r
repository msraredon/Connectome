#' Connectome Version 7
#'
#' Creates a connectomic edgelist from a Seurat 3.0 object
#'
#' @param object A Seurat object
#' @param use.supported_known Determines which set of ligand receptor pairs to use. Defaults to TRUE
#' @param use.supported If TRUE, overrides the above and includes all literature supported ligand receptor pairs.
#' @param use.all If TRUE, overrides the above and includes all possible ligand receptor pairs from the FANTOM5 database
#' @export
Connectome7 <- function(object,use.supported_known = F,use.supported = F,use.all = T,...){
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(plotrix)
  #load("AverageExpressionLog.Rdata") #Custom average expression function
  #load("AverageExpressionSE.Rdata") #Custom standard error function
  #load("PercentExpression.Rdata") #Custom percent expressing per cluster function

# Import ligand-receptor pairs and metadata
#ncomms8866 <- read.table("ncomms8866-s3.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
all.ligands <- ncomms8866[,2]
all.recepts <- ncomms8866[,4]
ncomms.supported <- ncomms8866[ncomms8866$Pair.Evidence == "literature supported",]
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

cluster.avgs_unfiltered <- cluster.avgs

cluster.avgs <- cluster.avgs[rowSums(cluster.avgs_unfiltered > 0) >= 1, ] #Removes the occasional zero that snuck through intiial filtering
cluster.avgs.raw <- cluster.avgs.raw[rowSums(cluster.avgs_unfiltered > 0) >= 1, ] #Removes the occasional zero that snuck through intiial filtering
cluster.avgs.scale <- cluster.avgs.scale[rowSums(cluster.avgs_unfiltered > 0) >= 1, ] #Removes the occasional zero that snuck through intiial filtering
cluster.avgs.SE <- cluster.avgs.SE[rowSums(cluster.avgs_unfiltered > 0) >= 1, ] #Removes the occasional zero that snuck through intiial filtering
cluster.avgs.raw.SE <- cluster.avgs.raw.SE[rowSums(cluster.avgs_unfiltered > 0) >= 1, ] #Removes the occasional zero that snuck through intiial filtering
cluster.avgs.scale.SE <- cluster.avgs.scale.SE[rowSums(cluster.avgs_unfiltered > 0) >= 1, ] #Removes the occasional zero that snuck through intiial filtering
cluster.pcts <- cluster.pcts[rowSums(cluster.avgs_unfiltered > 0) >= 1, ]


# Focus analysis on only genes where both ligand and receptor are expressed (excludes orphan L/R's)
# Note that single ligands can hit multiple receptors
ligand.indexes <- which(ligands %in% rownames(cluster.avgs))
recept.indexes <- which(recepts %in% rownames(cluster.avgs))
paired.indexes <- intersect(ligand.indexes,recept.indexes)
paired.ligands <- ligands[paired.indexes]
paired.recepts <- recepts[paired.indexes]

# Re-subset to acquire dimension-matched matrices
clust.ligands.paired <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
clust.recepts.paired <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 clust.ligands.paired[n,] <- as.matrix(cluster.avgs[paired.ligands[n],])
 clust.recepts.paired[n,] <- as.matrix(cluster.avgs[paired.recepts[n],])
}
clust.ligands.paired.raw <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
clust.recepts.paired.raw <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 clust.ligands.paired.raw[n,] <- as.matrix(cluster.avgs.raw[paired.ligands[n],])
 clust.recepts.paired.raw[n,] <- as.matrix(cluster.avgs.raw[paired.recepts[n],])
}
clust.ligands.paired.scale <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
clust.recepts.paired.scale <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 clust.ligands.paired.scale[n,] <- as.matrix(cluster.avgs.scale[paired.ligands[n],])
 clust.recepts.paired.scale[n,] <- as.matrix(cluster.avgs.scale[paired.recepts[n],])
}
##SE
clust.ligands.paired.SE <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
clust.recepts.paired.SE <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 clust.ligands.paired.SE[n,] <- as.matrix(cluster.avgs.SE[paired.ligands[n],])
 clust.recepts.paired.SE[n,] <- as.matrix(cluster.avgs.SE[paired.recepts[n],])
}
clust.ligands.paired.raw.SE <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
clust.recepts.paired.raw.SE <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 clust.ligands.paired.raw.SE[n,] <- as.matrix(cluster.avgs.raw.SE[paired.ligands[n],])
 clust.recepts.paired.raw.SE[n,] <- as.matrix(cluster.avgs.raw.SE[paired.recepts[n],])
}
clust.ligands.paired.scale.SE <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
clust.recepts.paired.scale.SE <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 clust.ligands.paired.scale.SE[n,] <- as.matrix(cluster.avgs.scale.SE[paired.ligands[n],])
 clust.recepts.paired.scale.SE[n,] <- as.matrix(cluster.avgs.scale.SE[paired.recepts[n],])
}
pct.ligands.paired <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
pct.recepts.paired <- matrix(NA, nrow=length(paired.ligands), ncol=ncol(cluster.avgs))
for (n in 1:length(paired.ligands)){
 pct.ligands.paired[n,] <- as.matrix(cluster.pcts[paired.ligands[n],])
 pct.recepts.paired[n,] <- as.matrix(cluster.pcts[paired.recepts[n],])
}
# Set row and column names of the above matrices
### These must be matrices since they have duplicate row names to allow a 1:1 matching btw ligand and receptor ###
rownames(clust.ligands.paired) <- paired.ligands
colnames(clust.ligands.paired) <- colnames(cluster.avgs)
rownames(clust.recepts.paired) <- paired.recepts
colnames(clust.recepts.paired) <- colnames(cluster.avgs)
rownames(clust.ligands.paired.raw) <- paired.ligands
colnames(clust.ligands.paired.raw) <- colnames(cluster.avgs.raw)
rownames(clust.recepts.paired.raw) <- paired.recepts
colnames(clust.recepts.paired.raw) <- colnames(cluster.avgs.raw)
rownames(clust.ligands.paired.scale) <- paired.ligands
colnames(clust.ligands.paired.scale) <- colnames(cluster.avgs.scale)
rownames(clust.recepts.paired.scale) <- paired.recepts
colnames(clust.recepts.paired.scale) <- colnames(cluster.avgs.scale)
#SE
rownames(clust.ligands.paired.SE) <- paired.ligands
colnames(clust.ligands.paired.SE) <- colnames(cluster.avgs.SE)
rownames(clust.recepts.paired.SE) <- paired.recepts
colnames(clust.recepts.paired.SE) <- colnames(cluster.avgs.SE)
rownames(clust.ligands.paired.raw.SE) <- paired.ligands
colnames(clust.ligands.paired.raw.SE) <- colnames(cluster.avgs.raw.SE)
rownames(clust.recepts.paired.raw.SE) <- paired.recepts
colnames(clust.recepts.paired.raw.SE) <- colnames(cluster.avgs.raw.SE)
rownames(clust.ligands.paired.scale.SE) <- paired.ligands
colnames(clust.ligands.paired.scale.SE) <- colnames(cluster.avgs.scale.SE)
rownames(clust.recepts.paired.scale.SE) <- paired.recepts
colnames(clust.recepts.paired.scale.SE) <- colnames(cluster.avgs.scale.SE)
#percents
rownames(pct.ligands.paired) <- paired.ligands
colnames(pct.ligands.paired) <- colnames(cluster.pcts)
rownames(pct.recepts.paired) <- paired.recepts
colnames(pct.recepts.paired) <- colnames(cluster.pcts)
# Generate full connectome
# This approach assumes that all ligands are felt by all cells capable of sensing them.
ligand_list <-row.names(clust.ligands.paired)
recept_list <-row.names(clust.recepts.paired)
### For each ligand-receptor pair, pull the source clusters and target clusters. Then, for each source cluster
### assume it touches all clusters expressing the matched receptor (target = targets).  Bind the ligand and
### receptor names to the dataframe, a paired name, and the ligand and receptor expression levels, and bind these together
### to form a single interactome.  Do this for every ligand-receptor pair and bind them together to get the complete connectome.
connectome = data.frame() #Initializes empty dataframe
for (n in 1:nrow(clust.ligands.paired)){
  #One signaling mechanism per row
  sources.data <- clust.ligands.paired[n,]
  targets.data <- clust.recepts.paired[n,]
  sources.data.raw <- clust.ligands.paired.raw[n,]
  targets.data.raw <- clust.recepts.paired.raw[n,]
  sources.data.scale <- clust.ligands.paired.scale[n,]
  targets.data.scale <- clust.recepts.paired.scale[n,]
  #SE
  sources.data.SE <- clust.ligands.paired.SE[n,]
  targets.data.SE <- clust.recepts.paired.SE[n,]
  sources.data.raw.SE <- clust.ligands.paired.raw.SE[n,]
  targets.data.raw.SE <- clust.recepts.paired.raw.SE[n,]
  sources.data.scale.SE <- clust.ligands.paired.scale.SE[n,]
  targets.data.scale.SE <- clust.recepts.paired.scale.SE[n,]
  # Percents
  sources.data.percents <- pct.ligands.paired[n,]
  targets.data.percents <- pct.recepts.paired[n,]

  sources.data_unfiltered <- sources.data
  targets.data_unfiltered <- targets.data

  sources.data <- sources.data[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data <- targets.data[targets.data_unfiltered>0] # This line removes zero expression values
  sources.data.raw <- sources.data.raw[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data.raw <- targets.data.raw[targets.data_unfiltered>0] # This line removes zero expression values
  sources.data.scale <- sources.data.scale[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data.scale <- targets.data.scale[targets.data_unfiltered>0] # This line removes zero expression values
  #SE
  sources.data.SE <- sources.data.SE[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data.SE <- targets.data.SE[targets.data_unfiltered>0] # This line removes zero expression values
  sources.data.raw.SE <- sources.data.raw.SE[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data.raw.SE <- targets.data.raw.SE[targets.data_unfiltered>0] # This line removes zero expression values
  sources.data.scale.SE <- sources.data.scale.SE[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data.scale.SE <- targets.data.scale.SE[targets.data_unfiltered>0] # This line removes zero expression values
  # percents
  sources.data.percents <- sources.data.percents[sources.data_unfiltered>0] # This line removes zero expression values
  targets.data.percents <- targets.data.percents[targets.data_unfiltered>0] # This line removes zero expression values

  sources <- names(sources.data)
  targets <- names(targets.data)
  interactome <- data.frame() #Initializes empty dataframe
  for (i in 1:length(sources)){
    int <- data.frame(source = sources[i],target = targets, ligand = ligand_list[n], receptor = recept_list[n],
      pair = paste(ligand_list[n],recept_list[n],sep="-"),
      ligand.expression = unname(sources.data[i]),
      recept.expression = unname(targets.data),
      ligand.scale = unname(sources.data.scale[i]),
      recept.scale = unname(targets.data.scale),
      ligand.raw = unname(sources.data.raw[i]),
      recept.raw = unname(targets.data.raw),
      ligand.exp.SE = unname(sources.data.SE[i]),
      recept.exp.SE = unname(targets.data.SE),
      ligand.scale.SE = unname(sources.data.scale.SE[i]),
      recept.scale.SE = unname(targets.data.scale.SE),
      ligand.raw.SE = unname(sources.data.raw.SE[i]),
      recept.raw.SE = unname(targets.data.raw.SE),
      percent.source = unname(sources.data.percents[i]),
      percent.target = unname(targets.data.percents)
    )
    interactome <- rbind(interactome,int)
  } #Generates the interactome matrix for a SINGLE ligand-receptor pair
  connectome <- rbind(connectome,interactome) #Consolidates into complete connectome
} #Consolidates into complete connectome
# Add weights
connectome$weight <- connectome$ligand.expression * connectome$recept.expression
connectome$weight_sc <- connectome$ligand.scale + connectome$recept.scale #Does this make sense?
connectome$weight_raw <- connectome$ligand.raw * connectome$recept.raw
return(connectome)
}
