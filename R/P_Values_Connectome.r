# This function calculates Wilcoxon-rank P-values for ligands and receptors of choice within a Seurat object.
# Requires a Seurat object and choice of L/R set
P_Values_Connectome <- function(object,use.supported_known = T,use.supported = F,use.all = F,...){
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(plotrix)
# Import ligand-receptor pairs and metadata
# ncomms8866 <- read.table("ncomms8866-s3.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)
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
genes.use.ligands <- rownames(object@raw.data)[rownames(object@raw.data) %in% ligands]
genes.use.recepts <- rownames(object@raw.data)[rownames(object@raw.data) %in% recepts]
genes.use = union(genes.use.ligands,genes.use.recepts)
# Create p-values of only these GOI
p_values <- FindAllMarkers(object,genes.use = genes.use,logfc.threshold = 0,min.pct = 0,return.thresh = Inf)
return(p_values)
}
