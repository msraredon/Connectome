#' CreateConnectome
#'
#' Creates an unfiltered connectomic edgelist from a Seurat 3.0 object.  Allows orphan ligand and receptors and fully non-expressed pairs.
#'
#' @param object A Seurat object
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param LR.database Accepts either 'fantom5' or 'custom'. If custom, a dataframe must be provided to argument custom.list with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations.
#' @param max.cells.per.ident Default NULL. If a value is input, input object will be downsampled to requested number of cells per identity. This can greatly improve run-time.
#' @param include.putative Default TRUE. Includes ligand-receptor pairs deemed putative in FANTOM5 database.
#' @param include.rejected Default FALSE. If TRUE, includes gene pairs labeled "EXCLUDED" in FANTOM5 database.  See ncomms8866 .rda file for qualifications for exclusion.
#' @param p.values Default TRUE. Runs a Wilcoxon Rank test to calculate adjusted p-value for ligand and receptor expression within the input object. Change to FALSE for decreased run-time.
#' @param weight.definition.norm Method of edgeweight definition from normalized data slot. Either 'sum','mean',or 'product'. Defaults to 'product'. 'Sum' adds values from sending and receiving clusters, 'mean' averages them, and 'product' multiplies them.
#' @param weight.definition.scale Method of edgeweight definition from scaled data slot. Either 'sum','mean',or 'product'. Defaults to 'mean'. 'Sum' adds values from sending and receiving clusters, 'mean' averages them, and 'product' multiplies them.
#' @param custom.list Optional. A dataframe for custom mapping, with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations. If modal categorizations are unknown, fill with 'UNCAT' or similar placeholder.
#' @export

CreateConnectome <- function(object,
                             LR.database = 'fantom5',
                             species = NULL,
                             include.putative = T,
                             include.rejected = F,
                             p.values = T,
                             max.cells.per.ident = NULL,
                             weight.definition.norm = 'product',
                             weight.definition.scale = 'mean',
                             custom.list = NULL,...){
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(plotrix)

  if (LR.database == 'fantom5' & is.null(species)){stop("\nPlease select species for fantom5 mapping. Allows 'human','mouse','rat', or 'pig' ")}

  # Downsample input object
  if (!is.null(max.cells.per.ident)){
    object <- Seurat::SubsetData(object = object,max.cells.per.ident = max.cells.per.ident)
    #object <- subset(object, cells = WhichCells(object, downsample = max.cells.per.ident))
  }

  if (LR.database == 'fantom5'){
    # Load ground-truth database (FANTOM5, species-converted as appropriate, per methodlogy in Raredon et al 2019, DOI: 10.1126/sciadv.aaw3851)
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
      ligands <- as.character(lit.put[,2]) #Determines the ligand list to use
      recepts <- as.character(lit.put[,4]) #Determines the receptor list to use
      modes <- as.character(lit.put[,'mode'])
    } else {
      ligands <- as.character(lit[,2]) #Determines the ligand list to use
      recepts <- as.character(lit[,4]) #Determines the receptor list to use
      modes <- as.character(lit[,'mode'])
    }
    if (include.rejected){
      ligands <- as.character(fantom5[,2])
      recepts <- as.character(fantom5[,4])
      modes <- as.character(fantom5[,'mode'])
    }
  }

  if (LR.database == 'custom'){
    if (is.null(custom.list)){stop("\nCustom mapping requested. Please provide custom list of Ligand-Receptor interactions")}
    ligands <- as.character(custom.list[,1])
    recepts <- as.character(custom.list[,2])
    modes <- as.character(custom.list[,3])
  }

  # Identify ligands and receptors expressed in the object
  ligands.use <- intersect(ligands,rownames(object@assays$RNA))
  recepts.use <- intersect(recepts,rownames(object@assays$RNA))
  genes.use = union(ligands.use,recepts.use)

  # Create averages and other relevant cluster-wise metrics, of only these GOI
  cluster.avgs <- AverageExpression(object,features = genes.use, assays = 'RNA')$RNA
  #cluster.avgs.raw <- AverageExpression(object,features = genes.use,slot = 'counts', assays = 'RNA')$RNA
  cluster.avgs.scale <- AverageExpression(object,features = genes.use,slot = "scale.data", assays = 'RNA')$RNA
  cluster.pcts <- PercentExpression_v2(object,features = genes.use,slot = 'counts')$RNA #percent cells in cluster expressing greater than zero

  # Include Wilcoxon Rank P-values?
  if (p.values){
    message(paste("\nCalculating p-values using Wilcoxon Rank"))
    cluster.p.values <- FindAllMarkers(object,assay = 'RNA',features = genes.use, test.use = 'wilcox',
                                       logfc.threshold = 0,min.pct = 0,return.thresh = Inf)
    cluster.p.values$cell.gene <- paste(cluster.p.values$cluster,cluster.p.values$gene,sep = ' - ')
    cluster.p.values$gene.cell <- paste(cluster.p.values$gene,cluster.p.values$cluster,sep = ' - ')
  }else{}

  # Generate full connectome
  sources <- colnames(cluster.avgs)
  targets <- colnames(cluster.avgs)
  message(paste("\nGenerating Connectome"))
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
                           #ligand.raw = cluster.avgs.raw[ligands,][,sources[i]],
                           #recept.raw = cluster.avgs.raw[recepts,][,targets[j]],
                           percent.source = cluster.pcts[ligands,][,sources[i]],
                           percent.target = cluster.pcts[recepts,][,targets[j]]
      )
      temp <- rbind(temp,vector)
    }
    connectome <- rbind(connectome,temp)
    Sys.sleep(0.5)
    setTxtProgressBar(pb,i)
  }
  # Log transform the expression values (since they were averaged outside of log space)
  connectome$ligand.expression <- log1p(connectome$ligand.expression)
  connectome$recept.expression <- log1p(connectome$recept.expression)

  # Add weights (normalized slot)
  if (weight.definition.norm == 'sum'){
    connectome$weight_norm<- connectome$ligand.expression + connectome$recept.expression
  }else{
    if (weight.definition.norm == 'mean'){
      connectome$weight_norm<- rowMeans(connectome[c('ligand.expression','recept.expression')])
    }else{
      if (weight.definition.norm == 'product'){
        connectome$weight_norm<- connectome$ligand.expression * connectome$recept.expression
      }else{
        message(paste("\nNo appropriate parameter specified for weight.definition.norm"))
    }
  }}
  # Add weights (scaled slot)
  if (weight.definition.scale == 'sum'){
    connectome$weight_sc <- connectome$ligand.scale + connectome$recept.scale
  }else{
    if (weight.definition.scale == 'mean'){
      connectome$weight_sc <- rowMeans(connectome[c('ligand.scale','recept.scale')])
    }else{
      if (weight.definition.scale == 'product'){
        connectome$weight_sc <- connectome$ligand.scale * connectome$recept.scale
      }else{
        message(paste("\nNo appropriate parameter specified for weight.definition.scale"))
    }
  }}

  #Additional columns
  connectome$vector <- paste(connectome$source,connectome$target,sep = ' - ')
  connectome$edge <- paste(connectome$source,connectome$ligand,connectome$receptor,connectome$target,sep = ' - ')
  connectome$source.ligand <- paste(connectome$source,connectome$ligand,sep = ' - ')
  connectome$receptor.target <- paste(connectome$receptor,connectome$target,sep = ' - ')
  if(!is.null(species)){connectome$species <- species}

  # Add p-values if required
  if (p.values){
    message(paste("\nMapping ligand p-values"))
    connectome <- base::merge(connectome, cluster.p.values[,c('p_val_adj','cell.gene')], by.x = "source.ligand", by.y = "cell.gene", all = TRUE)

    message(paste("\nMapping receptor p-values"))
    connectome <- base::merge(connectome, cluster.p.values[,c('p_val_adj','gene.cell')], by.x = "receptor.target", by.y = "gene.cell", all = TRUE)

    #Correct rows with non-sensical mappings
    connectome <- connectome[!is.na(connectome$source),]
    connectome <- connectome[!is.na(connectome$target),]

    # Fix column names
    names(connectome)[names(connectome) == 'p_val_adj.x'] <- 'p_val_adj.lig'
    names(connectome)[names(connectome) == 'p_val_adj.y'] <- 'p_val_adj.rec'
  }

  #Reorganize for presentation
  connectome <- connectome[,c('source','target',colnames(connectome)[!(colnames(connectome) %in% c('source','target'))])]

  message(paste("\nConnectome generation complete"))
  return(connectome)
}
