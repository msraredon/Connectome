#' SingleCellConnectome
#'
#' Takes a Seurat object and creates a single-cell connectome (rows are ligand-receptor mechanisms and columns are cell-cell vectors)
#'
#' @param object A Seurat object
#' @param species The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'
#' @param LR.database Accepts either 'fantom5' or 'custom'. If custom, a dataframe must be provided to argument custom.list with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations.
#' @param max.cells.per.ident Default NULL. If a value is input, input object will be downsampled to requested number of cells per identity. This can greatly improve run-time.
#' @param min.cells.per.ident Default NULL. If a value is input, only cell populations meeting this threshold will be included in network analysis. Can limit high-variation effects from small clusters.
#' @param include.putative Default TRUE. Includes ligand-receptor pairs deemed putative in FANTOM5 database.
#' @param include.rejected Default FALSE. If TRUE, includes gene pairs labeled "EXCLUDED" in FANTOM5 database.  See ncomms8866 .rda file for qualifications for exclusion.
#' @param custom.list Optional. A dataframe for custom mapping, with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations. If modal categorizations are unknown, fill with 'UNCAT' or similar placeholder.
#' @export

SingleCellConnectome <- function(object,
                             LR.database = 'fantom5',
                             species = NULL,
                             include.putative = T,
                             include.rejected = F,
                             #p.values = T,
                             max.cells.per.ident = NULL,
                             min.cells.per.ident = NULL,
                             slot.use = 'data',
                             weight.definition = 'product',
                             custom.list = NULL,
                             #calculate.DOR = F,
                             ...){
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

  # Limit object to cell populations larger than requested minimum
  if (!is.null(min.cells.per.ident)){
    idents.include <- names(table(Idents(object)))[table(Idents(object)) > min.cells.per.ident]
    object <- subset(object,idents = idents.include)
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
      fantom <- lit.put
    } else {
      fantom <- lit
    }
    if (include.rejected){
      fantom <- fantom5
    }
  }

  if (LR.database == 'custom'){
    if (is.null(custom.list)){stop("\nCustom mapping requested. Please provide custom list of Ligand-Receptor interactions")}
    ligands <- as.character(custom.list[,1])
    recepts <- as.character(custom.list[,2])
    modes <- as.character(custom.list[,3])
    fantom <- data.frame(Ligand.ApprovedSymbol = ligands, Receptor.ApprovedSymbol = recepts)
  }

  # Identify ligands and receptors expressed in the object
  fantom.specific <- subset(fantom,Ligand.ApprovedSymbol %in% rownames(object) & Receptor.ApprovedSymbol %in% rownames(object))
  ligands <- fantom.specific$Ligand.ApprovedSymbol
  receptors <- fantom.specific$Receptor.ApprovedSymbol

  # This works for now BUT might not work flawlessly to pull rows because of the known grep bug for shortened strings.  Need to do differently in final.
  # Also should we do this as columns or rows? Is one more efficient than the other?

  # Ligand Map
  lig.map <- object@assays$RNA@data[ligands,]
  lig.map2 <- do.call(cbind, replicate(ncol(lig.map), lig.map, simplify=FALSE))
  # Receptor Map
  rec.map <- object@assays$RNA@data[receptors,]
  rec.map2 <- matchingR:::repcol(rec.map,ncol(rec.map))
  # Merged map (can be done with any operator, here is multiplication (preserves zeroes, quantitative))
  sc.connectome <- lig.map2*rec.map2

  # Name
  rownames(sc.connectome) <- paste(rownames(lig.map),rownames(rec.map),sep = '-')
  colnames(sc.connectome) <- paste(colnames(lig.map2),rep(colnames(rec.map),each = length(colnames(rec.map))),sep = '-')

  message(paste("\nSingle-Cell Connectome generation complete"))
  return(sc.connectome)
}
