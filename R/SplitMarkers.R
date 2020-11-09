# object is a Seurat 3.0 object, function will operate on the active Identity slot
# this is formatted specifically for the HBEC object, which has a 'Condition' metadata slot
# GOI should be a list of genes of interest
# ... arguments to pass to FindMarkers

SplitMarkers <- function(object,GOI,...){
  # Stash idents and make new identities which identify each by their 'Condition' metadata
  celltypes <- as.character(unique(Idents(object)))
  # If using parcellations which are not conserved...
  celltypes <- names(which(rowSums(table(Idents(object),object$organ) != 0) == length(unique(object$organ))))
  
  celltypes.lung <- paste(celltypes, 'lung', sep = '_')
  celltypes.liver <- paste(celltypes, 'liver', sep = '_')
 
  object$celltype.condition <- paste(Idents(object), object$organ, sep = "_")
  object$celltype <- Idents(object)
  Idents(object) <- "celltype.condition"
  
  # Identify which genes, in each cell population, have an adjusted p-value < 0.05 based on a Wilcoxon rank test when compared to Mock
  genes <- GOI[GOI %in% rownames(object)]
  diff.p <- data.frame()
  for (i in 1:length(celltypes)){
    temp <- NULL
    try(
      temp <- FindMarkers(object,
                          ident.1 = celltypes.liver[i],
                          ident.2 = celltypes.lung[i],
                          verbose = FALSE,
                          features = genes,
                          min.cells.group = 0,...)
    )
    if(!is.null(temp)){
      temp2 <- subset(temp, p_val_adj < 0.05)
      if (nrow(temp2)>0){
        temp2$gene <- rownames(temp2)
        temp2$cell <- celltypes[i]
        diff.p <- rbind(diff.p, temp2)
      }
    }
  }
  diff.p$cell.gene <- paste(diff.p$cell,diff.p$gene,sep = '.')
 
  # Tag with condition
  diff.p$Comparison <- 'Liver.vs.Lung'
  
  differential.genes <- diff.p
  
  # Erase rownames which are nonsense
  rownames(differential.genes) <- NULL
  return(differential.genes)
}
