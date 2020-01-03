#' Convert_to_Homologues
#'
#' This function converts a MARKER LIST to a homologue object.  Can be used to convert in any direction between mouse, rat, pig, and human.
#' @param mark A marker list (the output from FindMarkers or FindAllMarkers in Seurat)
#' @param old.species The original species of the object (can be 'mouse','rat','pig',or 'human')
#' @param new.species The desired species mapping. Defaults to 'human'
#' @export

Convert_Markers_to_Homologues <- function(mark,old.species,new.species = "human"){
  hom <- subset(translator,translator$query.species %in% old.species & translator$match.species == new.species)
  hom.use <- names(which(table(hom$query.gene.symbol) == 1)) # 1:1 condition A
  hom.use2 <- names(which(table(hom$match.gene.symbol) == 1)) # 1:1 condtion B
  hom.list <- hom[hom$query.gene.symbol %in% hom.use & hom$match.gene.symbol %in% hom.use2,]
  #Raw data gene names
  old.genes <- mark$gene
  for (i in 1:length(old.genes)){
    if (old.genes[i] %in% hom.list$query.gene.symbol){
        old.genes[i] <- hom.list[hom.list$query.gene.symbol == old.genes[i],]$match.gene.symbol
    }
  }
  for (i in 1:length(old.genes)){
    if (toupper(old.genes[i]) %in% hom[hom$query.gene.symbol == old.genes[i],]$match.gene.symbol &
    old.genes[i] != toupper(old.genes[i])){
      index <-  hom[hom$query.gene.symbol == old.genes[i],]
      old.genes[i] <- unique(index[index$match.gene.symbol == toupper(old.genes[i]),]$match.gene.symbol)
    }
  }
  new.genes <- old.genes
  mark$gene <- new.genes
  return(mark)
}
