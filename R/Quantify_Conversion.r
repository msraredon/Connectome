#' Quantify_Conversion
#'
#' This function quantifies the homologue conversion of a GENE LIST to a homologue object.
#' @param gene.list A gene list
#' @param old.species The original species of the object (can be 'mouse','rat','pig',or 'human')
#' @param new.species The desired species mapping. Defaults to 'human'
#' @export

Quantify_Conversion <- function(gene.list,old.species,new.species = "human"){
  hom <- subset(translator,translator$query.species %in% old.species & translator$match.species == new.species)
  hom.use <- names(which(table(hom$query.gene.symbol) == 1)) # 1:1 condition A
  hom.use2 <- names(which(table(hom$match.gene.symbol) == 1)) # 1:1 condtion B
  hom.list <- hom[hom$query.gene.symbol %in% hom.use & hom$match.gene.symbol %in% hom.use2,]
  #Raw data gene names
  old.genes <- gene.list
  output <- data.frame(old.genes = gene.list,new.genes = 0,conversion = 0)
  for (i in 1:length(old.genes)){
    if (old.genes[i] %in% hom.list$query.gene.symbol){
        output[output$old.genes == old.genes[i],]$conversion <- 1
        #output[output$old.genes == old.genes[i],]$new.gene <- hom.list[hom.list$query.gene.symbol == old.genes[i],]$match.gene.symbol
    }
  }
  for (i in 1:length(old.genes)){
    if (toupper(old.genes[i]) %in% hom[hom$query.gene.symbol == old.genes[i],]$match.gene.symbol &
    old.genes[i] != toupper(old.genes[i])){
      index <-  hom[hom$query.gene.symbol == old.genes[i],]
      output[output$old.genes == old.genes[i],]$conversion <- 1
      #output[output$old.genes == old.genes[i],]$new.gene <- unique(index[index$match.gene.symbol == toupper(old.genes[i]),]$match.gene.symbol)
    }
  }
  return(output)
}
