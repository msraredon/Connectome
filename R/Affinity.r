#' This function makes an affinity matrix from a connectomic file
#' @param connectome A connectomic edgelist
#' @param meth.use What to sum. 'num' == number of edges, 'weight' == normalized expression weight, 'weight_sc' [DEFAULT] == scaled expression weight
#' @param per.exp Percent expression threshold. Defaults to 0.0
#' @param nodes.include Character vector of nodes to perform affinity analysis over
#' @export

Affinity <- function(connectome,meth.use = "weight_sc",per.exp = 0.0,nodes.include){
  require(dplyr)
  connectome$edge <- as.factor(paste(connectome$source,connectome$target))
  rel <- connectome
  #if (meth.use == 'weight_sc'){
    #rel <- rel[rel$ligand.scale>0,]
    #rel <- rel[rel$recept.scale>0,]
  #}
  rel <- rel[rel$percent.source > per.exp,]
  rel <- rel[rel$percent.target > per.exp,]


  nodes <- nodes.include
  nodes <- sort(as.character(nodes)) #Alphabetize



  aff <- matrix(,nrow=length(nodes),ncol=length(nodes))
  rownames(aff) <- nodes
  colnames(aff) <- nodes
  for (i in 1:length(nodes)){
    for (j in 1:length(nodes)){
      temp <- subset(rel, rel$source == nodes[i] & rel$target == nodes[j])
      if (nrow(temp) == 0 ){
        aff[i,j] <- 0
      }
      if (meth.use == "num"){
              aff[i,j] <- nrow(temp)
      }
      if (meth.use == "weight"){
              aff[i,j] <- sum(temp$weight)
      }
      if (meth.use == "weight_sc"){
              aff[i,j] <- sum(temp$weight_sc)
      }
    }
  }
  return(aff)
}
