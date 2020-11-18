#' EvenSplit
#'
#' Splits a Seurat object using SplitObject, but makes sure that for each identity, an even number of cells are sampled.
#'
#' @param connectome A connectomic edgelist
#' @param split.by Parameter by which to split the object
#' @export

EvenSplit <- function(object,split.by){
  # Number of cells grabb-able from each condition, per cell type
  nums <- table(Idents(object),object@meta.data[[split.by]])
  mins <- apply(nums, 1, FUN=min)
  split.object <- SplitObject(object,split.by = split.by)

  # Equalize the grabs across condition
  split.down <- list()
  for (i in 1:length(split.object)){
    pile <- list()
    for (j in 1:length(mins)){
      temp <- subset(split.object[[i]],idents = names(mins[j]))
      #temp <- SubsetData(temp, max.cells.per.ident = mins[j])
      temp <- subset(temp, cells = WhichCells(temp,downsample = mins[j]))
      pile[[j]] <- temp
    }
    split.down[[i]] <- pile
  }
  # Squash into a list of objects
  for (i in 1:length(split.down)){
    split.down[[i]] <- merge(split.down[[i]][[1]],split.down[[i]][2:length(split.down[[i]])])
  }
  names(split.down) <- names(split.object)
  return(split.down)
}
