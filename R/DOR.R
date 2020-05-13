#' DOR
#'
##' @param object A connectomic object
##' @param ident Cell identity of interest
##' @param features Features of interest to calculate DOR

DOR <- function(object,ident,features){
  require(dplyr)
  require(Seurat)

  cells <- WhichCells(object = object,idents = ident)
  features <- as.character(features)
  pseudo = 0.0
  if (length(features) > 1){
    TP <- rowSums(data.frame(object@assays$RNA@counts[features,cells]>0))
    FN <- rowSums(data.frame(object@assays$RNA@counts[features,cells] == 0))
    FP <- rowSums(data.frame(object@assays$RNA@counts[features,!colnames(object) %in% cells]>0))
    TN <- rowSums(data.frame(object@assays$RNA@counts[features,!colnames(object) %in% cells] == 0))
  }else{
    TP <- sum(object@assays$RNA@counts[features,cells]>0)
    FN <- sum(object@assays$RNA@counts[features,cells] == 0)
    FP <- sum(object@assays$RNA@counts[features,!colnames(object) %in% cells]>0)
    TN <- sum(object@assays$RNA@counts[features,!colnames(object) %in% cells] == 0)
  }
  DOR <- ((TP+pseudo)/(FP+pseudo))/((FN+pseudo)/(TN+pseudo))
  DOR <- log(DOR)
}
