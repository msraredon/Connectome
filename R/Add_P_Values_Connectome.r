#' This function takes the p-values calculated with 'P_Values_Connectome'
#' and adds them to an already-calculated connectome.  These two tasks are split up because the
#' unthresholded Wilcoxon rank test can take a long time
#' @param p_values The output from 'P_Values_Connectome'
#' @export

Add_P_Values_Connectome <- function(connectome,p_values){
  library(dplyr)
  library(tibble)
  library(plotrix)
# Identify things in connectome
sources <- as.character(unique(droplevels(connectome$source)))
ligands <- as.character(unique(droplevels(connectome$ligand)))
targets <- as.character(unique(droplevels(connectome$target)))
recepts <- as.character(unique(droplevels(connectome$receptor)))
connectome$lig.p <- 1
for (n in 1:length(sources)){
  for (m in 1:length(ligands)){
    p <- p_values[p_values$cluster == sources[n] & p_values$gene == ligands[m],]
    if (nrow(connectome[connectome$source == sources[n] & connectome$ligand == ligands[m],]) != 0 & nrow(p) == 1){
      connectome[connectome$source == sources[n] & connectome$ligand == ligands[m],]$lig.p <- p$p_val_adj
    }
  }
}
connectome$rec.p <- 1
for (n in 1:length(targets)){
  for (m in 1:length(recepts)){
    p <- p_values[p_values$cluster == targets[n] & p_values$gene == recepts[m],]
    if (nrow(connectome[connectome$target == targets[n] & connectome$receptor == recepts[m],]) != 0 & nrow(p) == 1){
      connectome[connectome$target == targets[n] & connectome$receptor == recepts[m],]$rec.p <- p$p_val_adj
    }
  }
}
return(connectome)
}
