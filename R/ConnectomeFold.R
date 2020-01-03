#' Connectome Fold Change Comparison
#'
#' Creates a fold-change connectome from two input connectomes, inlcuding MAX outputs for things present in only one
#'
#' @param connect.1 A connectome from a system
#' @param connect.2 A connectome from a different system, to be compared to connect.1

#' @export

ConnectomeFold <- function(connect.1, connect.2){
  base1 <- connect.1
  base2 <- connect.2
  # Make same orientation of rows
  base1 <- base1 %>% arrange(desc(edge))
  base2 <- base2 %>% arrange(desc(edge))
  # Build inital output
  out <- data.frame(source = base1$source,
                    target = base1$target,
                    ligand = base1$ligand,
                    receptor = base1$receptor,
                    pair = base1$pair,
                    edge = base1$edge,
                    ligand.norm.lfc = foldchange2logratio(foldchange(base2$ligand.expression,base1$ligand.expression)),
                    recept.norm.lfc = foldchange2logratio(foldchange(base2$recept.expression,base1$recept.expression)),
                    pct.source.lfc = foldchange2logratio(foldchange(base2$percent.source,base1$percent.source)),
                    pct.target.lfc = foldchange2logratio(foldchange(base2$percent.target,base1$percent.target)),
                    pct.source.1 = base1$percent.source,
                    pct.source.2 = base2$percent.source,
                    pct.target.1 = base1$percent.target,
                    pct.target.2 = base2$percent.target)
  out$score <- abs(out$ligand.norm.lfc) + abs(out$recept.norm.lfc)
  # Correct NA for ligand fold change
  #out[which(is.na(out$ligand.norm.lfc) & !is.na(out$pct.source.1)),]$ligand.norm.lfc <- 'NEG'
  #out[which(is.na(out$ligand.norm.lfc) & !is.na(out$pct.source.2)),]$ligand.norm.lfc <- 'POS'
  # Correct NA for receptor fold change
  #out[which(is.na(out$recept.norm.lfc) & !is.na(out$pct.target.1)),]$recept.norm.lfc <- 'NEG'
  #out[which(is.na(out$recept.norm.lfc) & !is.na(out$pct.target.2)),]$recept.norm.lfc <- 'POS'
return(out)
}
