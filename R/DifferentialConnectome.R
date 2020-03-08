#' DifferentialConnectome
#'
#' Currently in beta testing. Creates a fold-change connectome from two input connectomes, generally unfiltered.  Must be node-aligned, from the same reference mapping, and unfiltered. ('edge' columns must contain identical entries, though not necessarily in the same order.)
#'
#' @param connect.1 A connectome from a system
#' @param connect.2 A connectome from a different system, to be compared to connect.1
#' @param min.pct Default 0.1. Threshold to return clusterwise observations for both ligand and receptor. Only needs to be satisfied in connect.1 OR in connect.2.

#' @export

DifferentialConnectome <- function(connect.1, connect.2,min.pct = 0.1){
  require(gtools)
  base1 <- connect.1
  base2 <- connect.2
  # Make same orientation of rows
  base1 <- base1 %>% arrange(desc(edge))
  base2 <- base2 %>% arrange(desc(edge))
  # Test to make sure identical
  if (sum(base1$edge != base2$edge)){stop("\nConnectomes do not have identical edges.")}
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
# Score
out$score <- abs(out$ligand.norm.lfc) * abs(out$recept.norm.lfc)

# Remove nonsense values (0 to 0 and or non-mapped 'NA' values - not useful for differential testing)
out <- subset(out,pct.source.1 > 0 | pct.source.2 > 0)
out <- subset(out, pct.target.1 > 0 | pct.target.2 > 0)

# Subset based on min.pct
if (!is.null(min.pct)){
  out <- subset(out,pct.source.1 > min.pct | pct.source.2 > min.pct)
  out <- subset(out, pct.target.1 > min.pct | pct.target.2 > min.pct)
}

return(out)
}
