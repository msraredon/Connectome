#' Annotate
#'
#' This function adds a column called 'mode' with a category for each L-R signaling interaction.
#' These categorizations are based on literature review and should be updated and broken down with further specificity as more information is learned.
#' Many of these molecules can likely be classified under multiple modes, currently this function does not support this.
#'
#' @param connectome A connectomic edgelist
#' @export

Annotate <- function(connectome){
  connectome$mode <- NA
  try(connectome[grep("EDA|ANXA1",connectome$pair),]$mode <- "Dev & Diff")
  try(connectome[grep("WNT|RSPO|DKK",connectome$pair),]$mode <- "Wnt")
  try(connectome[grep("SHH",connectome$pair),]$mode <- "Shh")
  try(connectome[grep("NOTCH|JAG",connectome$pair),]$mode <- "Notch")
  try(connectome[grep("KIT",connectome$pair),]$mode <- "Kit")
  try(connectome[grep("SEMA|EFN|SLIT|CCL|CX|RARR|HEBP|FPR2|RGMA",connectome$pair),]$mode <- "Chemotaxis")
  try(connectome[grep("SLIT",connectome$pair),]$mode <- "Slit")
  try(connectome[grep("SEMA",connectome$pair),]$mode <- "Semaphorins")
  try(connectome[grep("EFN",connectome$pair),]$mode <- "Ephrins")
  try(connectome[grep("COL|LAM|FN1|NPNT|CYR|HSP|CHAD|RELN|THBS|OMG|NTN|NID|HAS2|SPP1",connectome$pair),]$mode <- "Matrix-Cell Interactions")
  try(connectome[grep("ANGPT|GF|BMP|INHB|BDNF|RTN|PTN|CSF|MDK|EREG",connectome$pair),]$mode <- "Growth factors")
  try(connectome[grep("ADAM|NAMPT|GSTP|PLD|PSAP|SERPIN|A2M|PLAU|MMP",connectome$pair),]$mode <- "Enzymes & Inhibition")
  try(connectome[grep("APOE|ASIP|CALR|GPI|TF|AGRN|LIPH",connectome$pair),]$mode <- "Paracrine Signals")
  try(connectome[grep("ADIPOQ|OSTN|INHA",connectome$pair),]$mode <- "Hormone")
  try(connectome[grep("ADM|EDN|APL|VIP|AGT|TFPI|PLAT",connectome$pair),]$mode <- "Vasoactive & Coag")
  try(connectome[connectome$ligand %in% c("F8","F7"),]$mode <- "Vasoactive & Coag")
  try(connectome[grep("IL|OSM|PF4|TNF|C5|SERPING|C3|HLA|LIF",connectome$pair),]$mode <- "Cytokine & Immune")
  try(connectome[grep("NLGN|GAS6|NRG|ICAM|VCAM|SELL|CDH|NCAM|NXPH|SELP|SELE|PROS1",connectome$pair),]$mode <- "Cell-cell adhesion")
  try(connectome[grep("LRPAP|PSEN1",connectome$pair),]$mode <- "Other")
  connectome$mode <- as.factor(connectome$mode)
  return(connectome)
}
