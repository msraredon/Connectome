#' Annotate
#'
#' This function adds a column called 'mode' with a category for each L-R signaling interaction.
#' These categorizations are based on literature review and expert opinion.
#' Many of these molecules can likely be classified under multiple modes --  currently this function does not support this.
#'
#' @param connectome A connectomic edgelist
#' @export

Annotate2 <- function(connectome){
  connectome$mode <- NA
# By pair:
  try(connectome[grep("NLGN|GAS6|NRG|ICAM|VCAM|SELL|CDH|NCAM|NXPH|SELP|SELE|PROS1",connectome$pair),]$mode <- "Cell-cell adhesion")
  try(connectome[grep("NPNT|CYR|CHAD|NID|HAS2|SPP1|CTGF",connectome$pair),]$mode <- "Matrix (assorted)")
  try(connectome[grep("CALR|RTN4|PSAP|LRPAP1|PSEN1|SCARB1|JMJD6",connectome$pair),]$mode <- "Intracellular trafficking")
  try(connectome[grep("WNT|RSPO|DKK",connectome$pair),]$mode <- "WNT")
# Ligand fine grain:
  try(connectome[grep("ADAM",connectome$ligand),]$mode <- "ADAM")
  try(connectome[grep("ADIPOQ",connectome$ligand),]$mode <- "Adiponectin")
  try(connectome[grep("APOE",connectome$ligand),]$mode <- "APOE")
  try(connectome[grep("SERPIN|A2M",connectome$ligand),]$mode <- "Protease inhibition")
  try(connectome[grep("ADM|AGT|EDN1|VIP|APLN|TFPI|PLAT|PLAU",connectome$ligand),]$mode <- "Vasoactive")
  try(connectome[connectome$ligand %in% c("F2","F8","F7"),]$mode <- "Vasoactive")
  try(connectome[grep("FPR",connectome$ligand),]$mode <- "FPR")
  try(connectome[grep("RELN|AGRN|THBS|OMG|HSPG",connectome$ligand),]$mode <- "Matrix glycoproteins")
  try(connectome[grep("AREG|EREG|HBEGF|TGFA",connectome$ligand),]$mode <- "EGF")
  try(connectome[grep("ARTN|BDNF|NRTN|NTF3",connectome$ligand),]$mode <- "Neurotrophins")
  try(connectome[grep("BMP|RGMA",connectome$ligand),]$mode <- "BMP")
  try(connectome[grep("BTLA|EDA|FASLG|TNF",connectome$ligand),]$mode <- "TNF")
  try(connectome[grep("C3|C5|CFH",connectome$ligand),]$mode <- "Complement")
  try(connectome[grep("CCL",connectome$ligand),]$mode <- "CCL")
  try(connectome[grep("CSF",connectome$ligand),]$mode <- "CSF")
  try(connectome[grep("SHH",connectome$ligand),]$mode <- "SHH")
  try(connectome[grep("EFN",connectome$ligand),]$mode <- "Ephrins")
  try(connectome[grep("FGF",connectome$ligand),]$mode <- "FGF")
  try(connectome[grep("VEGF|PGF",connectome$ligand),]$mode <- "VEGF")
  try(connectome[grep("TGFB|INHBB|INHBA",connectome$ligand),]$mode <- "TGFB")
  try(connectome[grep("VEGF",connectome$ligand),]$mode <- "VEGF")
  try(connectome[grep("RARRES",connectome$ligand),]$mode <- "RARRES")
  try(connectome[grep("MDK|PTN",connectome$ligand),]$mode <- "NEGF")
  try(connectome[grep("IGF",connectome$ligand),]$mode <- "IGF")
  try(connectome[grep("IL|LIF|OSM|CTF1",connectome$ligand),]$mode <- "Interleukins")
  try(connectome[grep("CXCL|PF4",connectome$ligand),]$mode <- "CXCL")
  try(connectome[grep("ANGPT",connectome$ligand),]$mode <- "ANGPT")
  try(connectome[grep("MMP",connectome$ligand),]$mode <- "MMP")
  try(connectome[grep("PDGF",connectome$ligand),]$mode <- "PDGF")
  try(connectome[grep("SLIT",connectome$ligand),]$mode <- "SLIT")
  try(connectome[grep("SEMA",connectome$ligand),]$mode <- "Semaphorins")
  try(connectome[grep("COL",connectome$ligand),]$mode <- "Collagens")
  try(connectome[grep("LAM|NTN",connectome$ligand),]$mode <- "Laminins")
  try(connectome[grep("FN1",connectome$ligand),]$mode <- "Fibronectin")
  try(connectome[grep("PENK",connectome$ligand),]$mode <- "Opiod")
  # Receptor fine grain:
  try(connectome[grep("NOTCH|JAG",connectome$receptor),]$mode <- "NOTCH")
  try(connectome[grep("MET",connectome$receptor),]$mode <- "MET")
  try(connectome[grep("LPAR1",connectome$receptor),]$mode <- 'Lysophosphatidic acid')
  try(connectome[grep("AMFR|FPR1|FPR2|ATRN",connectome$receptor),]$mode <- 'Chemotaxis')
  try(connectome[grep("C5AR1",connectome$receptor),]$mode <- 'Complement')
  try(connectome[grep("KIT",connectome$receptor),]$mode <- "KIT")
  try(connectome[grep("TRAF2",connectome$receptor),]$mode <- "TNF")
  try(connectome[grep("TFRC",connectome$receptor),]$mode <- "Transferrin")
# Retracted
  connectome <- connectome[!connectome$pair == 'NAMPT-INSR',]
# Uncategorized
  #connectome[is.na(connectome$mode),]$mode <- 'UNCAT'
# Make factor
  connectome$mode <- as.factor(connectome$mode)
# Return
  return(connectome)
}
