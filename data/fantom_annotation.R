setwd("~/Documents/GitHub/connectome/data")
load("~/Documents/GitHub/connectome/data/ncomms8866_orig.rda")

Annotate3 <- function(fantom){
  fantom$mode <- NA
  # By pair:
  try(fantom[grep("NLGN|GAS6|NRG|ICAM|VCAM|SELL|CDH|NCAM|NXPH|SELP|SELE|PROS1",fantom$Pair.Name),]$mode <- "Cell-cell adhesion")
  try(fantom[grep("NPNT|CYR|CHAD|NID|HAS2|SPP1|CTGF",fantom$Pair.Name),]$mode <- "Matrix (assorted)")
  try(fantom[grep("CALR|RTN4|PSAP|LRPAP1|PSEN1|SCARB1|JMJD6",fantom$Pair.Name),]$mode <- "Intracellular trafficking")
  try(fantom[grep("WNT|RSPO|DKK",fantom$Pair.Name),]$mode <- "WNT")
  # Ligand fine grain:
  try(fantom[grep("ADAM",fantom$Ligand.ApprovedSymbol),]$mode <- "ADAM")
  try(fantom[grep("ADIPOQ",fantom$Ligand.ApprovedSymbol),]$mode <- "Adiponectin")
  try(fantom[grep("APOE",fantom$Ligand.ApprovedSymbol),]$mode <- "APOE")
  try(fantom[grep("SERPIN|A2M",fantom$Ligand.ApprovedSymbol),]$mode <- "Protease inhibition")
  try(fantom[grep("ADM|AGT|EDN1|VIP|APLN|TFPI|PLAT|PLAU|C1Q|KNG1",fantom$Ligand.ApprovedSymbol),]$mode <- "Vasoactive")
  try(fantom[fantom$Ligand.ApprovedSymbol %in% c("F2","F8","F7",'VWF'),]$mode <- "Vasoactive")
  try(fantom[grep("RELN|AGRN|THBS|OMG|HSPG|VTN",fantom$Ligand.ApprovedSymbol),]$mode <- "Matrix glycoproteins")
  try(fantom[grep("AREG|EREG|HBEGF|TGFA",fantom$Ligand.ApprovedSymbol),]$mode <- "EGF")
  try(fantom[fantom$Ligand.ApprovedSymbol %in% c("EGF"),]$mode <- "EGF")
  try(fantom[grep("ARTN|BDNF|NRTN|NTF3",fantom$Ligand.ApprovedSymbol),]$mode <- "Neurotrophins")
  try(fantom[grep("BMP|RGMA",fantom$Ligand.ApprovedSymbol),]$mode <- "BMP")
  try(fantom[grep("BTLA|EDA|FASLG|TNF",fantom$Ligand.ApprovedSymbol),]$mode <- "TNF")
  try(fantom[grep("C3|C5|CFH",fantom$Ligand.ApprovedSymbol),]$mode <- "Complement")
  try(fantom[grep("CCL",fantom$Ligand.ApprovedSymbol),]$mode <- "CC")
  try(fantom[grep("CSF",fantom$Ligand.ApprovedSymbol),]$mode <- "CSF")
  try(fantom[grep("SHH|DHH",fantom$Ligand.ApprovedSymbol),]$mode <- "Hedgehog")
  try(fantom[grep("EFN",fantom$Ligand.ApprovedSymbol),]$mode <- "Ephrins")
  try(fantom[grep("FGF",fantom$Ligand.ApprovedSymbol),]$mode <- "FGF")
  try(fantom[grep("VEGF|PGF",fantom$Ligand.ApprovedSymbol),]$mode <- "VEGF")
  try(fantom[grep("TGFB|INHBB|INHBA",fantom$Ligand.ApprovedSymbol),]$mode <- "TGFB")
  try(fantom[grep("RARRES",fantom$Ligand.ApprovedSymbol),]$mode <- "RARRES")
  try(fantom[grep("MDK|PTN",fantom$Ligand.ApprovedSymbol),]$mode <- "NEGF")
  try(fantom[grep("IGF",fantom$Ligand.ApprovedSymbol),]$mode <- "IGF")
  try(fantom[grep("IL|LIF|OSM|CTF1",fantom$Ligand.ApprovedSymbol),]$mode <- "Interleukins")
  try(fantom[grep("CXCL|PF4",fantom$Ligand.ApprovedSymbol),]$mode <- "CXCL")
  try(fantom[grep("ANGPT",fantom$Ligand.ApprovedSymbol),]$mode <- "ANGPT")
  try(fantom[grep("MMP|TIMP",fantom$Ligand.ApprovedSymbol),]$mode <- "MMP")
  try(fantom[grep("PDGF",fantom$Ligand.ApprovedSymbol),]$mode <- "PDGF")
  try(fantom[grep("SLIT",fantom$Ligand.ApprovedSymbol),]$mode <- "SLIT")
  try(fantom[grep("SEMA",fantom$Ligand.ApprovedSymbol),]$mode <- "Semaphorins")
  try(fantom[grep("COL",fantom$Ligand.ApprovedSymbol),]$mode <- "Collagens")
  try(fantom[grep("LAM|NTN",fantom$Ligand.ApprovedSymbol),]$mode <- "Laminins")
  try(fantom[grep("FN1",fantom$Ligand.ApprovedSymbol),]$mode <- "Fibronectin")
  try(fantom[grep("PENK",fantom$Ligand.ApprovedSymbol),]$mode <- "Opiod")
  try(fantom[grep("HLA-",fantom$Ligand.ApprovedSymbol),]$mode <- "MHC")
  try(fantom[grep("UBA52",fantom$Ligand.ApprovedSymbol),]$mode <- "Ubiquitin")
  # Receptor fine grain:
  try(fantom[grep("NOTCH|JAG",fantom$Receptor.ApprovedSymbol),]$mode <- "NOTCH")
  try(fantom[grep("MET",fantom$Receptor.ApprovedSymbol),]$mode <- "MET")
  try(fantom[grep("LPAR1",fantom$Receptor.ApprovedSymbol),]$mode <- 'Lysophosphatidic acid')
  try(fantom[grep("AMFR|FPR1|FPR2|ATRN",fantom$Receptor.ApprovedSymbol),]$mode <- 'Chemotaxis')
  try(fantom[grep("C5AR1",fantom$Receptor.ApprovedSymbol),]$mode <- 'Complement')
  try(fantom[grep("KIT",fantom$Receptor.ApprovedSymbol),]$mode <- "KIT")
  try(fantom[grep("BMPR",fantom$Receptor.ApprovedSymbol),]$mode <- "BMP")
  try(fantom[grep("TRAF2",fantom$Receptor.ApprovedSymbol),]$mode <- "TNF")
  try(fantom[grep("TFRC",fantom$Receptor.ApprovedSymbol),]$mode <- "Transferrin")
  try(fantom[grep("EGFR",fantom$Receptor.ApprovedSymbol),]$mode <- "EGF")
  try(fantom[grep("FGFR",fantom$Receptor.ApprovedSymbol),]$mode <- "FGF")
  try(fantom[grep("TLR",fantom$Receptor.ApprovedSymbol),]$mode <- "TLR")
  try(fantom[grep("CCR",fantom$Receptor.ApprovedSymbol),]$mode <- "CC")
  try(fantom[grep("THBD",fantom$Receptor.ApprovedSymbol),]$mode <- 'Vasoactive')

  # Retracted
  fantom <- fantom[!fantom$Pair.Name == 'NAMPT-INSR',]
  # Uncategorized
  fantom[is.na(fantom$mode),]$mode <- 'UNCAT'
  # Make factor
  fantom$mode <- as.factor(fantom$mode)
  # Return
  return(fantom)
}

ncomms8866_human <- Annotate3(fantom = ncomms8866)

#ncomms8866_human[ncomms8866_human$mode == 'UNCAT',]$Pair.Name
save(ncomms8866_human, file = 'ncomms8866_human.rda')
