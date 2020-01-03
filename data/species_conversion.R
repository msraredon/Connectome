setwd("~/Documents/GitHub/connectome/data")
load("~/Documents/GitHub/connectome/data/translator.rda")
load("~/Documents/GitHub/connectome/data/ncomms8866_human.rda")

Convert_GOI_to_Homologues <- function(GOI,old.species,new.species){
  hom <- subset(translator,translator$query.species %in% old.species & translator$match.species == new.species)
  hom.use <- names(which(table(hom$query.gene.symbol) == 1)) # 1:1 condition A
  hom.use2 <- names(which(table(hom$match.gene.symbol) == 1)) # 1:1 condition B
  hom.list <- hom[hom$query.gene.symbol %in% hom.use & hom$match.gene.symbol %in% hom.use2,]
  #Raw data gene names
  old.genes <- GOI
  for (i in 1:length(old.genes)){
    if (old.genes[i] %in% hom.list$query.gene.symbol){
      old.genes[i] <- hom.list[hom.list$query.gene.symbol == old.genes[i],]$match.gene.symbol
    }
  }
  if (new.species == 'pig'){
  for (i in 1:length(old.genes)){
    if (toupper(old.genes[i]) %in% hom[hom$query.gene.symbol == old.genes[i],]$match.gene.symbol &
        old.genes[i] != toupper(old.genes[i])){
      index <-  hom[hom$query.gene.symbol == old.genes[i],]
      old.genes[i] <- unique(index[index$match.gene.symbol == toupper(old.genes[i]),]$match.gene.symbol)
    }
  }
  }
  else{
    for (i in 1:length(old.genes)){
      if (paste(substring(old.genes[i], 1,1), tolower(substring(old.genes[i], 2)), sep="", collapse=" ") %in% hom[hom$query.gene.symbol == old.genes[i],]$match.gene.symbol &
          old.genes[i] != paste(substring(old.genes[i], 1,1), tolower(substring(old.genes[i], 2)), sep="", collapse=" ")){
        index <-  hom[hom$query.gene.symbol == old.genes[i],]
        old.genes[i] <- unique(index[index$match.gene.symbol == paste(substring(old.genes[i], 1,1), tolower(substring(old.genes[i], 2)), sep="", collapse=" "),]$match.gene.symbol)
      }
    }  
  }
  new.genes <- old.genes
  return(new.genes)
}
 # Test the above function
GOI <- c('VEGFA','WNT5A','PTPRC','RLN1','RLN3')
Convert_GOI_to_Homologues(GOI,old.species = 'human',new.species = 'rat')
Convert_GOI_to_Homologues(GOI,old.species = 'human',new.species = 'mouse')
Convert_GOI_to_Homologues(GOI,old.species = 'human',new.species = 'pig')

# Perform conversions
ncomms8866_rat <- ncomms8866_human
ncomms8866_rat$Ligand.ApprovedSymbol <- Convert_GOI_to_Homologues(ncomms8866_rat$Ligand.ApprovedSymbol,old.species = 'human',new.species = 'rat')
ncomms8866_rat$Receptor.ApprovedSymbol <- Convert_GOI_to_Homologues(ncomms8866_rat$Receptor.ApprovedSymbol,old.species = 'human',new.species = 'rat')
save(ncomms8866_rat,file = 'ncomms8866_rat.rda')


ncomms8866_mouse <- ncomms8866_human
ncomms8866_mouse$Ligand.ApprovedSymbol <- Convert_GOI_to_Homologues(ncomms8866_mouse$Ligand.ApprovedSymbol,old.species = 'human',new.species = 'mouse')
ncomms8866_mouse$Receptor.ApprovedSymbol <- Convert_GOI_to_Homologues(ncomms8866_mouse$Receptor.ApprovedSymbol,old.species = 'human',new.species = 'mouse')
save(ncomms8866_mouse,file = 'ncomms8866_mouse.rda')

ncomms8866_pig <- ncomms8866_human
ncomms8866_pig$Ligand.ApprovedSymbol <- Convert_GOI_to_Homologues(ncomms8866_pig$Ligand.ApprovedSymbol,old.species = 'human',new.species = 'pig')
ncomms8866_pig$Receptor.ApprovedSymbol <- Convert_GOI_to_Homologues(ncomms8866_pig$Receptor.ApprovedSymbol,old.species = 'human',new.species = 'pig')
save(ncomms8866_pig,file = 'ncomms8866_pig.rda')

