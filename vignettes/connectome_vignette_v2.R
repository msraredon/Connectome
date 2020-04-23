# Load data
load("~/Box Sync/Science Response/Base objects for computational work (Use Human GEN2 for revision work)/objects/panc8_cca 2018-08-04 .Robj")

# Attach packages
require(Seurat)
require(SeuratData)
require(connectome)

# Load data
InstallData('panc8')
data('panc8')
table(Idents(panc8))
Idents(panc8) <- panc8[['celltype']]
table(Idents(panc8))

# Scale all genes (can also, if necessary, only scale ligands and receptors)
panc8 <- ScaleData(panc8,features = rownames(panc8))

# Make connectome
panc8.con <- CreateConnectome(panc8,species = 'human',p.values = F,calculate.DOR = F)

# Filter to edges of interest (many options here)
panc8.con2 <- FilterConnectome(panc8.con,min.pct = 0.1,min.z = 0.25,remove.na = T)

# Network plot
NetworkPlot(panc8.con2,'Vegfa',min.pct = 0.75,weight.attribute = 'weight_sc')

# ModalDotPlot (all modes)
ModalDotPlot(panc8.con2,modes.include = NULL,min.z = NULL,weight.attribute = 'weight_sc')

# ModalDotPlot (limited modes)
ModalDotPlot(panc8.con2,modes.include = c('VEGF','WNT','Semaphorins','NOTCH','FGF'),weight.attribute = 'weight_sc',min.z = 0)

# CellCellScatter
CellCellScatter(panc8.con2,sources.include = 'ATI',targets.include = 'Fib_Col13a1+',
                label.threshold = 1,
                weight.attribute = 'weight_sc',min.pct = 0.25,min.z = 0)

# SignalScatter
SignalScatter(panc8.con2, features = 'Wnt3a',label.threshold = 1,weight.attribute = 'weight_sc')

# Four kinds of CircosPlots
test <- FilterConnectome(panc8.con2,sources.include = c('ATI','ATII'),targets.include = c('Fib_Col13a1+','SMCS','Fib_Col14a1+'))
test <- data.frame(test %>% group_by(vector) %>% top_n(5,weight_sc))
CircosPlot(test,weight.attribute = 'weight_norm')
CircosPlot(test,weight.attribute = 'weight_sc')
CircosPlot(test,weight.attribute = 'weight_norm',balanced.edges = F)
CircosPlot(test,weight.attribute = 'weight_sc',balanced.edges = F)

CircosPlot(panc8.con2,min.z = 1,
           targets.include = 'EC_cap')

#### Differential demo ####

# Load in Seupanc8 Differential Demo data
require(Seupanc8Data)
InstallData('ifnb')
data('ifnb')
table(Idents(ifnb))
Idents(ifnb) <- ifnb[['seupanc8_annotations']]
table(Idents(ifnb))

# Split by condition
ifnb.list <- SplitObject(ifnb,split.by = 'stim')

#  Make connectomes
ifnb.con.list <- list()
for (i in 1:length(ifnb.list)){
  ifnb.list[[i]] <- ScaleData(ifnb.list[[i]],features = rownames(ifnb.list[[i]]))
  ifnb.con.list[[i]] <- CreateConnectome(ifnb.list[[i]],species = 'human',p.values = F)
}
names(ifnb.con.list) <- names(ifnb.list)

# Make differential connectome
diff <- DifferentialConnectome(ifnb.con.list[[1]],ifnb.con.list[[2]])

# Initial Scoring Plot
pdf(file = 'Scoring Plot Demo.pdf',width = 50,height=8)
DifferentialScoringPlot(diff,min.score = 10,min.pct = 0.1,infinity.to.max = T)
dev.off()

# Differential Circos Plot (edgeweights here are the perturbation score)
# All edges meeting thresholds:
CircosDiff(diff,min.score = 10,min.pct = 0.1,infinity.to.max = T,lab.cex = 0.4)
# Just specific cell interactions:
CircosDiff(diff,min.score = 10,min.pct = 0.1,infinity.to.max = T,lab.cex = 0.4,
           sources.include = c('pDC','CD8 T','B'),targets.include = c('CD16 Mono','CD14 Mono'))
# Differential signaling in a specific cellular 'niche' due to perturbation:
CircosDiff(diff,min.score = 5,min.pct = 0.1,infinity.to.max = T,lab.cex = 0.4,
           targets.include = c('CD14 Mono'))
# Can look at differential signaling hitting a single receptor:
CircosDiff(diff,min.pct = 0.1,infinity.to.max = T,lab.cex = 0.4,
           targets.include = c('CD14 Mono'),features = c('CCR5'))
