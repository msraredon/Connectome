# Connectome v0.2.2

## scRNAseq connectomics

`Connectome` is an R toolkit to explore cell-cell connectivity patterns based on ligand and receptor data in heterogeneous single-cell datasets. It is designed to work with Seurat from Satija Lab.

This software compiles and extends the methods described in Raredon MSB et al (2019) <doi:10.1126/sciadv.aaw3851>.

Currently capable of creating mappings for human, mouse, rat, and pig, against the FANTOM5 ligand-receptor data found in Ramilowksi JA et al (2015) <doi:10.1038/ncomms8866>, or against any user-provided list of paired ligand-receptor interactions.

<p align="center">
<img src="/man/figures/Big_Connectome.png" alt="BigConnectome"
	title="BigConnectome" width="300" height="300" />
</p>

## Installation
To install `Connectome` in R, you may run:

```
library(devtools)
install_github('msraredon/connectome', ref = 'master')
```

## Functions to analyze a single tissue system:

`CreateConnectome` takes as input a Seurat 3.0 object and uses the active identity slot to define nodes for network analysis. The output of this function is an edgelist connecting pairs of nodes via specific ligand-receptor mechanisms, with many edge attributes allowing downstream processing and filtration.

`FilterConnectome` is a streamlined way of reducing the above edgelist (which is generally large) to a smaller subset of edges more likely to be of biological and statistical interest.

`NetworkPlot` provides a simple way to visualize networks of interest.  It is a wrapper for the R package `igraph`, and allows filtration on all arguments which can be passed to `FilterConnectome`.

`Centrality` creates a paired centrality plot allowing identification of cell types dominating production or reception of specific modes of signaling. Allows filtration on all arguments which can be passed to `FilterConnectome`.

`SignalScatter` aids in identification of top cell-cell signaling vectors within a network of interest.

`CellCellScatter` aids in identification of top signaling mechanisms between a specified source and target cell of interest.

`CircosPlot` replaces the hive plots displayed in the original manuscript.  CircosPlot is a versatile plotting function which includes many adaptable parameters and yields easy-to interpret quantitative graphs using the R package `circlize`. Allows filtration on all arguments which can be passed to `FilterConnectome`.

## Functions to compare cell-cell signaling across two tissue systems:

`DifferentialConnectome` Allows direct quantitative comparison of two connectomes.

`DifferentialScoringPlot` Provides a comprehensive heatmap view of differences of interest between two connectomes.

`CircosDiff` Creates a CircosPlot of a differential connectome, leveraging an always positive 'perturbation score'

## Functions to compare cell-cell signaling across two or more tissue systems:

`CompareCentrality` Takes any list of connectomes and compares sending- and receiving- centrality, side-by-side, for a given network subset.

## Reference
Please cite Raredon, Micha Sam Brickman, et al. "Single-cell connectomic analysis of adult mammalian lungs." Science Advances 5.12 (2019): eaaw3851. <doi:10.1126/sciadv.aaw3851>
