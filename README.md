# connectome v0.2.0

## scRNAseq connectomics

`Connectome` is a toolkit to explore cell-cell connectivity patterns based on ligand and receptor data in heterogeneous single-cell datasets. See Raredon MSB et al (2019) <doi:10.1126/sciadv.aaw3851> for more details.

Currently capable of creating mappings for human, mouse, rat, and pig, against the FANTOM5 ligand-receptor data found in Ramilowksi JA et al (2015) <doi:10.1038/ncomms8866>, or against any user-provided list of paired ligand-receptor interactions.

<p align="center">
<img src="https://github.com/msraredon/connectome/blob/develop/Big%20Connectome.png" alt="BigC"
	title="BigC" width="150" height="150" align="center" />
</p>

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
