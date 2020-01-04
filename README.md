# connectome v0.1.0

scRNAseq connectomics

A toolkit to explore cell-cell connectivity patterns based on ligand and receptor data in heterogeneous single-cell datasets. See Raredon MSB et al (2019) <doi:10.1126/sciadv.aaw3851> for more details.

Currently capable of creating mappings for human, mouse, rat, and pig, against the FANTOM5 ligand-receptor data found in Ramilowksi JA et al (2015) <doi:10.1038/ncomms8866>.

CreateConnectome takes as input a Seurat 3.0 object and uses the active identity slot to define nodes for network analysis. The output of this function is an edgelist connecting pairs of nodes via specific ligand-receptor mechanisms, with many edge attributes allowing downstream processing and filtration.

FilterConnectome is a streamlined way of reducing the above edgelist (which is generally large) to a smaller subset of edges more likely to be of biological interest.

NetworkPlot provides a simple way to visualize networks of interest.  Wrapper for igraph package.

SignalScatter is a simple ggplot wrapper allowing identification of top cell-cell signaling vectors involving a specific gene of interest.

HiveOrganization will provide output to make hive plots via http://wodaklab.org/hivegraph/graph
