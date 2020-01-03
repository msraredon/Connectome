# connectome
scRNAseq connectomics

A R package containing functions and files ready for beta distribution to collaborators. Intended to cover all the basics of connectomic analysis. Will evolve with time, and will be curated to continuously update.

A good working order is:

1) Generate Seurat object
2) Convert to Homologues
3) Generate Connectome
4) Calculate and add P-values to Connectome
5) Annotate Connectome with signaling modes

6) Filter to nodes of interest, p-value thresh, percent expression, if desired
7) Affinity + Network comparisons / visualization
8) Hive plot visualizations
9) Nearest Neighbor, Modal NN, and Modal Dot Plot visualizations
