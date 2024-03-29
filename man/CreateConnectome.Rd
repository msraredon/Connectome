% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateConnectome.r
\name{CreateConnectome}
\alias{CreateConnectome}
\title{CreateConnectome}
\usage{
CreateConnectome(
  object,
  LR.database = "fantom5",
  species = NULL,
  include.putative = T,
  include.rejected = F,
  p.values = T,
  max.cells.per.ident = NULL,
  min.cells.per.ident = NULL,
  weight.definition.norm = "product",
  weight.definition.scale = "mean",
  custom.list = NULL,
  calculate.DOR = F,
  assay = "RNA",
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{LR.database}{Accepts either 'fantom5' or 'custom'. If custom, a dataframe must be provided to argument custom.list with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations.}

\item{species}{The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'}

\item{include.putative}{Default TRUE. Includes ligand-receptor pairs deemed putative in FANTOM5 database.}

\item{include.rejected}{Default FALSE. If TRUE, includes gene pairs labeled "EXCLUDED" in FANTOM5 database.  See ncomms8866 .rda file for qualifications for exclusion.}

\item{p.values}{Default TRUE. Runs a Wilcoxon Rank test to calculate adjusted p-value for ligand and receptor expression within the input object. Change to FALSE for decreased run-time.}

\item{max.cells.per.ident}{Default NULL. If a value is input, input object will be downsampled to requested number of cells per identity. This can greatly improve run-time.}

\item{min.cells.per.ident}{Default NULL. If a value is input, only cell populations meeting this threshold will be included in network analysis. Can limit high-variation effects from small clusters.}

\item{weight.definition.norm}{Method of edgeweight definition from normalized data slot. Either 'sum','mean',or 'product'. Defaults to 'product'. 'Sum' adds values from sending and receiving clusters, 'mean' averages them, and 'product' multiplies them.}

\item{weight.definition.scale}{Method of edgeweight definition from scaled data slot. Either 'sum','mean',or 'product'. Defaults to 'mean'. 'Sum' adds values from sending and receiving clusters, 'mean' averages them, and 'product' multiplies them.}

\item{custom.list}{Optional. A dataframe for custom mapping, with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations. If modal categorizations are unknown, fill with 'UNCAT' or similar placeholder.}

\item{calculate.DOR}{Whether to include log normalized Diagnostic Odds Ratio in the output.  Default FALSE.}

\item{assay}{Which assay to employ to create connectomic mapping. Default = 'RNA', also accepts 'SCT' and others}
}
\description{
Creates an unfiltered connectomic edgelist from a Seurat 3.0 object.  Allows orphan ligand and receptors and fully non-expressed pairs.
}
