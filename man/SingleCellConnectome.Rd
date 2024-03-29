% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SingleCellConnectome.R
\name{SingleCellConnectome}
\alias{SingleCellConnectome}
\title{SingleCellConnectome}
\usage{
SingleCellConnectome(
  object,
  LR.database = "fantom5",
  species = NULL,
  include.putative = T,
  include.rejected = F,
  max.cells.per.ident = NULL,
  min.cells.per.ident = NULL,
  slot.use = "data",
  weight.definition = "product",
  custom.list = NULL,
  ...
)
}
\arguments{
\item{object}{A Seurat object}

\item{LR.database}{Accepts either 'fantom5' or 'custom'. If custom, a dataframe must be provided to argument custom.list with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations.}

\item{species}{The species of the object that is being processed.  Only required if LR.database = 'fantom5', and allows 'human','mouse','rat', or 'pig'}

\item{include.putative}{Default TRUE. Includes ligand-receptor pairs deemed putative in FANTOM5 database.}

\item{include.rejected}{Default FALSE. If TRUE, includes gene pairs labeled "EXCLUDED" in FANTOM5 database.  See ncomms8866 .rda file for qualifications for exclusion.}

\item{max.cells.per.ident}{Default NULL. If a value is input, input object will be downsampled to requested number of cells per identity. This can greatly improve run-time.}

\item{min.cells.per.ident}{Default NULL. If a value is input, only cell populations meeting this threshold will be included in network analysis. Can limit high-variation effects from small clusters.}

\item{custom.list}{Optional. A dataframe for custom mapping, with the first column equal to ligands, second column equal to associated receptors, and third column equal to desired modal categorizations. If modal categorizations are unknown, fill with 'UNCAT' or similar placeholder.}
}
\description{
Takes a Seurat object and creates a single-cell connectome (rows are ligand-receptor mechanisms and columns are cell-cell vectors)
}
