% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EvenSplit.R
\name{EvenSplit}
\alias{EvenSplit}
\title{EvenSplit}
\usage{
EvenSplit(object, split.by)
}
\arguments{
\item{split.by}{Parameter by which to split the object}

\item{connectome}{A connectomic edgelist}
}
\description{
Splits a Seurat object using SplitObject, but makes sure that for each identity, an even number of cells are sampled.
}
