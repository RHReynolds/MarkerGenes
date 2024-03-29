% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{run_ewce_controlled}
\alias{run_ewce_controlled}
\title{Deprecated function: Run EWCE with control for GC content and transcript
length.}
\usage{
run_ewce_controlled(
  list_of_genes,
  ...,
  celltypeLevel = c(1, 2, "both"),
  reps,
  genelistSpecies = c("human"),
  sctSpecies = c("mouse", "human"),
  mouse_to_human = NULL
)
}
\arguments{
\item{list_of_genes}{Named list containing vectors with gene IDs
corresponding to each gene list to be tested.}

\item{...}{Names of input ctds. An individual ctd, should have been generated
using \code{generate.celltype.data} from the EWCE package, and contains
mean expression and specificity matrices from a scRNAseq study.}

\item{celltypeLevel}{Choose whether to run only with level 1 cell types,
level 2 cell types or both. 1 for only level 1, or "both" for both.}

\item{reps}{Number of repeats that should be performed in bootstrapping.}

\item{genelistSpecies}{Has to be human, as geneSizeControl is set to TRUE,
and this requires extracting GC content and transcript length from the
human ENSEMBL.}

\item{sctSpecies}{Either 'mouse' or 'human' depending on the ctd datasets
being used. Species must be the same across all input ctds.}

\item{mouse_to_human}{Default is to use the default dataframe supplied by
within MarkerGenes packages, which contains all Human->Mouse orthologs for
all human genes. If user wishes to supply their own dataframe, must contain
the columns, 'hgnc_symbol' and 'MGI_symbol'.}
}
\value{
A dataframe of EWCE results for each of the gene lists in each study.
}
\description{
This function is deprecated on account of major changes to EWCE (see EWCE
news: \url{https://nathanskene.github.io/EWCE/news/index.html}), which also
include development of
[MultiEWCE](\url{https://github.com/neurogenomics/MultiEWCE/}), which permits
running multiple gene lists in parallel. Perform EWCE bootstrapping across
several gene lists in more than one specificity matrix (aka. ctd).
}
\keyword{internal}
