% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{mouse_to_human_orthologs}
\alias{mouse_to_human_orthologs}
\title{Table of Human-->Mouse orthologs for all human genes}
\format{
A data frame with nine variables:
\describe{
\item{\code{ensembl_gene_id}}{Human ensembl id}
\item{\code{hgnc_symbol}}{HGNC symbol}
\item{\code{entrezgene}}{Human NCBI gene ID}
\item{\code{mmusculus_homolog_ensembl_gene}}{Mouse ensembl id}
\item{\code{mmusculus_homolog_associated_gene_name}}{Mouse gene symbol}
\item{\code{mmusculus_homolog_orthology_type}}{The ortholog type (one2one, one2many, many2many)}
\item{\code{mmusculus_homolog_orthology_confidence}}{Orthology confidence}
}
}
\usage{
mouse_to_human_orthologs
}
\description{
A dataset containing mouse orthologs of human genes.
}
\keyword{datasets}
