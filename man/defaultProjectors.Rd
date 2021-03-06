% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infra.R
\name{defaultProjectors}
\alias{defaultProjectors}
\title{schematize a sequence of projection methods}
\usage{
defaultProjectors(warmdim = 15)
}
\arguments{
\item{warmdim}{numeric(1) defaults to 15 -- defines the dimension
reduction using irlba()$u for 'warming up' prior to UMAP or t-SNE runs}
}
\value{
a list with two components, 'projectors'
and 'retrievers'.  'projectors' is a named list of
 functions with parameters 'x' and 'ncomp',
where the intended usage has 'x' inheriting from SingleCellExperiment, 
and 'ncomp' numeric defining the number of projection components
to be retained.  'retrievers' is a named list of functions
with parameter 'x', which extracts the ncomp x nsamp
projection for each method.
}
\description{
schematize a sequence of projection methods
}
\note{
The default sequence is irlba::irlba, umap::umap,
Rtsne::Rtsne; see the note on `warmdim` parameter.
Prior to pcmp version 0.3.0, TSNE and UMAP projections were
computed 'in one shot' with high target numbers of dimensions.
Beginning in 0.3.0, TSNE and UMAP projections are computed
targeting 2 and 3 dimensions separately for TSNE, and 2-4
dimensions for UMAP.
}
\examples{
defaultProjectors
}
