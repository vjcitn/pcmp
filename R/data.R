#' demonstration SingleCellExperiment
#' @importFrom utils data
#' @docType data
#' @format SingleCellExperiment
#' @source see metadata component
#' @note Developed by stratified random sampling, with strata
#' defined by donor (3 levels) and brain region (2 levels, ACC and VIS).
#' Retains 300 cells chosen randomly within each stratum.
#' @examples
#' suppressPackageStartupMessages({
#'  requireNamespace("SingleCellExperiment")
#' })
#' pcmp::sce300xx
"sce300xx"

#' collection of 6 selections from reduced sce300xx
#' @note gene set used for limmaTabs formed by requiring MAD over samples > 0
#' @format SingleCellExperiment instance
#' @examples
#' data(acc4vis2)
#' plotSelMap(acc4vis2, "TSNE2")
"acc4vis2"

#' collection of selections from sce300xx, corresponding to vignette displays
#' @format SingleCellExperiment instance, produced by pcmpApp
#' @examples
#' data(vigAccum)
#' head(metadata(vigAccum)$limmaTabs,3)
"vigAccum"

#' annotation of immgen mouse immune cell resource
#' @format data.frame
#' @source \url{https://gist.github.com/nachocab/3d9f374e0ade031c475a}
#' @examples
#' head(immgenAnno)
"immgenAnno"

#' Antibody-capture "protein" data for the 10x PBMC 10k experiment
#' @format data.frame
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_10k_protein_v3} 
#' imported using DropletUtils, subset to the antibody feature measurements,
#' reducedDims added using addProjections/defaultProjectors modified
#' to use svd directly with this small featureset, and clustered
#' using kmeans for 9 clusters (cluster indicator added to colData).
#' @examples
#' Ab10k
"Ab10k"
