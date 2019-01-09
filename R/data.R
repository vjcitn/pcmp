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

#' collection of 6 selections from sce300xx
#' @format PcmpSels instance
#' @examples
#' data(acc4vis2)
#' data(sce300xx)
#' rd = reducedDims(sce300xx)
#' plot(rd$TS[,3], rd$TS[,4])
#' inds1 = acc4vis2@cellSets[[1]]
#' points(rd$TS[ inds1, 3], rd$TS[ inds1, 4], pch=19, col="red")
"acc4vis2"

#' collection of selections from sce300xx, corresponding to vignette displays
#' @format PcmpSels instance
#' @examples
#' data(vigAccum)
#' head(geneTable(vigAccum),3)
"vigAccum"

#' annotation of immgen mouse immune cell resource
#' @format data.frame
#' @source \url{https://gist.github.com/nachocab/3d9f374e0ade031c475a}
#' @examples
#' head(immgenAnno)
"immgenAnno"
