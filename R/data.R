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
