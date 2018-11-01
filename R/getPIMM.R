#' add the 'percent immature mRNA' estimates to a SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @param exsce SingleCellExperiment with exon counts
#' @param intsce SingleCellExperiment with intron counts
#' @return SingleCellExperiment exsce with an additional 
#' colData column PIMM with the average value over all genes for
#' each cell of the quantity (intronCount/(intronCount+exonCount))
#' @export
addPIMM = function(exsce, intsce) {
 stopifnot(all.equal(colnames(exsce), colnames(intsce)))
 den = (assay(exsce)+assay(intsce))
 pct = ifelse(den != 0, assay(intsce)/den, NA)
 exsce$PIMM = apply(pct,2,mean, na.rm=TRUE)
 exsce
}
