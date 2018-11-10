#' organize record of selections in pcmp
#' @importFrom graphics plot points
#' @importFrom stats as.formula biplot prcomp
#' @name PcmpSels
#' @rdname PcmpSels
#' @aliases PcmpSels-class
#' @slot cellSets list() of cell identifiers for each selection
#' @slot geneTable a single data.frame with a variable `sel` indicating
#' the selections to which rows belong
#' @exportClass PcmpSels
setClass("PcmpSels", representation(cellSets="list",
    geneTable="data.frame"))
#' construct pcmp selections object
#' @param obj a list with elements 'cells' and 'limmaTab'
#' @export
PcmpSels = function(obj) new("PcmpSels", cellSets=obj$cells,
   geneTable=obj$limmaTab)

setMethod("show", "PcmpSels", function(object) {
 cat(sprintf("PcmpSels instance, collecting %d selections from pcmpApp.\n",
    length(object@cellSets)))
 cat("  use cellSets(), geneTable() or replay() to extract information.\n")
})

#' extract list of cells selected
#' @param x instance of PcmpSels
#' @export cellSets
setGeneric("cellSets", function(x)standardGeneric("cellSets"))

#' @rdname cellSets
setMethod("cellSets", "PcmpSels", function(x) x@cellSets)

#' extract limma topTable
#' @param x instance of PcmpSels
#' @export geneTable
setGeneric("geneTable", function(x)standardGeneric("geneTable"))

#' @rdname geneTable
setMethod("geneTable", "PcmpSels", function(x) x@geneTable)

#' replay a selection by plotting the selected cells
#' @param sce a SingleCellExperiment instance to which pcmpApp was applied
#' @param sels an instance of PcmpSels
#' @param whichProj character(1) name of reducedDims list element to use
#' @param dim1 numeric(1) x axis of replay
#' @param dim2 numeric(1) y axis of replay
#' @examples
#' replay(pcmp::sce300xx, pcmp::acc4vis2, "TS", 3, 4)
#' @export
replay = function (sce, sels, whichProj, dim1=1, dim2=2)
{
    stopifnot(is(sce, "SummarizedExperiment"), is(sels, "PcmpSels"))
    stopifnot(whichProj %in% names(reducedDims(sce)))
    rddat = reducedDims(sce)[[whichProj]]
    nsels = length(cs <- cellSets(sels))
    plot(rddat[, dim1], rddat[, dim2], col = "gray", xlab = paste0(whichProj, 
        "(", dim1, ")"), ylab = paste0(whichProj, "(", dim2, 
        ")"))
    for (i in 1:nsels) points(rddat[cs[[i]], dim1], rddat[cs[[i]], 
        dim2], col = i + 1)
    legx = min(rddat[,dim1], na.rm=TRUE)
    legy = max(rddat[,dim2], na.rm=TRUE)
    legend(1.05*legx, .95*legy, legend=paste0("sel", 1:nsels),
      col=2:(nsels+1), pch=19)
}

replayOld = function(sce, sels, whichProj, dim1, dim2) {
 stopifnot(is(sce, "SummarizedExperiment"), is(sels, "PcmpSels"))
 stopifnot(whichProj %in% names(reducedDims(sce)))
 rddat = reducedDims(sce)[[whichProj]]
 nsels = length(cs <- cellSets(sels))
 plot(rddat[,dim1], rddat[,dim2], col="gray", xlab=paste0(whichProj, "(",dim1,")"),
       ylab=paste0(whichProj, "(", dim2, ")"))
 for (i in 1:nsels)
   points(rddat[cs[[i]],dim1], rddat[cs[[i]],dim2], col=i+1)
 }


#' produce biplot of top genes for a selection
#' @param sce a SingleCellExperiment instance to which pcmpApp was applied
#' @param sels an instance of PcmpSels
#' @param which a numeric(1) telling which of the 
#' selections is used to form a biplot.  The 
#' order used is that of `geneSets(sels)`
#' @param ntopgenes numeric(1) number of 
#' (most differentially-expressed, selection vs all others) 
#' genes to use in PCA for selected samples; the
#' assay data are transformed as log(assay+1) before PCA is carried out
#' @param cex passed to biplot.prcomp
#' @param \dots passed to biplot.prcomp
#' @examples
#' op = par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' for (i in 1:4) biplotSel(pcmp::sce300xx, pcmp::acc4vis2, which=i)
#' par(op)
#' @export
biplotSel = function (sce, sels, which = 1, ntopgenes=6, cex=c(2,1), ...) 
{
    cs = cellSets(sels)
    gt = geneTable(sels)
    gts = split(gt, gt[, "selnum"])
    gns = lapply(gts, function(x) x$featid)
    curpc = prcomp(t(log(assay(sce[gns[[which]][1:ntopgenes], cs[[which]]])+1)))
    biplot(curpc, xlabs=rep(".", length(cs[[which]])), expand=.8, cex=cex, ...)
}
