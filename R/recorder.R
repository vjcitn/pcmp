#' organize record of selections in pcmp
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

#' @export cellSets
setGeneric("cellSets", function(x)standardGeneric("cellSets"))
setMethod("cellSets", "PcmpSels", function(x) x@cellSets)
#' @export geneTable
setGeneric("geneTable", function(x)standardGeneric("geneTable"))
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
replay = function(sce, sels, whichProj, dim1, dim2) {
 stopifnot(is(sce, "SummarizedExperiment"), is(sels, "PcmpSels"))
 stopifnot(whichProj %in% names(reducedDims(sce)))
 rddat = reducedDims(sce)[[whichProj]]
 nsels = length(cs <- cellSets(sels))
 plot(rddat[,dim1], rddat[,dim2], col="gray", xlab=paste0(whichProj, "(",dim1,")"),
       ylab=paste0(whichProj, "(", dim2, ")"))
 for (i in 1:nsels)
   points(rddat[cs[[i]],dim1], rddat[cs[[i]],dim2], col=i+1)
 }


#' produce biplots of top genes for selections
#' @param sce a SingleCellExperiment instance to which pcmpApp was applied
#' @param sels an instance of PcmpSels
#' @param which a numeric(4) collection telling which of the selections occupy a 2 x 2 biplot array
#' @param ntopgenes numeric(1) number of (most differentially-expressed, selection vs all others) genes to use in PCA for selected samples
#' @examples
#' op = par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' for (i in 1:4) biplots(pcmp::sce300xx, pcmp::acc4vis2, which=i)
#' par(op)
#' @export
biplots = function (sce, sels, which = 1, ntopgenes=6) 
{
#    opar = par(no.readonly = TRUE)
#    on.exit(par(opar))
#    par(mfrow = c(2, 2))
    cs = cellSets(sels)
    gt = geneTable(sels)
    gts = split(gt, gt[, "selnum"])
    gns = lapply(gts, rownames)
    curpc = prcomp(t(log(assay(sce[gns[[which]][1:ntopgenes], cs[[which]]])+1)))
    biplot(curpc, xlabs=rep(".", length(cs[[which]])), expand=.8)
}
