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
PcmpSels = function(obj) new("PcmpSels", cellSets=obj$cells,
   geneTable=obj$limmaTab)

setMethod("show", "PcmpSels", function(object) {
 cat(sprintf("PcmpSels instance, collecting %d selections from pcmpApp.\n",
    length(object@cellSets)))
})
