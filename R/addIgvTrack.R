#' wrapper to add a weights track to an igvR instance 
#' 
#' @param   x   SingleCellExperiment with weights in mcols(x)$trk
#' @param   trk name of the column in mcols(x) with weights to plot
#' @param   igv an igvR object (genome will be set to match unique(genome(x)))
#' @param   co  color for the track (violet)
#' @param   a   autoscale the track? (TRUE)
#' @param   ... additional parameters passed to igvR::GRangesQuantitativeTrack
#' 
#' @return  none, called for side effects 
#' 
#' @examples
#' library(igvR)
#' igv <- igvR()
#' mcols(SCE)$aTrack <- rpois(nrow(SCE), lambda=50)
#' addIgvTrack(SCE, "aTrack", igv) 
#'
#' @seealso igvR::GRangesQuantitativeTrack
#' @seealso igvR::displayTrack
#' @seealso igvR::igvR
#'
#' @import  igvR
#' 
addIgvTrack <- function(x, trk, igv, co="violet", a=TRUE, ...) {

  stopifnot(is(igv, "igvR"))
  rr <- rowRanges(x)
  setGenome(igv, unique(genome(rr)))
  mcols(rr) <- mcols(rr)[, trk, drop=FALSE]
  names(mcols(rr)) <- "score"
  displayTrack(igv, 
               GRangesQuantitativeTrack(trk, rr, color=co, autoscale=a, ...))

}
