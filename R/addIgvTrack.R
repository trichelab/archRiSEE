#' wrapper to add a weights track to an igvR instance 
#' 
#' @param   x   SingleCellExperiment with weights in mcols(rowRanges)[[trk]]
#' @param   trk name of the column in mcols(rowRanges(x)) with weights to plot
#' @param   igv igvR object with the appropriate genome to match x 
#' @param   co  color for the track (violet)
#' @param   a   autoscale the track? (TRUE)
#' @param   ... additional parameters passed to igvR::displayTrack
#' 
#' @return  none, called for side effects 
#' 
#' @import  igvR
#' 
addIgvTrack <- function(x, trk, igv, co="violet", a=TRUE, ...) {

  rr <- rowRanges(sce)
  setGenome(igv, unique(genome(rr)))
  mcols(rr) <- mcols(rr)[, trk, drop=FALSE]
  names(mcols(rr)) <- "score"
  displayTrack(igv, GRangesQuantitativeTrack(trk, rr, color=co, autoscale=a))

}
