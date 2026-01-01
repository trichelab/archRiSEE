#' wrapper to add a weights track to an igvR instance 
#' 
#' @param   x   SingleCellExperiment with weights in mcols(x)$trk
#' @param   trk name of the column in mcols(x) with weights to plot
#' @param   igv an igvR object with an appropriate genome to match genome(x)
#' @param   co  color for the track (violet)
#' @param   a   autoscale the track? (TRUE)
#' @param   ... additional parameters passed to igvR::GRangesQuantitativeTrack
#' 
#' @return  nothing useful (called for side effects)
#' 
#' @examples
#' library(igvR)
#' igv <- igvR()
#' setGenome(igv, "hg38")
#' mcols(SCE)$aTrack <- rpois(nrow(SCE), lambda=50)
#' addIgvTrack(SCE, "aTrack", igv) 
#'
#' @seealso igvR::GRangesQuantitativeTrack
#' @seealso igvR::displayTrack
#' @seealso igvR::setGenome
#' @seealso igvR::igvR
#'
#' @import  igvR
#' 
#' @export
#'
addIgvTrack <- function(x, trk, igv, co="violet", a=TRUE, ...) {

  stopifnot(is(igv, "igvR"))
  g <- unique(genome(x))
  if (is.na(g)) {
    warning("No genome currently set for x; fix this to avoid major mishaps.")
  } else { 
    stopifnot(g %in% getSupportedGenomes(igv))
    message('Please be sure to setGenome(igv, "', g, '") before adding tracks.')
  }
  rr <- .makeIgvTrackFromRowRanges(x, trk)  
  displayTrack(igv, 
               GRangesQuantitativeTrack(trk, rr, color=co, autoscale=a, ...))

}


# helper fn
.makeIgvTrackFromRowRanges <- function(x, trk) { 

  rr <- rowRanges(x)
  mcols(rr) <- mcols(rr)[, trk, drop=FALSE]
  names(mcols(rr)) <- "score"
  return(rr)

}
