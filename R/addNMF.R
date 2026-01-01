#' add an NMF decomposition to the reducedDims and rowRanges of an SCE
#'
#' By default, the weights for each factor are copied to mcols(rowRanges(x)).
#'
#' @param x             a SingleCellExperiment, usually from archRtoSCE
#' @param k             rank for the NMF decomposition on logCounts (30) 
#' @param colDat        add NMF scores to colData(x)? (TRUE)
#' @param rowDat        add NMF weights to mcols(rowRanges(x))? (TRUE)
#' @param ...           additional arguments to pass to RcppML::nmf()
#'
#' @return              SingleCellExperiment with reducedDim(x, "NMF")
#'
#' @details             the mcols(rowRanges(x)) weights are for igvR plotting 
#'
#' @seealso             archRtoSCE
#' @seealso             iSEEarchR
#' @seealso             addIgvTrack
#'
#' @import RcppML
#' @import scuttle
#' @import SingleCellExperiment
#'
#' @export
#'
addNMF <- function(x, k=30, colDat=TRUE, rowDat=TRUE, ...) { 

  if (! "logcounts" %in% assayNames(x)) x <- logNormCounts(x) 
  
  orig <- options()[["RcppML.verbose"]]
  message("Fitting rank-", k, " NMF model on logCounts(x)...")
  options("RcppML.verbose" = TRUE)
  metadata(x)$NMF <- RcppML::nmf(logcounts(x), k=k, ...)
  message("Copying NMF hat matrix to reducedDim(x, 'NMF')...")
  reducedDim(x, "NMF") <- t(metadata(x)$NMF@h)
  NMFdims <- ncol(reducedDim(x, "NMF"))
  names(reducedDim(x, "NMF")) <- paste0("NMF", seq_len(NMFdims))
  if (colDat) { 
    message("Copying NMF hat columns to colData() for iSEE visualization...")
    for (i in colnames(reducedDim(x, "NMF"))) {
      message("Adding ", i, " as colData(x)$", toupper(i))
      colData(x)[, toupper(i)] <- reducedDim(x, "NMF")[, i]
    }
  }
  if (rowDat) {
    message("Copying NMF weights to mcols(rowRanges(x)) for igvR plotting...")
    for (i in colnames(metadata(x)$NMF@w)) { 
      message("Adding NMF@w[, ", i, "] as mcols(rowRanges(x))$", toupper(i))
      mcols(x)[, toupper(i)] <- metadata(x)$NMF@w[, i]
    }
  }
  options("RcppML.verbose" = orig)
  return(x) 

}
