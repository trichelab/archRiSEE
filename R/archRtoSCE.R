#' turn an ArchR project into a SingleCellExperiment 
#'
#' The default tile size is 500bp, just like ArchR's default. 
#'
#' @param proj          an ArchRproject 
#' @param tile          tile size (500)
#' @param addNMF        add an NMF decomposition of logCounts? (FALSE)
#' @param k             rank for the NMF decomposition on logCounts (30) 
#' @param colDat        add LSI and/or NMF scores to colData(SCE)? (FALSE)
#' @param ...           additional arguments to pass to ArchR::addTileMatrix
#'
#' @return              a SingleCellExperiment 
#'
#' @details presumes LSI and UMAP. you will almost certainly want to tweak this.
#'
#' @import ArchR
#' @import scater
#' @import RcppML
#' @import GenomicRanges
#' @import SingleCellExperiment
#'
#' @export
#'
archRtoSCE <- function(proj, tile=500, addNMF=FALSE, k=30, colDat=FALSE, ...) { 

  message("Adding TileMatrix to ArchR project (this may take a while)...")
  proj <- addTileMatrix(proj, force=TRUE, binarize=FALSE, tileSize=tile, ...)
  message("Converting to SingleCellExperiment...")
  SCE <- as(getMatrixFromProject(proj, "TileMatrix"), "SingleCellExperiment") 
  rowData(SCE)$assay<- "FragmentCounts" # for binding to any other altExps
  message("Log-normalizing fragment counts...") 
  SCE <- logNormCounts(SCE, assay.type="TileMatrix")
  mainExpName(SCE) <- "FragmentCounts"
  rowData(SCE)$end <- rowData(SCE)$start + (tile - 1) 
  rowRanges(SCE) <- as(rowData(SCE), "GRanges")
  rownames(SCE) <- as.character(rowRanges(SCE))
  colData(SCE) <- proj@cellColData
  reducedDim(SCE, "UMAP") <- proj@embeddings$UMAP$df
  names(reducedDim(SCE, "UMAP")) <- c("UMAP1", "UMAP2")

  # NMF should just be a clone of this 
  message("Copying LSI scores to reducedDim(SCE, 'LSI')...")
  reducedDim(SCE, "LSI") <- proj@reducedDims$IterativeLSI$matSVD
  LSIdims <- ncol(proj@reducedDims$IterativeLSI$matSVD)
  names(reducedDim(SCE, "LSI")) <- paste0("LSI", seq_len(LSIdims))
  if (colDat) { 
    message("Copying LSI dimensions to colData(SCE) for iSEE visualization...")
    for (i in names(reducedDim(SCE, "LSI"))) colData(SCE)[, i] <- i
  }

  if (addNMF) {
    message("Fitting rank-", k, " NMF model to logCounts(SCE)...")
    metadata(SCE)$NMF <- RcppML::nmf(logcounts(SCE), k=k)
    message("Copying NMF hat matrix to reducedDim(SCE, 'NMF')...")
    reducedDim(SCE, "NMF") <- t(metadata(SCE)$NMF@h)
    NMFdims <- ncol(reducedDim(SCE, "NMF"))
    names(reducedDim(SCE, "NMF")) <- paste0("NMF", seq_len(NMFdims))
    if (colDat) { 
      message("Copying NMF dimensions to colData() for iSEE visualization...")
      for (i in names(reducedDim(SCE, "NMF"))) colData(SCE)[, i] <- i
    }
  }
     
  message("Done.")
  return(SCE)

}
