#' turn an ArchR project into a SingleCellExperiment 
#'
#' The default tile size is 500bp, just like ArchR's default. 
#'
#' @param proj          an ArchRproject 
#' @param tile          tile size (500)
#' @param addNMF        add an NMF decomposition of logCounts? (FALSE)
#' @param k             rank for the NMF decomposition on logCounts (30) 
#' @param ...           additional arguments to pass to ArchR::addTileMatrix
#'
#' @return              a SingleCellExperiment 
#'
#' @details presumes LSI and UMAP. you will almost certainly want to tweak this.
#'
#' @import ArchR
#' @import scater
#' @import GenomicRanges
#' @import SingleCellExperiment
#'
#' @export
#'
archRtoSCE <- function(proj, tile=500, addNMF=FALSE, k=30) { 

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
  reducedDim(SCE, "LSI") <- proj@reducedDims$IterativeLSI$matSVD
  LSIdims <- proj@reducedDims$IterativeLSI$matSVD
  names(reducedDim(SCE, "LSI")) <- paste0("LSI", seq_len(LSIdims))
  for (i in names(reducedDim(SCE, "LSI"))) colData(SCE)[, i] <- i

  if (addNMF) warning("NMF support hasn't been bolted on yet, yell at Tim")

  return(SCE)

}
