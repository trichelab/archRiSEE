#' convenience wrapper for iSEE on an ArchR-derived singlecellexperiment
#'
#' Note that this is a quick and dirty affair so don't expect much. 
#'
#' @param x             a SingleCellExperiment with UMAP and LSI reducedDims
#' @param colorColumn   the colData column to color plots ("mitoCluster")
#'
#' @details presumes NMF and UMAP. you will almost certainly want to tweak this.
#'
#' @import iSEE
#'
#' @export
#'
iSEEarchR <- function(x, colorColumn = "Clusters") { 

  # reason for this will become obvious 
  stopifnot(c("UMAP","LSI") %in% reducedDimNames(x))

  # if LSI columns aren't in colData already, add them 
  # could do the same for NMF if it is in reducedDimNames(x)... 
  LSIdim <- ncol(reducedDim(ewsSub, "LSI"))
  if (sum(grepl("LSI", names(colData(ewsSub)), ignore.case=TRUE)) < LSIdim) { 
    for (i in seq_len(LSIdim)) { 
      colData(x)[, paste0("LSI", i)] <- reducedDim(x, "LSI")[, i]
    }
  }

  # launch
  iSEE(x,
     initial = 
       list(UMAP = new("ReducedDimensionPlot",
                       FontSize = 1.5,
                       PointSize = 4,
                       VisualBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = colorColumn,
                       Type = "NMF_UMAP"
                       ),
            FACTORS = new("ColumnDataPlot",
                       YAxis = "TSSEnrichment",
                       XAxis = "Column data",
                       XAxisColumnData = "Sample",
                       FontSize = 1.5,
                       PointSize = 4,
                       DataBoxOpen = TRUE,
                       SelectionBoxOpen = TRUE,
                       ColorBy = "Column data",
                       ColorByColumnData = colorColumn,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       ),
            COLDAT = new("ColumnDataTable",
                       SelectionBoxOpen = TRUE,
                       ColumnSelectionDynamicSource = TRUE,
                       ColumnSelectionSource = "ReducedDimensionPlot1"
                       )
            )
     )
}
