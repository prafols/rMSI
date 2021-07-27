#' groupComparisonByPixelsClosestToROIsMedianSpectra
#'
#' @param intensityMatrix Intensity matrix where rows represent pixels and columns represent peaks.
#' @param ROIs Vector with a code value for each pixel in the intensity matrix.
#' @param numberOfPixels Number of pixels closest to the median specrum of each ROI to use in the comparison.
#' @param testMethod Test method to use in the group comparison. Options are: Wilcoxon, Students or Kruskal.
#'
#' @return List containing the p.value and FC for all ions at each group comparision.
#' @export
#'
#' @examples
groupComparisonByPixelsClosestToROIsMedianSpectra <- function(intensityMatrix, ROIs, numberOfPixels, testMethod = "Kruskal")
{
  if(!any(testMethod == c("Wilcoxon","Students","Kruskal")))
  {
    writeLines("testMethod must be: Wilcoxon, Students or Kruskal \n")
    stop()
  }
  
  pixelsPerROI <- list()
  for(i in unique(ROIs))
  {
    median_spectrum <- apply(intensityMatrix[ROIs==i,],2, median) #Median spectrum of the ROI
    pixelsPerROI[[i]] <- order(apply(intensityMatrix[ROIs==i, ], 1, function(x) return(cor(x, median_spectrum))), decreasing = TRUE)[1:numberOfPixels]
  }
  
  
  results <- list()
  for(labelA in unique(ROIs))
  {
    rowsA <- which(ROIs == labelA)[pixelsPerROI[[labelA]]]
    for(labelB in unique(ROIs))
    {
      if(labelA != labelB)
      {
        if(!any(names(results) == paste0(labelB,"_vs_",labelA)))
        {
          rowsB <- which(ROIs == labelB)[pixelsPerROI[[labelB]]]
          
          lg2FC <- c()
          lg10pV <- c()
          
          for(ion in 1:ncol(intensityMatrix))
          {
            lg2FC <- c(lg2FC, log2(median(intensityMatrix[rowsA,ion])/median(intensityMatrix[rowsB,ion])))
            if(testMethod == "Wilcoxon") lg10pV <- c(lg10pV, -log10(wilcox.test(x = intensityMatrix[rowsA,ion],y = intensityMatrix[rowsB,ion])$p.value))
            if(testMethod == "Students") lg10pV <- c(lg10pV, -log10(t.test(x = intensityMatrix[rowsA,ion],y = intensityMatrix[rowsB,ion])$p.value))
            if(testMethod == "Kruskal") lg10pV <- c(lg10pV, -log10(kruskal.test(x = c(intensityMatrix[rowsA,ion], intensityMatrix[rowsB,ion]),g = rep(c("a","b"), each = numberOfPixels))$p.value))
          }
          results[[paste0(labelA,"_vs_",labelB)]] <- list(pV = lg10pV, FC = lg2FC)
        }
      }
    }
  }
  
  for(i in 1:length(results))
  {
    if(sum(is.infinite(results[[i]]$pV)) != 0)
    {
      results[[i]]$pV[is.infinite(results[[i]]$pV)] <- max(results[[i]]$pV[!is.infinite(results[[i]]$pV)])
    }
  }

  return(results)
}