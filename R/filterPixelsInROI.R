#' filterPixelsInROI
#'
#' @param intensityMatrix Intensity matrix where rows represent pixels and columns represent peaks
#' @param ROIs Vector with a code value for each pixel in the intensity matrix.
#' @param topPercentage Percentage of top intensity pixels to remove.
#' @param lowPercentage Percentage of lower intensity pixels to remove.
#' @param magnitude Magnitude to use in filtering. RMS and TIC are the options allowed.
#' @param plotResults Flag to plot the results of the filtering in a violin plot.
#'
#' @return
#' @export
#'
#' @examples
filterPixelsInROI <- function(intensityMatrix, ROIs, topPercentage = 0.01, lowPercentage = 0.1, magnitude = "RMS", plotResults = FALSE)
{
  pixelsToUse <- c()
  if(magnitude == "TIC")
  {
    TIC <- apply(intensityMatrix, 1, sum)
    for(label in as.vector(unique(ROIs)))
    {
      ROIpixels <- which(ROIs == label)
      topPixels <- round(length(ROIpixels)*topPercentage)
      bottomPixels <- round(length(ROIpixels)*lowPercentage)
      ROIpixelsByOrder <- order(TIC[ROIpixels],decreasing = T)
      pixelsToUse <- c(pixelsToUse, ROIpixels[ROIpixelsByOrder[-c(1:topPixels, 1+length(ROIpixels)-(bottomPixels:1))]])
    }
    
    fROIs <- ROIs[pixelsToUse]
    fTIC <- apply(intensityMatrix[pixelsToUse,], 1, sum)
    if(plotResults)
    {
      g1 <- ggplot2::ggplot() + ggplot2::geom_violin(ggplot2::aes(x = c(ROIs,fROIs), y = c(TIC,fTIC), fill = rep(c("All pixels","Filtered pixels"), times = c(length(TIC), length(fTIC))))) +
        ggplot2::theme_bw() + ggplot2::labs(x="", y = "TIC", fill = "", title = paste0("ROI pixels filtering by ",magnitude)) + ggplot2::theme(axis.text.x = ggplot2::element_blank())
      print(g1)
    }
  }
  else if(magnitude == "RMS")
  {
    RMS <- apply(intensityMatrix, 1, function(x) sqrt(sum(x*x)/length(x)))
    for(label in as.vector(unique(ROIs)))
    {
      ROIpixels <- which(ROIs == label)
      topPixels <- round(length(ROIpixels)*topPercentage)
      bottomPixels <- round(length(ROIpixels)*lowPercentage)
      ROIpixelsByOrder <- order(RMS[ROIpixels],decreasing = T)
      pixelsToUse <- c(pixelsToUse, ROIpixels[ROIpixelsByOrder[-c(1:topPixels, 1+length(ROIpixels)-(bottomPixels:1))]])
    }
    
    fROIs <- ROIs[pixelsToUse]
    fRMS <- apply(intensityMatrix[pixelsToUse,], 1, function(x) sqrt(sum(x*x)/length(x)))
    if(plotResults)
    {
      g1 <- ggplot2::ggplot() + ggplot2::geom_violin(ggplot2::aes(x = c(ROIs,fROIs), y = c(RMS,fRMS), fill = rep(c("All pixels","Filtered pixels"), times = c(length(RMS), length(fRMS))))) +
        ggplot2::theme_bw() + ggplot2::labs(x="", y = "RMS", fill = "", title = paste0("ROI pixels filtering by ",magnitude)) + ggplot2::theme(axis.text.x = ggplot2::element_blank())
      print(g1)
    }
  } 
  else 
  {
    writeLines("Magnitude must be: TIC or RMS")
    stop()
  }
  return(pixelsToUse)
}