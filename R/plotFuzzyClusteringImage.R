
#' plotFuzzyClusteringImageG
#'
#' @param peakMatrix rMSIprocPeakMatrix object.
#' @param clusters  Integer vector that indicates at which cluster belongs each pixel.
#' @param membership Numeric vector that contains the membership value of the pixels to the cluster they belong in clusters. 
#' @param title Title for the plot if desider. By default no title.
#' @param contrastAlpha If TRUE maps the range of memberships to the full range of alpha (0 to 1), increasing the contrast between pixels with different memberships.
#' @param plotMembershipMap If TRUE plots the membership image and the hard clustering at different images
#'
#' @return
#' @export
#'
#' @examples
plotFuzzyClusteringImageG <- function(peakMatrix, clusters, membership, title = "", contrastAlpha = FALSE, plotMembershipMap = FALSE)
{
  rasterData <- rMSIproc:::ArrangeMultipleImg2Plot(peakMat = peakMatrix,
                                                   values = clusters,
                                                   nrow = 1,
                                                   ncol = 2,byrow = T)
  
  pltDf <- data.frame( x = rasterData$pos[,"x"], y = rasterData$pos[,"y"], intensity = rasterData$values, alpha = membership)
  pltDf$y <- max(pltDf$y) - pltDf$y
  
  if(!plotMembershipMap)
  {
    if(contrastAlpha)
    {
      pltRas01 <- ggplot2::ggplot(data = pltDf, mapping = ggplot2::aes(x,y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = as.factor(intensity), alpha = alpha)) + ggplot2::coord_equal() +
        ggplot2::scale_fill_brewer(palette = "Accent", name = "cluster") +
        ggplot2::scale_alpha_continuous(name = "membership") +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.line=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks=ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       panel.border=ggplot2::element_blank(),
                       panel.background=ggplot2::element_rect(fill="white",colour = "white")) +
        ggplot2::labs(title = paste(title))
      print(pltRas01)
    }
    else
    {
      pltRas01 <- ggplot2::ggplot(data = pltDf, mapping = ggplot2::aes(x,y)) +
        ggplot2::geom_tile(ggplot2::aes(fill = as.factor(intensity), alpha = alpha)) + ggplot2::coord_equal() +
        ggplot2::scale_fill_brewer(palette = "Accent", name = "cluster") +
        ggplot2::scale_alpha_continuous(range = c(0,1), limits = c(0,1), name = "membership") +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       axis.line=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.text.y=ggplot2::element_blank(),
                       axis.ticks=ggplot2::element_blank(),
                       axis.title.x=ggplot2::element_blank(),
                       axis.title.y=ggplot2::element_blank(),
                       panel.border=ggplot2::element_blank(),
                       panel.background=ggplot2::element_rect(fill="white",colour = "white")) +
        ggplot2::labs(title = paste(title))
      print(pltRas01)
    }
  }
  else
  {
    pltRas01 <- ggplot2::ggplot(data = pltDf, mapping = ggplot2::aes(x,y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = as.factor(intensity))) + ggplot2::coord_equal() +
      ggplot2::scale_fill_brewer(palette = "Accent", name = "cluster") +
      ggplot2::scale_alpha_continuous(name = "membership") +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     panel.border=ggplot2::element_blank(),
                     panel.background=ggplot2::element_rect(fill="white",colour = "white")) +
      ggplot2::labs(title = paste(title))
    print(pltRas01)
    
    pltRas02 <- ggplot2::ggplot(data = pltDf, mapping = ggplot2::aes(x,y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = alpha)) + ggplot2::coord_equal() +
      ggplot2::scale_fill_viridis_c(name = "membership") +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                     panel.grid.minor = ggplot2::element_blank(),
                     axis.line=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     panel.border=ggplot2::element_blank(),
                     panel.background=ggplot2::element_rect(fill="white",colour = "white")) +
      ggplot2::labs(title = paste(title))
    print(pltRas02)
  }
}
