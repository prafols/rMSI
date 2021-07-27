#' plotTopIonsInROIsComparaison
#'
#' @param topIonsList List returned by the function "getTopIonsFromGroupComparison"
#' @param title Title of the plots
#'
#' @return
#' @export
#'
#' @examples
plotTopIonsInROIsComparaison <- function(topIonsList, title = "")
{
  library(ggplot2)
  topIonsList$Up$siginfScore <- sqrt((((topIonsList$Up$`log2(FC)`)/max((topIonsList$Up$`log2(FC)`)))^2) + (((topIonsList$Up$`-log10(p.value)`)/max((topIonsList$Up$`-log10(p.value)`)))^2))
  topIonsList$Up$tag <- rep("Up-regulated", times = length(topIonsList$Up$siginfScore))
  topIonsList$Up <- topIonsList$Up[order(topIonsList$Up$siginfScore,decreasing = F),]
  
  topIonsList$Down$siginfScore <- -sqrt((((topIonsList$Down$`log2(FC)`)/max(topIonsList$Down$`log2(FC)`))^2) + (((topIonsList$Down$`-log10(p.value)`)/max(topIonsList$Down$`-log10(p.value)`))^2))
  topIonsList$Down$tag <- rep("Down-regulated", times = length(topIonsList$Down$siginfScore))
  topIonsList$Down <- topIonsList$Down[order(topIonsList$Down$siginfScore,decreasing = F),]
  
  topIonsCombined <- rbind(topIonsList$Down, topIonsList$Up)
  topIonsCombined$`m/z` <- factor(format(topIonsCombined$`m/z`,digits = 7), format(topIonsCombined$`m/z`,digits = 7))
  topIonsCombined$tag <- factor(topIonsCombined$tag, unique(topIonsCombined$tag)[2:1])
  print(ggplot2::ggplot() + ggplot2::geom_col(ggplot2::aes(y = topIonsCombined$`m/z`, x = topIonsCombined$`log2(FC)`, fill = topIonsCombined$`-log10(p.value)` )) + 
          ggplot2::theme_bw() + ggplot2::labs( y = "m/z", x = "log2(FC)", fill = "-log10(p.Value)", title = title, subtitle = "Significative ions") + 
          ggplot2::scale_fill_viridis_c())
  return(topIonsCombined)
}