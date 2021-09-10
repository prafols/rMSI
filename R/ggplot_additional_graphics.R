#########################################################################
#     rMSIproc - R package for MSI data processing
#     Copyright (C) 2019 Pere Rafols Soler
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
############################################################################


#' theme_black.
#' custom ggplot theme to display ion map images with a black background
#'
#' @param base_size 
#' @param base_family 
#'
theme_black <- function(base_size = 12, base_family = "") {
  
  ggplot2::theme_grey(base_size = base_size, base_family = base_family) + 
    
    ggplot2::theme(
      # Specify axis options
      axis.line = ggplot2::element_blank(),  
      axis.text.x = ggplot2::element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = ggplot2::element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = ggplot2::element_line(color = "white", size  =  0.2),  
      axis.title.x = ggplot2::element_text(size = base_size, color = "white", margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.title.y = ggplot2::element_text(size = base_size, color = "white", angle = 90, margin = ggplot2::margin(0, 10, 0, 0)),  
      axis.ticks.length = ggplot2::unit(0.3, "lines"),   
      # Specify legend options
      legend.background = ggplot2::element_rect(color = NA, fill = "black"),  
      legend.key = ggplot2::element_rect(color = "white",  fill = "black"),  
      legend.key.size = ggplot2::unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = ggplot2::element_text(size = base_size*0.8, color = "white"),  
      legend.title = ggplot2::element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "bottom",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      #legend.direction = "horizontal",  
      legend.box = NULL, 
      legend.spacing = ggplot2::unit(1,"line"),
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0, unit="line"),
      # Specify panel options
      panel.background = ggplot2::element_rect(fill = "black", color  =  NA),  
      panel.border = ggplot2::element_rect(fill = NA, color = "white"),  
      panel.grid.major = ggplot2::element_line(color = "grey35"),  
      panel.grid.minor = ggplot2::element_line(color = "grey20"),  
      panel.spacing = ggplot2::unit(0.5, "lines"),   
      # Specify facetting options
      strip.background = ggplot2::element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = ggplot2:: element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = ggplot2:: element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = ggplot2::element_rect(color = "black", fill = "black"),  
      plot.title = ggplot2::element_text(size = base_size*1.2, color = "white", hjust = 0.5,  margin = ggplot2::margin(0,0,10,0)),  
      plot.margin = ggplot2::unit(rep(1.2, 4), "lines")
      
    )
}

#' plotMassDriftImageG.
#' 
#' plots an image map displaying the mass shift at each raster position.
#'
#' @param peakMatrix an rMSIproc peak matrix.
#' @param peakList an rMSIproc peak list (no binning here!).
#' @param target_mass the target mass to represent.
#' @param error_range_ppm an error range in ppm to be represented around the target mass.
#' @param plot_rows number of rows to arrange multiple images in the plotting area.
#' @param plot_cols number of columns to arrange multiple images in the plotting area.
#' @param plot_byrow a boolean idicating if the plotted images must be sorted by rows.
#' @param plot_rotations a vector with the rotation in degree to apply to each image.
#' @param plot_mirror_X a vector of booleans idicatinc if each image must be flipped horizontally.
#' @param plot_mirror_Y a vector of booleans idicatinc if each image must be flipped vertically.
#' @param plot_margin a numeric value that determines the separation between images.
#' @param plot_labels text labels to be used for each image.
#' @param title_label Text label for the plot main title.
#' @param fixed_aspect_ratio set this flag to true to fix the aspect ratio of the ion images.
#'
#' @return a ggplot2 object.
#' @export
#' 
plotMassDriftImageG <- function(peakMatrix, peakList, target_mass, error_range_ppm,
                               plot_rows = 2, plot_cols = 2, plot_byrow = T, plot_rotations = rep(0, length(peakMatrix$names)),
                               plot_mirror_X =  rep(F, length(peakMatrix$names)), plot_mirror_Y =  rep(F, length(peakMatrix$names)),
                               plot_margin = 40, plot_labels = peakMatrix$names, title_label = "", fixed_aspect_ratio = F)
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
  {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  mass_drift <- c()
  intens <- c()
  for( i in 1:length(peakList))
  {
    iMZ <- which.min(abs(target_mass - peakList[[i]]$mass))
    if(length(iMZ) == 1)
    {
      mass_drift <- c(mass_drift, 1e6*(target_mass - peakList[[i]]$mass[iMZ])/target_mass) 
      intens <- c(intens, peakList[[i]]$intensity[iMZ]) #/sum(peakList[[i]]$intensity)) just set peakList TIC norm before entring here
    }
    else
    {
      mass_drift <- c(mass_drift, NA)
      intens <- c(intens, 0)
    }
    
  }
  
  mass_drift[mass_drift > error_range_ppm] <- NA
  mass_drift[mass_drift < -error_range_ppm] <- NA
  
  pltDf <- data.frame( x = peakMatrix$pos[,"x"], y = peakMatrix$pos[,"y"], mass_shift = mass_drift, intensity = intens) 
  
  #Calculate acq orders
  pltDf$acq_order <- rep(NA, nrow(pltDf))
  iStart <- 1
  for( i in 1:length(peakMatrix$numPixels))
  {
    iStop <- iStart + peakMatrix$numPixels[i] - 1
    if(i == 1)
    {
      pltDf$acq_order[iStart:iStop] <- rMSI::SortIDsByAcquisition( list(pos= pltDf[iStart:iStop, c("x", "y")] ) )
    }
    else
    {
      pltDf$acq_order[iStart:iStop] <- rMSI::SortIDsByAcquisition( list(pos= pltDf[iStart:iStop, c("x", "y")] ) ) + sum(peakMatrix$numPixels[1:(i-1) ]) 
    }
    iStart <- iStop + 1
  }
  
  #Apply Y offsets
  text_labels_yOffsets <- c(min(peakMatrix$pos[ 1:(peakMatrix$numPixels[1] - 1) ,"y"] ))
  offsetY <- 0
  if(length(peakMatrix$numPixels) > 1)
  {
    iStart <- 1
    for( i in 1:(length(peakMatrix$numPixels)-1)) 
    {
      iStop <- iStart + peakMatrix$numPixels[i] - 1
      offsetY <- 20 + max( peakMatrix$pos[ iStart:iStop ,"y"] ) + offsetY
      text_labels_yOffsets <- c(text_labels_yOffsets, offsetY)
      iStart <- iStop + 1
      pltDf$y[ iStart:(iStart + peakMatrix$numPixels[i+1] -1 ) ] <-  pltDf$y[ iStart:(iStart + peakMatrix$numPixels[i+1] -1 ) ] + offsetY
    }
  
  }
  
  pltDf$dataset <- unlist(lapply(1:length(peakMatrix$numPixels), function(i){ rep(paste0("img",i), peakMatrix$numPixels[i] )  }  ))
  pltDf <- pltDf[order(pltDf$acq_order), ]
  
  #Invert Y values
  maxAuxY <- max(pltDf$y)
  pltDf$y <- maxAuxY - pltDf$y
  text_labels_yOffsets <- maxAuxY - text_labels_yOffsets
  
  
  pltRas <- plotValuesImageG(peakMatrix = peakMatrix, pixel_values = pltDf$mass_shift, 
                            scale_label = sprintf("m/z shift at %.4f [ppm]", target_mass ), title_label = title_label,
                            plot_rows = plot_rows, plot_cols = plot_cols, plot_byrow = plot_byrow, plot_rotations = plot_rotations,
                            plot_mirror_X = plot_mirror_X, plot_mirror_Y = plot_mirror_Y, plot_margin = plot_margin, plot_labels = plot_labels,
                            gradient_scale_colours = rev(rainbow(n = 100, start = 0.1, end = 0.6)), 
                            gradient_scale_limits = c(-error_range_ppm, error_range_ppm), fixed_aspect_ratio = fixed_aspect_ratio)

  return( pltRas)
}


#' plotMassDriftG.
#' 
#' Plot the mass shift observed at a target mass.
#'
#' @param peakMatrix an rMSIproc peak matrix.
#' @param peakList an rMSIproc peak list (no binning here!).
#' @param target_mass the target mass to represent.
#' @param error_range_ppm an error range in ppm to be represented around the target mass.
#' @param min_SNR the minimum signal to noise ratio to be represented.
#' @param mass_offset a mass offset to shift the plot mass axis relative to the target mass.
#' @param title_label a string to be used as plot title. Chemical forumlas will be parsed to produce better results.
#' @param normalization a vector containing the normalization value for each pixel or NA if no normalization should be applied.
#' @param visible_legend a boolean specfing if a legend detailing the included MS images must be displayed.
#' @param legend_title the title for the legend.
#' @param N numbe of dots to display for each pixel, only the N highest SNR will be displayed.
#'
#' @return a ggplot2 object.
#' @export
#'
plotMassDriftG <- function(peakMatrix, peakList, target_mass, error_range_ppm, min_SNR = 10, mass_offset = 0, 
                           title_label="", normalization = NA, visible_legend = T, legend_title = "", N = 1)
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
  {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  mass_range <- (error_range_ppm*(target_mass+mass_offset))/1e6
  pkmass <- c()
  intens <- c()
  snr <- c()
  pixelID<- c()
  pixelX <- c()
  pixelY <- c()
  
  #apply the normalization
  if(length(normalization) == length(peakList))
  {
    for(i in 1:length(peakList))
    {
      if(normalization[i] > 0 )
      {
        peakList[[i]]$intensity <-  peakList[[i]]$intensity/normalization[i]
      }
      else
      {
        peakList[[i]]$intensity <- 0
      }
    }
  }
  
  cat("Getting peaks in the mass range...\n")
  pb <- txtProgressBar(min = 0, max = length(peakList), initial = 0, style = 3)
  for( i in 1:length(peakList))
  {
    setTxtProgressBar(pb, i)
    iMZ_min <- which.min(abs((target_mass+mass_offset-mass_range) - peakList[[i]]$mass))
    iMZ_max <- which.min(abs((target_mass+mass_offset+mass_range) - peakList[[i]]$mass))
    
    if( length(iMZ_min) > 0 && length(iMZ_max) > 0 )
    {
      if( iMZ_max >= iMZ_min )
      {
        keep_id <- c()
        for( iMZ in iMZ_min:iMZ_max  )
        {
          #if(peakList[[i]]$mass[iMZ] >= (target_mass+mass_offset-mass_range) && peakList[[i]]$mass[iMZ] <= (target_mass+mass_offset+mass_range))
          
          if(is.null(peakList[[i]]$SNR)) #No SNR info so just keep the peak
          {
            keep_id <- c(keep_id, iMZ)
          }
          else
          {
            if(peakList[[i]]$SNR[iMZ] >= min_SNR && peakList[[i]]$mass[iMZ] >= (target_mass+mass_offset-mass_range) && peakList[[i]]$mass[iMZ] <= (target_mass+mass_offset+mass_range))
            {
              keep_id <- c(keep_id, iMZ)
            }
          }
        }
        
        #Filter the dataframe to retain only N of the highest SNR
        if(length(keep_id) > 0)
        {
          if(is.null(peakList[[i]]$SNR)) #No SNR info so just keep the peak
          {
            dfSnrFilter <- data.frame( id = keep_id, SNR = rep(1, length(peakList[[i]]$mass[keep_id])))
          }
          else
          {
            dfSnrFilter <- data.frame( id = keep_id, SNR = peakList[[i]]$SNR[keep_id] )
          }
          if(nrow(dfSnrFilter) > N)
          {
            dfSnrFilter <- dfSnrFilter[order(dfSnrFilter$SNR, decreasing = T),] 
            dfSnrFilter <- dfSnrFilter[1:N, ]
          }
          
          pkmass <- c(pkmass, peakList[[i]]$mass[dfSnrFilter$id])
          intens <- c(intens, peakList[[i]]$intensity[dfSnrFilter$id])
          
          if(is.null(peakList[[i]]$SNR)) #No SNR info so just keep the peak
          {
            snr <- c(snr, 1) 
          }
          else
          {
            snr <- c(snr, peakList[[i]]$SNR[dfSnrFilter$id])
          }
          pixelID <- c(pixelID, rep(i, nrow(dfSnrFilter)))
          pixelX <- c(pixelX, rep( peakMatrix$pos[i, "x"], nrow(dfSnrFilter)))
          pixelY <- c(pixelY, rep( peakMatrix$pos[i, "y"], nrow(dfSnrFilter)))
        }
      }
    }
  }
  close(pb)
  
  #Saturate snr at min_SNR
  # sat_id <- which(snr > min_SNR)
  # if(length(sat_id) > 0)
  # {
  #   snr[sat_id] <- min_SNR
  # }

  pltDf <- data.frame( mass = pkmass, intensity = intens, SNR = snr, pixel = pixelID, x = pixelX, y = pixelY) 
  rm(pkmass)
  rm(intens)
  rm(snr)
  rm(pixelID)
  rm(pixelX)
  rm(pixelY)
  
  #Calculate acq orders
  pltDf$dataset <- rep(NA, nrow(pltDf))
  pxIDStart <- 1
  for( i in 1:length(peakMatrix$numPixels))
  {
    pxIDStop <- pxIDStart + peakMatrix$numPixels[i] - 1
    currRows <- which( pltDf$pixel %in% pxIDStart:pxIDStop)
    if(length(currRows) > 0)
    {
      pltDf$dataset[currRows] <- rep(paste0("img",i), length(currRows) )
      if(length(currRows) > 1)
      {
        acq_ordered_partial_rows <- rMSI::SortIDsByAcquisition( list(pos= pltDf[currRows, c("x", "y")] ) )
      }
      else
      {
        acq_ordered_partial_rows <- 1
      }
      pltDf[currRows, ] <- pltDf[ min(currRows) + acq_ordered_partial_rows - 1, ]
    }
    pxIDStart <- pxIDStop + 1
  }
  pltDf$acq_order <- 1:nrow(pltDf)
  
  pltPts <- ggplot2::ggplot(pltDf, ggplot2::aes(mass, -acq_order))
  pltPts <- pltPts + ggplot2::theme_bw(base_size = 15 )
  pltPts <- pltPts + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                           axis.ticks.y=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           legend.position=c(1, 1),
                           legend.justification = c(1, 1),
                           legend.direction = "vertical",
                           legend.background = ggplot2::element_rect(fill="transparent"),
                           legend.key = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_text(hjust=0.5, size = 12),
                           plot.title = ggplot2::element_text(hjust=0.5)
                           )
  pltPts <- pltPts + ggplot2::geom_point(ggplot2::aes(colour = dataset, alpha = SNR), size = 1, shape = 16)
  pltPts <- pltPts + ggplot2::labs( x = expression(~italic(m/z)~""), title = parse(text= chemFormula2Expression(title_label) ))
  pltPts <- pltPts +  ggplot2::scale_x_continuous(limits=c(target_mass+mass_offset-mass_range, target_mass+mass_offset+mass_range),
                                         sec.axis = ggplot2::sec_axis(~1e6*( (.) - target_mass)/(target_mass) , name = expression(~italic(m/z)~" shift [ppm]"),
                                                             #breaks = seq(  from = -error_range_ppm, to = error_range_ppm, by = round(error_range_ppm/5)) ))
                                                             breaks = seq( from = round(1e6*(mass_offset - mass_range)/(target_mass)), 
                                                                           to = round(1e6*(mass_offset + mass_range)/(target_mass)),
                                                                           by = round(error_range_ppm/5)) ))
  if(visible_legend)
  {
    pltPts <- pltPts + ggplot2::guides(alpha = "none", colour = ggplot2::guide_legend( ncol = 1, title = legend_title, override.aes = list(size=4)))   
  }
  else
  {
    pltPts <- pltPts + ggplot2::guides(alpha =  "none", colour =  "none") 
  }
  pltPts <- pltPts + ggplot2::geom_vline( linetype = 3, xintercept = target_mass, colour = "blue",  size=0.4)
  
  return(pltPts)
}



#' plotPeakImageG.
#' 
#' plots and ion map image using the rMSIproc peak matrix.
#' This function is equivalent to the rMSIproc::plotPeakImage but used ggplot to produce more publication-ready results.
#'
#' @param peakMatrix an rMSIproc peak matrix.
#' @param mass the ion mass to be plot.
#' @param normalization a vector containing the normalization value for each pixel or NA if no normalization should be applied.
#' @param plot_rows number of rows of the plotted matrix layout (only used if byrow == F).
#' @param plot_cols number of cols of the plotted matrix layout (only used if byrow == T).
#' @param plot_byrow a bool specifing if images must be arranged in rows.
#' @param plot_rotations a vector with the rotation of each images specified in degrees.
#' @param plot_mirror_X a bool vector specifing if each image must be flipped or not in X direction prior to rotation.
#' @param plot_mirror_Y a bool vector specifing if each image must be flipped or not in Y direction prior to rotation.
#' @param plot_margin the separation between plotted images.
#' @param plot_labels an alternative character vector with the labels to be displayed for each image.
#' @param title_label the main title of the plot
#' @param fixed_aspect_ratio set this flag to true to fix the aspect ratio of the ion images. 
#' @param display_colorbar set if the colour bar must be displayed.
#'
#' @return a ggplot2 object.
#' @export
#' 
plotPeakImageG <- function(peakMatrix, mass, normalization = NA,
                          plot_rows = 2, plot_cols = 2, plot_byrow = T, plot_rotations = rep(0, length(peakMatrix$names)),
                          plot_mirror_X =  rep(F, length(peakMatrix$names)), plot_mirror_Y =  rep(F, length(peakMatrix$names)),
                          plot_margin = 40, plot_labels = peakMatrix$names, title_label = "", fixed_aspect_ratio = F, display_colorbar = T)
{
  icol <- which.min(abs(mass-peakMatrix$mass))
  plotValues <- peakMatrix$intensity[,icol]

  #apply the normalization
  if(length(normalization) == length(plotValues))
  {
    plotValues <- plotValues/normalization
    plotValues[normalization == 0] <- 0
  }
  return (plotValuesImageG(peakMatrix =  peakMatrix, pixel_values = plotValues,
                          scale_label = sprintf("m/z %.4f", peakMatrix$mass[icol]), title_label = title_label, 
                          plot_rows = plot_rows, plot_cols = plot_cols, plot_byrow = plot_byrow, plot_rotations = plot_rotations, 
                          plot_mirror_X = plot_mirror_X, plot_mirror_Y =  plot_mirror_Y, plot_margin =  plot_margin,
                          plot_labels = plot_labels, fixed_aspect_ratio = fixed_aspect_ratio, display_colorbar = display_colorbar))
}

#' plotClusterImageG.
#' 
#' This function is just for convenience since calling plotValuesImageG with a factor produces the same results.
#' This function is equivalent to rMSIproc::plotClusterImage() but using the ggplot2.
#' 
#' @param peakMatrix an rMSIproc peak matrix.
#' @param clusters a vector with integer number according the cluster of each pixel.
#' @param title_label Text label for the plot main title.
#' @param plot_rows number of rows to arrange multiple images in the plotting area.
#' @param plot_cols number of columns to arrange multiple images in the plotting area.
#' @param plot_byrow a boolean idicating if the plotted images must be sorted by rows.
#' @param plot_rotations a vector with the rotation in degree to apply to each image.
#' @param plot_mirror_X a vector of booleans idicatinc if each image must be flipped horizontally.
#' @param plot_mirror_Y a vector of booleans idicatinc if each image must be flipped vertically. 
#' @param plot_margin a numeric value that determines the separation between images.
#' @param plot_labels text labels to be used for each image.
#' @param fixed_aspect_ratio set this flag to true to fix the aspect ratio of the ion images. 
#'
#' @return a ggplot2 object.
#' @export
#' 
plotClusterImageG <- function(peakMatrix, clusters, 
                              plot_rows = 2, plot_cols = 2, plot_byrow = T, plot_rotations = rep(0, length(peakMatrix$names)),
                              plot_mirror_X =  rep(F, length(peakMatrix$names)), plot_mirror_Y =  rep(F, length(peakMatrix$names)),
                              plot_margin = 40, plot_labels = peakMatrix$names, title_label = "", fixed_aspect_ratio = F)
{
  return (plotValuesImageG(peakMatrix =  peakMatrix, pixel_values = as.factor(clusters),
                           title_label = title_label, 
                           plot_rows = plot_rows, plot_cols = plot_cols, plot_byrow = plot_byrow, plot_rotations = plot_rotations, 
                           plot_mirror_X = plot_mirror_X, plot_mirror_Y =  plot_mirror_Y, plot_margin =  plot_margin,
                           plot_labels = plot_labels, fixed_aspect_ratio = fixed_aspect_ratio))
}

#' plotValuesImageG.
#'
#' Plot a raster image of arbitrary pixel values. This cna be use to plot PCA results for example.
#' This function is equivalent to the rMSIproc::plotValuesImage() but using ggplot2.
#' If the pixel_values is a "factor" object then discrete colors are used. This is useful to plot clustering results.
#'
#' @param peakMatrix an rMSIproc peak matrix.
#' @param pixel_values a vector with the pixel values, factor object can be used for a discrete image.
#' @param scale_label Text label for the plot scale bar.
#' @param title_label Text label for the plot main title.
#' @param plot_rows number of rows to arrange multiple images in the plotting area.
#' @param plot_cols number of columns to arrange multiple images in the plotting area.
#' @param plot_byrow a boolean idicating if the plotted images must be sorted by rows.
#' @param plot_rotations a vector with the rotation in degree to apply to each image.
#' @param plot_mirror_X a vector of booleans idicatinc if each image must be flipped horizontally.
#' @param plot_mirror_Y a vector of booleans idicatinc if each image must be flipped vertically. 
#' @param plot_margin a numeric value that determines the separation between images.
#' @param plot_labels text labels to be used for each image.
#' @param gradient_scale_colours alternative color scale for the image.
#' @param gradient_scale_limits alternative limits for the pixel values.
#' @param fixed_aspect_ratio set this flag to true to fix the aspect ratio of the ion images.
#' @param display_colorbar set if the colour bar must be displayed.
#'
#' @return a ggplot2 object.
#' @export
#'
plotValuesImageG <- function(peakMatrix, pixel_values, scale_label = "", title_label = "",
                          plot_rows = 2, plot_cols = 2, plot_byrow = T, plot_rotations = rep(0, length(peakMatrix$names)),
                          plot_mirror_X =  rep(F, length(peakMatrix$names)), plot_mirror_Y =  rep(F, length(peakMatrix$names)),
                          plot_margin = 40,  plot_labels = peakMatrix$names,
                          gradient_scale_colours = rev(rainbow(n = 100, start = 0, end = 0.6)), 
                          gradient_scale_limits = c(min(pixel_values), max(pixel_values)), 
                          fixed_aspect_ratio = F, display_colorbar = T)
{
  if (!requireNamespace("ggplot2", quietly = TRUE)) 
  {
    stop("Package ggplot2 needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  rasterData <- ArrangeMultipleImg2Plot( peakMat = peakMatrix, values = pixel_values,
                                                    nrow = plot_rows, ncol = plot_cols, byrow =  plot_byrow, margin = plot_margin,
                                                    rotations = plot_rotations, mirror_x = plot_mirror_X, mirror_y = plot_mirror_Y)  
  
  pltDf <- data.frame( x = rasterData$pos[,"x"], y = rasterData$pos[,"y"], intensity = rasterData$values)
  pltDf$y <- max(pltDf$y) - pltDf$y

  pltRas <- ggplot2::ggplot(pltDf, ggplot2::aes(x, y))
  pltRas <- pltRas + theme_black(base_size = 15 )
  pltRas <- pltRas +  ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                           panel.grid.minor = ggplot2::element_blank(),
                           axis.line=ggplot2::element_blank(),
                           axis.text.x=ggplot2::element_blank(),
                           axis.text.y=ggplot2::element_blank(),
                           axis.ticks=ggplot2::element_blank(),
                           axis.title.x=ggplot2::element_blank(),
                           axis.title.y=ggplot2::element_blank(),
                           panel.border=ggplot2::element_blank()
                           )
  if(class(pixel_values) == "factor")
  {
    pltRas <- pltRas +  ggplot2::geom_raster(ggplot2::aes(fill = as.factor(intensity))) 
    pltRas <- pltRas +  ggplot2::scale_fill_manual(values = rainbow(n=length(levels(pixel_values))), name = "Cluster:")
  }
  else
  {
    pltRas <- pltRas +  ggplot2::geom_raster(ggplot2::aes(fill = intensity, alpha = log(abs(intensity)+0.01)) ) #Sum 0.01 to avoid NaNs on log(intensity)
    pltRas <- pltRas +  ggplot2::scale_fill_gradientn(scale_label, na.value = "black", 
                                            colours = gradient_scale_colours, 
                                            limits = gradient_scale_limits,
                                            breaks = pretty(gradient_scale_limits,n=10))
    if(display_colorbar)
    {
    pltRas <- pltRas +  ggplot2::guides(alpha = "none", fill = ggplot2::guide_colourbar(title.position = "top",
                                                                                   title.hjust = 0.5,
                                                                                   barheight = 1,
                                                                                   barwidth = 25,
                                                                                   frame.colour = "white"))
    }
    else
    {
      pltRas <- pltRas +  ggplot2::guides(alpha = "none", fill = F)
    }
  }

  pltRas <- pltRas +  ggplot2::annotate("text", label = plot_labels,
                              size = 4, colour = "white", hjust = 0.5, vjust = 0,
                              x =  rasterData$lab_x,
                              y = rasterData$lab_y - 0.9*plot_margin
                              )
  pltRas <- pltRas +  ggplot2::labs(title = parse(text= chemFormula2Expression(title_label) ))
  if(fixed_aspect_ratio)
  {
    pltRas <- pltRas +  ggplot2::coord_fixed() #This is used to get fixed asspect ratios but it is causing some issues 
  }
  
  return(pltRas)
}


#' chemFormula2Expression.
#' 
#' Converts a string containing a chemical forumla to an R expression that allows plotting the formula with sub-indices properly.
#'
#' @param strFormula a string containing a chemical formula.
#'
#' @return a formated R expression that must be used with parse().
#'
chemFormula2Expression <- function(strFormula)
{
  numPos <- unlist(gregexpr('\\(?[0-9,.]', strFormula))
  rawFormulaStr <- ""
  antWasNumber <- F
  for( i in 1:nchar(strFormula))
  {
    if( i %in% numPos)
    {
      #Parse number
      rawFormulaStr <- paste0(rawFormulaStr, "[", substr(strFormula, start = i, stop = i), "]")
      antWasNumber <- T
    }
    else
    {
      #Parse letter
      if( antWasNumber)
      {
        rawFormulaStr <- paste0(rawFormulaStr, "*")
      }
      rawFormulaStr <- paste0(rawFormulaStr, substr(strFormula, start = i, stop = i))
      antWasNumber <- F
    }
  }
  
  #return(parse(text = rawFormulaStr))
  return(rawFormulaStr)
}



#' plotMassDriftCompared.
#' 
#' Creates a plot in a grid to compare the spectral alignment of an ion before and after the rMSIproc processing.
#' For a given mass this function will display in a single figure:
#'  - The ion map.
#'  - The alignment before the processing (RAW).
#'  - The alignment after the processing (Aligned).
#'
#' @param peakMatrix_RAW an rMSIproc peak matrix obtained without the alginment routine (RAW).
#' @param peakMatrix_ALNG an rMSIproc peak matrix obtained with the alginment routine (Aligned).
#' @param peaklist_RAW a peak list generated with rMSIproc without the alginment routine (RAW).
#' @param peaklist_ALNG a peak list generated with rMSIproc with the alginment routine (Aligned).
#' @param mass the ion mass to be plot.
#' @param error_range_ppm a mass range to be plot specified in ppm.
#' @param min_SNR peaks with a signal to noise ratio below this parameter will be discarded.
#' @param mass_offset_RAW a mass offset applied to the RAW data to manually compensate possible mass shifts.
#' @param title_label the main title of the plot.
#' @param normalization_RAW  a vector containing the normalization value for each pixel in the RAW data or NA if no normalization should be applied.
#' @param normalization_ALNG a vector containing the normalization value for each pixel in the aligned data or NA if no normalization should be applied.
#' @param ion_map_plot_cols number of columns to arrange multiple images in the plotting area.
#' @param ion_map_plot_rows number of rows to arrange multiple images in the plotting area.
#' @param ion_map_plot_byrow a boolean idicating if the plotted images must be sorted by rows.
#' @param ion_map_plot_rotations a vector with the rotation in degree to apply to each image.
#' @param ion_map_plot_labels text labels to be used for each image.
#' @param ion_map_plot_margin a numeric value that determines the separation between images.
#' @param ion_map_plot_mirror_X a vector of booleans idicatinc if each image must be flipped horizontally.
#' @param ion_map_plot_mirror_Y a vector of booleans idicatinc if each image must be flipped vertically. 
#' @param ion_map_fixed_aspect_ratio set this flag to true to fix the aspect ratio of the ion images.
#' 
#' @return a ggplot2 object.
#' @export
#' 
plotMassDriftComparedG <- function(peakMatrix_RAW, peakMatrix_ALNG, peaklist_RAW, peaklist_ALNG, mass,
                                  error_range_ppm=400, min_SNR = 10, mass_offset_RAW = 0, 
                                  title_label="", normalization_RAW = NA, normalization_ALNG = NA,
                                  ion_map_plot_cols = 2, ion_map_plot_rows = 2, ion_map_plot_byrow = T,
                                  ion_map_plot_rotations = rep(0, length(peakMatrix_ALNG$names)), 
                                  ion_map_plot_labels = peakMatrix_ALNG$names,
                                  ion_map_plot_margin = 40, 
                                  ion_map_plot_mirror_X =  rep(F, length(peakMatrix_ALNG$names)), 
                                  ion_map_plot_mirror_Y =  rep(F, length(peakMatrix_ALNG$names)), 
                                  ion_map_fixed_aspect_ratio = F)
{
  if (!requireNamespace("gridExtra", quietly = TRUE)) 
  {
    stop("Package gridExtra needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  pltDriftRaw <- plotMassDriftG(peakMatrix_RAW, peaklist_RAW, mass, error_range_ppm, min_SNR = min_SNR, 
                               mass_offset = mass_offset_RAW, title_label = "", normalization = normalization_RAW, visible_legend = T, legend_title = "RAW")
  
  pltDriftAlng <- plotMassDriftG(peakMatrix_ALNG, peaklist_ALNG, mass, error_range_ppm, min_SNR = min_SNR, 
                                mass_offset = 0, title_label = "", normalization = normalization_ALNG, visible_legend = T, legend_title = "Aligned")
  
  pltIon<-plotPeakImageG(peakMatrix = peakMatrix_ALNG, mass = mass, normalization = normalization_ALNG, title_label = title_label,
                        plot_cols = ion_map_plot_cols, plot_rows = ion_map_plot_rows, plot_byrow = ion_map_plot_byrow, 
                        plot_rotations = ion_map_plot_rotations, plot_labels = ion_map_plot_labels,
                        plot_margin = ion_map_plot_margin, plot_mirror_X =  ion_map_plot_mirror_X, plot_mirror_Y =  ion_map_plot_mirror_Y, fixed_aspect_ratio = ion_map_fixed_aspect_ratio )
  
  
  pGrid<-gridExtra::arrangeGrob(pltIon, pltDriftRaw, pltDriftAlng,
                                 layout_matrix = rbind(c(1, 1),
                                                       c(2, 3))
                                )
  
  return(pGrid)
}

#' plotMassDriftImageComparedG.
#' 
#' plots the image map displaying the mass shift at each raster position of two datasets.
#' this allows to compare the results of the alignment routine (raw vs. aliged).
#'
#' @param peakMatrix_RAW an rMSIproc peak matrix obtained without the alginment routine (RAW).
#' @param peakMatrix_ALNG an rMSIproc peak matrix obtained with the alginment routine (Aligned).
#' @param peaklist_RAW a peak list generated with rMSIproc without the alginment routine (RAW).
#' @param peaklist_ALNG a peak list generated with rMSIproc with the alginment routine (Aligned).
#' @param mass the ion mass to be plot.
#' @param error_range_ppm  a mass range to be plot specified in ppm.
#' @param mass_offset_RAW a mass offset applied to the RAW data to manually compensate possible mass shifts.
#' @param title_label the main title of the plot.
#' @param ion_map_plot_cols number of columns to arrange multiple images in the plotting area.
#' @param ion_map_plot_rows number of rows to arrange multiple images in the plotting area.
#' @param ion_map_plot_byrow a boolean idicating if the plotted images must be sorted by rows.
#' @param ion_map_plot_rotations a vector with the rotation in degree to apply to each image.
#' @param ion_map_plot_labels text labels to be used for each image.
#' @param ion_map_plot_margin a numeric value that determines the separation between images.
#' @param ion_map_plot_mirror_X a vector of booleans idicatinc if each image must be flipped horizontally.
#' @param ion_map_plot_mirror_Y a vector of booleans idicatinc if each image must be flipped vertically. 
#' @param ion_map_fixed_aspect_ratio set this flag to true to fix the aspect ratio of the ion images.
#'
#' @return a ggplot2 object.
#' @export
#'
plotMassDriftImageComparedG <- function(peakMatrix_RAW, peakMatrix_ALNG, peaklist_RAW, peaklist_ALNG, mass,
                                        error_range_ppm=400, mass_offset_RAW = 0, 
                                        title_label="",
                                        ion_map_plot_cols = 2, ion_map_plot_rows = 2, ion_map_plot_byrow = T,
                                        ion_map_plot_rotations = rep(0, length(peakMatrix_ALNG$names)), 
                                        ion_map_plot_labels = peakMatrix_ALNG$names,
                                        ion_map_plot_margin = 40, 
                                        ion_map_plot_mirror_X =  rep(F, length(peakMatrix_ALNG$names)), 
                                        ion_map_plot_mirror_Y =  rep(F, length(peakMatrix_ALNG$names)), 
                                        ion_map_fixed_aspect_ratio = F)
{
  if (!requireNamespace("gridExtra", quietly = TRUE)) 
  {
    stop("Package gridExtra needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  pltRaw <- plotMassDriftImageG(peakMatrix = peakMatrix_RAW, peakList = peaklist_RAW,
                                              target_mass = mass + mass_offset_RAW, error_range_ppm = error_range_ppm, 
                                              plot_rows = ion_map_plot_rows, 
                                              plot_cols = ion_map_plot_cols,
                                              plot_byrow = ion_map_plot_byrow, 
                                              plot_rotations = ion_map_plot_rotations, 
                                              plot_mirror_X = ion_map_plot_mirror_X, 
                                              plot_mirror_Y = ion_map_plot_mirror_Y,
                                              plot_margin = ion_map_plot_margin, 
                                              plot_labels = ion_map_plot_labels, 
                                              title_label = paste("RAW:",title_label), 
                                              fixed_aspect_ratio = ion_map_fixed_aspect_ratio )
  
  pltAlng <- plotMassDriftImageG(peakMatrix = peakMatrix_ALNG, peakList = peaklist_ALNG,
                                              target_mass = mass, error_range_ppm = error_range_ppm, 
                                              plot_rows = ion_map_plot_rows, 
                                              plot_cols = ion_map_plot_cols,
                                              plot_byrow = ion_map_plot_byrow, 
                                              plot_rotations = ion_map_plot_rotations, 
                                              plot_mirror_X = ion_map_plot_mirror_X, 
                                              plot_mirror_Y = ion_map_plot_mirror_Y,
                                              plot_margin = ion_map_plot_margin, 
                                              plot_labels = ion_map_plot_labels, 
                                              title_label = paste("Aligned:",title_label), 
                                              fixed_aspect_ratio = ion_map_fixed_aspect_ratio )
  
  pGrid<-gridExtra::arrangeGrob(pltRaw, pltAlng, nrow = 1 )
  return(pGrid)

}
