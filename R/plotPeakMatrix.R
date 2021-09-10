#########################################################################
#     rMSIproc - R package for MSI data processing
#     Copyright (C) 2017 Pere Rafols Soler
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

#' plotPeakImage.
#' 
#' plot the ion image map of a given mass or column of a rMSIproc peak matrix object.
#' If the peak matrix contains data from various datasets the images will be layout horizontally or vertically.
#' At leas mz or column must be specified.
#'
#' @param peakMatrix the peak matrix in an rMSIproc object.
#' @param mz the peak mass to plot, the nearest peak mass will be ploted.
#' @param column the column of the peak matrix to plot.
#' @param matrix the name of the peak matrix to plot.
#' @param normalization the name of normalization to use.
#' @param nrow number of rows of the plotted matrix layout (only used if byrow == F)
#' @param ncol number of cols of the plotted matrix layout (only used if byrow == T)
#' @param byrow a bool specifing if images must be arranged in rows.
#' @param margin the separation between plotted images.
#' @param img_names a character vector with img names in the desierd plotting order.
#' @param labels an alternative character vector with the labels to be displayed for each image.
#' @param rotations a vector with the rotation of each images specified in degrees.
#' @param mirror_x a bool vector specifing if each image must be flipped or not in X direction prior to rotation.
#' @param mirror_y a bool vector specifing if each image must be flipped or not in Y direction prior to rotation.
#' @param pixel_size_um the pixel resolution in um.
#' @param light the lighting of the plotted image.
#' @export
#'
plotPeakImage <- function(peakMatrix, mz = NULL, column = NULL, matrix = "intensity", normalization = NULL, 
                          nrow = 2, ncol = 2, byrow = T, 
                          margin = 20, img_names = peakMatrix$names, 
                          labels = img_names,
                          rotations = rep(0, length(img_names)),
                          mirror_x = rep(F, length(img_names)), mirror_y = rep(F, length(img_names)), 
                          pixel_size_um = 100, light = 8)
{
  if(is.null(peakMatrix[[matrix]]))
  {
    stop("The selected peak matrix does not exist.\n")
  }
  
  if(is.null(mz) && is.null(column))
  {
    stop("At least mass or a column must be specified")
  }
  
  if( !is.null(mz))
  {
    column <- which.min(abs(peakMatrix$mass - mz))
  }
  
  if( !is.null(normalization))
  {
    if( is.null(peakMatrix$normalizations))
    {
      stop("No normalizations avaliable for the supplied peak matrix")
    }
    
    normVals <- peakMatrix$normalizations[[normalization]]
    if(is.null(normVals))
    {
      stop(paste0("No normalization found with the name:", normalization))
    }
  }
  else
  {
    normVals <- rep(1, nrow(peakMatrix$pos))
  }
  
  rasterData <- ArrangeMultipleImg2Plot(peakMatrix, peakMatrix[[matrix]][,column]/normVals, nrow, ncol, byrow, margin, img_names, rotations, mirror_x, mirror_y)
  rMSI::PlotValues(rasterData$pos, rasterData$values, rotate = 0, scale_title = sprintf("m/z %.4f", peakMatrix$mass[column]), pixel_size_um = pixel_size_um, vlight = light,
                   labels_x = rasterData$lab_x, labels_y = rasterData$lab_y, labels_text = labels)
}

#' plotClusterImage.
#'
#' Plot a segmentation image with the user-given clusters.
#'
#' @param peakMatrix the peak matrix in an rMSIproc object.
#' @param clusters a vector with integer number according the cluster of each pixel.
#' @param nrow number of rows of the plotted matrix layout (only used if byrow == F)
#' @param ncol number of cols of the plotted matrix layout (only used if byrow == T)
#' @param byrow a bool specifing if images must be arranged in rows.
#' @param margin the separation between plotted images.
#' @param img_names a character vector with img names in the desierd plotting order.
#' @param labels an alternative character vector with the labels to be displayed for each image.
#' @param rotations a vector with the rotation of each images specified in degrees.
#' @param mirror_x a bool vector specifing if each image must be flipped or not in X direction prior to rotation.
#' @param mirror_y a bool vector specifing if each image must be flipped or not in Y direction prior to rotation.
#' @param pixel_size_um the pixel resolution in um.
#'
#' @return a vector with the used color for each cluster sorted according clustering numering in assending order.
#' @export
#'
plotClusterImage <- function(peakMatrix, clusters, 
                             nrow = 2, ncol = 2, byrow = T, 
                             margin = 20, img_names = peakMatrix$names, 
                             labels = img_names,
                             rotations = rep(0, length(img_names)),
                             mirror_x = rep(F, length(img_names)), mirror_y = rep(F, length(img_names)),
                             pixel_size_um = 100)
{
  rasterData <- ArrangeMultipleImg2Plot(peakMatrix, clusters, nrow, ncol, byrow, margin, img_names, rotations, mirror_x, mirror_y)
  clusColors <- rMSI::PlotClusterImage(rasterData$pos, rasterData$values, rotate = 0, pixel_size_um, 
                                       labels_x = rasterData$lab_x, labels_y = rasterData$lab_y, labels_text = labels)
  legend( "topright", col = clusColors, legend =  paste0("Clus_",sort(unique(clusters), decreasing = F)), pch = 20, horiz = F,  bty = "n", text.col = "white" )
  return(clusColors)
}

#' plotValuesImage
#'
#' Plot arbitrary values in a MS images raster positions.
#'
#' @param peakMatrix the peak matrix in an rMSIproc object.
#' @param values a vector with the values to be plotted in each pixel.
#' @param nrow number of rows of the plotted matrix layout (only used if byrow == F)
#' @param ncol number of cols of the plotted matrix layout (only used if byrow == T)
#' @param byrow a bool specifing if images must be arranged in rows.
#' @param margin the separation between plotted images.
#' @param img_names a character vector with img names in the desierd plotting order.
#' @param labels an alternative character vector with the labels to be displayed for each image.
#' @param rotations a vector with the rotation of each images specified in degrees.
#' @param mirror_x a bool vector specifing if each image must be flipped or not in X direction prior to rotation.
#' @param mirror_y a bool vector specifing if each image must be flipped or not in Y direction prior to rotation.
#' @param pixel_size_um the pixel resolution in um.
#' @param scale_title the label for the color scale.
#' @param light the lighting of the plotted image.
#' @export
#'
plotValuesImage <- function(peakMatrix, values, 
                             nrow = 2, ncol = 2, byrow = T, 
                             margin = 20, img_names = peakMatrix$names, 
                             labels = img_names,
                             rotations = rep(0, length(img_names)),
                             mirror_x = rep(F, length(img_names)), mirror_y = rep(F, length(img_names)),
                             pixel_size_um = 100, scale_title = "", light = 8)
{
  rasterData <- ArrangeMultipleImg2Plot(peakMatrix, values, nrow, ncol, byrow, margin, img_names, rotations, mirror_x, mirror_y)
  rMSI::PlotValues(rasterData$pos, rasterData$values, rotate = 0, scale_title = scale_title, pixel_size_um = pixel_size_um, vlight = light,
                   labels_x = rasterData$lab_x, labels_y = rasterData$lab_y, labels_text = labels)
}



#' ArrangeMultipleImg2Plot
#' 
#' Prepare multiple images to be plotted with arbitrary data that must be specified as a vector using the values parameter.
#' The values vector must be sorted according the pos matrix of the provided peak matrix object (peakMat).
#' The images will be plotted laidout in a matrix according the parameters nrow, ncol and byrow which follows the
#' same conventions as the R matrix() function.
#' A different order of the images can be specified with the img_name parameter. Just provide the images names to plotted in the desierd order.
#' Additionally each single image can be flip horizontally and vertically and/or rotated to any angle using the parameters mirror_x, mirror_y and rotations.
#' 
#' @param peakMat an rMSIproc peak matrix.
#' @param values a vector of pixel values following the same order as pixels appear in the peak matrix.
#' @param nrow number of rows of the plotted matrix layout (only used if byrow == F)
#' @param ncol number of cols of the plotted matrix layout (only used if byrow == T)
#' @param byrow a bool specifing if images must be arranged in rows.
#' @param margin the separation between plotted images.
#' @param img_names a character vector with img names in the desierd plotting order.
#' @param rotations a vector with the rotation of each images specified in degrees.
#' @param mirror_x a bool vector specifing if each image must be flipped or not in X direction prior to rotation.
#' @param mirror_y a bool vector specifing if each image must be flipped or not in Y direction prior to rotation.
#'
#' @return a list object conaining the arranged position matrix, the pixel values and the labels positions.
#'
ArrangeMultipleImg2Plot <- function( peakMat, values, nrow, ncol, byrow = T, margin = 20, img_names = peakMat$names, 
                                rotations = rep(0, length(img_names)),
                                mirror_x = rep(F, length(img_names)), mirror_y = rep(F, length(img_names)) )
{
  #Re-Order data according the image names
  lstRefac <- list()
  iRow <- 1
  iCol <- 1
  nextXOffset <- 0
  nextYOffset <- 0
  partialOffset <- 0
  for(i in 1:length(img_names))
  {
    iloc <- which(img_names[i] == peakMat$names)[1]
    if( length(iloc) > 0  )
    {
      if( iloc > 1)
      {
        istart <- sum(peakMat$numPixels[ 1:(iloc-1) ]) + 1 
      }
      else
      {
        istart <- 1
      }
      istop <- istart + peakMat$numPixels[iloc] - 1
      lstRefac[[length(lstRefac) + 1]] <- list( name = img_names[i], pos = peakMat$pos[istart:istop,], vals = values[istart:istop] )
      
      #Apply horizontal mirror
      if(mirror_x[i])
      {
        lstRefac[[length(lstRefac)]]$pos[,"x"] <- (-1)*lstRefac[[length(lstRefac)]]$pos[,"x"] + max(lstRefac[[length(lstRefac)]]$pos[,"x"])
      }
      
      #Apply vertical mirror
      if(mirror_y[i])
      {
        lstRefac[[length(lstRefac)]]$pos[,"y"] <- (-1)*lstRefac[[length(lstRefac)]]$pos[,"y"] + max(lstRefac[[length(lstRefac)]]$pos[,"y"])
      }
      
      #Apply the rotation
      x_aux <- lstRefac[[length(lstRefac)]]$pos[,"x"]
      y_aux <- lstRefac[[length(lstRefac)]]$pos[,"y"]
      lstRefac[[length(lstRefac)]]$pos[,"x"] <- x_aux*cos(rotations[i]*pi/180) + y_aux*sin(rotations[i]*pi/180)
      lstRefac[[length(lstRefac)]]$pos[,"y"] <- -x_aux*sin(rotations[i]*pi/180) + y_aux*cos(rotations[i]*pi/180)
      
      #Offset each pos matrix to 0,0
      lstRefac[[length(lstRefac)]]$pos[,"x"] <-  lstRefac[[length(lstRefac)]]$pos[,"x"] - min(lstRefac[[length(lstRefac)]]$pos[,"x"])
      lstRefac[[length(lstRefac)]]$pos[,"y"] <-  lstRefac[[length(lstRefac)]]$pos[,"y"] - min(lstRefac[[length(lstRefac)]]$pos[,"y"])
    }
  }
  
  #Apply the raster offsets
  nextXOffset <- 0
  nextYOffset <- 0
  if(byrow)
  {
    nrow <- ceiling(length(lstRefac)/ncol)
    for(icol in 1:ncol)
    {
      ind_offsets <- seq(from = icol, by = ncol, length.out = nrow)
      ind_offsets <- ind_offsets[(ind_offsets<=length(lstRefac))]
      for(ilst in ind_offsets)
      {
        lstRefac[[ilst]]$pos[,"x"] <- lstRefac[[ilst]]$pos[,"x"] + nextXOffset
      }
      nextXOffset <- margin + max( unlist( lapply(lstRefac[ind_offsets], function(x){ x$pos[,"x"] }) ) )
    }
    for( irow in 1:nrow)
    {
      ind_offsets <- irow*ncol + (1:ncol) - ncol
      ind_offsets <- ind_offsets[(ind_offsets<=length(lstRefac))]
      for(ilst in ind_offsets)
      {
        lstRefac[[ilst]]$pos[,"y"] <- lstRefac[[ilst]]$pos[,"y"] + nextYOffset
      }
      nextYOffset <- margin + max( unlist( lapply(lstRefac[ind_offsets], function(x){ x$pos[,"y"] }) ) )
    }
  }
  else
  {
    ncol <- ceiling(length(lstRefac)/nrow)
    for(irow in 1:nrow)
    {
      ind_offsets <- seq(from = irow, by = nrow, length.out = ncol)
      ind_offsets <- ind_offsets[(ind_offsets<=length(lstRefac))]
      for(ilst in ind_offsets)
      {
        lstRefac[[ilst]]$pos[,"y"] <- lstRefac[[ilst]]$pos[,"y"] + nextYOffset
      }
      nextYOffset <- margin + max( unlist( lapply(lstRefac[ind_offsets], function(x){ x$pos[,"y"] }) ) )
    }
    for( icol in 1:ncol)
    {
      ind_offsets <- icol*nrow + (1:nrow) - nrow
      ind_offsets <- ind_offsets[(ind_offsets<=length(lstRefac))]
      for(ilst in ind_offsets)
      {
        lstRefac[[ilst]]$pos[,"x"] <- lstRefac[[ilst]]$pos[,"x"] + nextXOffset
      }
      nextXOffset <- margin + max( unlist( lapply(lstRefac[ind_offsets], function(x){ x$pos[,"x"] }) ) )
    }
  }

    #Abort if no data was selected
  if(length(lstRefac) == 0)
  {
    stop("Error: Empty data matrix, no matching image name was selected.\n")
  }
  
  #Raster matrix positon
  rasterValues <- c()
  labelsPos <- matrix( nrow = length(lstRefac), ncol = 2)
  colnames(labelsPos) <- c("x", "y")
  for( i in 1:length(lstRefac))
  {
    if(class(lstRefac[[i]]$vals) == "factor")
    {
      rasterValues <- c(rasterValues, as.character( lstRefac[[i]]$vals) )
    }
    else
    {
      rasterValues <- c(rasterValues, lstRefac[[i]]$vals) 
    }
    if( i > 1)
    {
      rasterPos <- rbind(rasterPos, lstRefac[[i]]$pos)
    }
    else
    {
      rasterPos <- lstRefac[[i]]$pos
      colnames(rasterPos) <- colnames(lstRefac[[i]]$pos)
    }
  }
  
  #Global offsets
  rasterPos[,"x"] <- rasterPos[,"x"] + 1
  rasterPos[,"y"] <- rasterPos[,"y"] + margin
  
  #Round the coordinates to ensure non-float numbers since them produce errors in many plot subsystems
  rasterPos[,"x"] <- round(rasterPos[,"x"])
  rasterPos[,"y"] <- round(rasterPos[,"y"])
  
  #Fill label positions, yes I need another loop because rasterPos is not completelly filled during previous loop
  for( i in 1:length(lstRefac))
  {
    labelsPos[i,"x"] <-  min(lstRefac[[i]]$pos[,"x"]) + 0.5*(max( lstRefac[[i]]$pos[,"x"] ) - min(lstRefac[[i]]$pos[,"x"])) 
    labelsPos[i,"y"] <-  max(rasterPos[,"y"]) - min(lstRefac[[i]]$pos[,"y"])
  }
  
  #And return the plot data
  return(list(pos = rasterPos, values = rasterValues, lab_x = labelsPos[,"x"], lab_y = labelsPos[,"y"]) )
}
