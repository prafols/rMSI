#########################################################################
#     rMSI - R package for MSI data handling and visualization
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

#' buildImgIdVectorFromPeakMatrix.
#' 
#' Builds a integer vector containing all rMSI objects ID's accroding the peak matrix row order.
#' The resulting ID vector can be used to locate a spectrum ID in a peak matrix of multiple data sets.
#'
#' @param pkMat an rMSIproc peak matrix object.
#'
#' @return an integer vector with all ID's of each image in the pkMat object.
#' @export
#'
buildImgIdVectorFromPeakMatrix <- function(pkMat)
{
  imgIDs <- rep(NA, sum(pkMat$numPixels))
  istart <- 1
  for( i in 1:length(pkMat$numPixels))
  {
    istop <- istart + pkMat$numPixels[i] - 1
    imgIDs[istart:istop] <- 1:pkMat$numPixels[i]
    istart <- istop + 1
  }
  return(as.integer(imgIDs))
}

#' getImgIdsFromPeakMatrixRows.
#' 
#' Calculate the rMSI object ID's corresponding to a subset of row of a rMSIproc peak matrix.
#'
#' @param pkMat an rMSIproc peak matrix object.
#' @param rows a vector of peak matrix rows.
#'
#' @return a list with the rMSI object ID's and image that correpond the specified rows.
#' @export
#'
getImgIdsFromPeakMatrixRows <- function(pkMat, rows)
{
  #Sort rows assending and remove duplicates for fater performance
  rows <- sort(unique(rows), decreasing = F)
  
  id_lst <- list()
  allIds <- buildImgIdVectorFromPeakMatrix(pkMat)
  iRowStart <- 1
  for( i in 1:length(pkMat$names))
  {
    iRowStop <- iRowStart + pkMat$numPixels[i] - 1
    iSubRows<- which(rows %in% iRowStart:iRowStop )
    if( length(iSubRows) > 0)
    {
      id_lst[[i]] <- list( name = pkMat$names[i], uuid = pkMat$uuid[i], id = allIds[rows[iSubRows]] , pkMatRow = rows[iSubRows])
    }
    else
    {
      id_lst[[i]] <- list( name = pkMat$names[i], uuid = pkMat$uuid[i], id = c() , pkMatRow = c())
    }
    iRowStart <- iRowStop + 1
  }
  return(id_lst)
}

#' getPeakMatrixRowsFromImgIds.
#' 
#' Obtains a vector of rMSIproc peak matrix rows corresponding to a vector of rMSI obj ID's.
#'
#' @param pkMat an rMSIproc peak matrix object.
#' @param img_num the number of the data set object to select in pkMat.
#' @param ids a vector of rMSI obj ID's.
#'
#' @return a vector containing the rows of peak matrix that correspond to the selected rMSI object ID's.
#' @export
#'
getPeakMatrixRowsFromImgIds <- function( pkMat, img_num, ids )
{
  if(length(pkMat$numPixels) < img_num)
  {
    stop(paste("Error: pkMat does not containg", img_num, "images.\n"))
  }
  
  if(img_num > 1)
  {
    startRow <- sum(pkMat$numPixels[1:(img_num-1)]) + 1 
  }
  else
  {
    startRow <- 1
  }
  stopRow <- startRow + pkMat$numPixels[img_num] - 1
  selRows <- ids + startRow - 1
  
  if(max(selRows) > stopRow || min(selRows) < startRow)
  {
    stop("Error: Id selection is out of range.\n")
  }
  
  return(selRows)
}

#' Subsetting operator for rMSIproc peak matrices.
#'
#' @param x rMSIproc peak matrix object.
#' @param pixels the selected rows to retain in the peak matrix. 
#'               This argument can be an integer, a boolean or a character.
#' @param columns an integer vector with the peak matrix columns to retain.
#'
#' @return rMSIproc peak matrix object.
#'
#' @examples
#' #For the following example we will load an rMSIproc peak matrix in the pks variable:
#' pks <- rMSIproc::LoadPeakMatrix("/path/to/my/peak/matrix.zip")
#' 
#' #Subsetting a peak marix by pixel ID's:
#' pks1_100 <- pks[1:100, ]
#' 
#' #Subsetting a peak matrix using a boolean expression:
#' clus <- kmeans(pks$intensity/pks$normalizations$TIC, centers = 5) #perform a kmeans clustering with the whole peak matrix (TIC normalized)
#' clus2SubImg <- pks[clus$cluster==2, ] #Creat a subdataset with only pixels belonging to cluster 2
#' 
#' #Subsetting a peak matrix using image names:
#' pks_brain1 <- pks["Brain_img1", ]
#' 
#' #Subsetting to a specific column range:
#' pks_massrange <- pks[, 10:50]
#' 
#' #Subsetting by columns and rows:
#' pks_brain1 <- pks["Brain_img1", 10:50]
#'
#' @export
`[.rMSIprocPeakMatrix` <- function(x, pixels = 1:sum(x$numPixels), columns = 1:length(x$mass) )
{
  if(class(pixels) == "character")
  {
    if(any(!(pixels %in% x$names )))
    {
      stop("Invalid pixels: some of the provided names do not exist in the peak matrix.")
    }
    pixels <- (unlist(lapply(1:length(x$numPixels), function(iimg){ rep(x$names[iimg], x$numPixels[iimg])  }))) %in% pixels
    #Now pixels has been converted to logical so this is the first step, nect logical will be converter to integers
  }
  
  if(class(pixels) == "logical")
  {
    pixels <- as.integer(which(pixels))
  }
  
  if( (any(pixels) == 0) || (max(pixels) > sum(x$numPixels)) || (min(pixels) < (-sum(x$numPixels))))
  {
    stop("Error: invalid pixel ID's: out of range.")  
  }
  
  if( (any(columns) == 0) || (max(columns) > length(x$mass)) || (min(columns) < (-length(x$mass))))
  {
    stop("Error: invalid columns: out of range.")  
  }
  
  
  # Columns subsetting
  if(all(columns < 0))
  {
    columns <- unique(sort(columns, decreasing = T))
    columns <- (1:length(x$mass))[columns]
  } 
  else
  {
    columns <- unique(sort(columns)) #columns index must be always in accending order
  }
  
  
  # Pixels subsetting
  if(all(pixels < 0))
  {
    pixels <- unique(sort(pixels, decreasing = T)) 
    pixels <- (1:nrow(x$intensity))[pixels]
  } 
  else
  {
    pixels <- unique(sort(pixels)) #Pixel index must be always in accending order
  }
  
  
  x$mass <- x$mass[columns]  
  x$intensity <- x$intensity[pixels, columns]
  x$area <- x$area[pixels, columns]
  x$SNR <- x$SNR[pixels, columns]
  
  
  firstID <- 1
  originalNumPixelsSubImages <- x$numPixels
  originalNameSubImages <- x$names
  originalUUIDSubImages <- x$uuid
  x$numPixels <- c() #Start with empty images
  x$names <- c() #Start with empty images
  x$uuid <- c() #Start with empty images
  for(i in 1:length(originalNumPixelsSubImages))
  {
    lastID <- firstID + originalNumPixelsSubImages[i] - 1
    IDinSubImg <- which(firstID:lastID %in% pixels)
    if(length(IDinSubImg) > 0)
    {
      x$numPixels <- c(x$numPixels, length(IDinSubImg))
      x$names <- c(x$names, originalNameSubImages[i])
      x$uuid <- c(x$uuid, originalUUIDSubImages[i])
    }
    firstID <- lastID + 1
  }
  
  x$posMotors <- x$posMotors[pixels,]
  x$pos <- x$pos[pixels,]
  x$normalizations <- x$normalizations[pixels,]
  
  return(x)
}

#' Generic plot method for rMSIproc peak matrix.
#' @md 
#' @param x rMSIproc peak matrix object.
#' @param values the values used by the plot. 
#'               The behaviour of this parameter is controlled by the 'method' argument.
#' @param method a method used by the plot. Available options are: "mz", "values" and "clusters".
#' @param use_ggplot a boolean specifing if a ggplot2 backed must be used for plotting.
#'
#' @details 
#' This generic plot method allows to create different graphics from an rMSIproc peak matrix.
#' The plot type is controlled by the 'method' argument with the following option:\itemize{
#' \item __mz__: produce an ion map using the 'values' argument as the target m/z channel to display (this is the default behaviour).
#' \item __values__: produce an image using the values given in the 'values' argument as pixel intensities. 
#'   This method is useful to display the results of user's calculation directly over the image.
#' \item __clusters__: the 'values' argument contains a vector of integers indicating at which cluster belong each pixel in the image. 
#'   This method produce an image automatically coloured according the clusters vector. The colours codes are returned to be reused in further graphics. 
#'   }
#'
#' @examples
#' #For the following example we will load an rMSIproc peak matrix in the pks variable:
#' pks <- rMSIproc::LoadPeakMatrix("/path/to/my/peak/matrix.zip")
#' 
#' #Plot the m/z 848.7 distribution:
#' plot(pks, 848.7)
#' 
#' #Plot the TIC value of each pixel using the 'values' method:
#' plot(pks, values = pks$normalizations$TIC, method = "values")
#' 
#' #Cluster the peak matrix using kmeans and display each cluster on the image:
#' clus <- kmeans(pks$intensity/pks$normalizations$TIC, centers = 5)
#' plot(pks, clus$cluster, method = "clusters")
#' 
#' @export
plot.rMSIprocPeakMatrix <- function(x, values, method = "mz", use_ggplot=FALSE)
{
  if(method == "mz")
  {
    #TODO if mz values in a vector produce multiple plots
    values <- values[1] #Currently only one ion image is allowed
    if(use_ggplot)
    {
      return(plotPeakImageG(x, values))
    }
    else
    {
      plotPeakImage(x, values)  
    }
  }
  else if(method == "values")
  {
    if(use_ggplot)
    {
      return(plotValuesImageG(x, values))
    }
    else
    {
      plotValuesImage(x, values)
    }
  }
  else if(method == "clusters")
  {
    if(use_ggplot)
    {
      return(plotClusterImageG(x, values))
    }
    else
    {
      return(plotClusterImage(x, values))
    }
  }
  else
  {
    stop("Error: Invalid plot method. Valid methods are: \"mz\", \"values\" and \"clusters\".")
  }
}

#' Displays a summary of a peak matrix.
#'
#' @param x rMSIproc peak matrix object.
#'
#' @export
summary.rMSIprocPeakMatrix <- function(x)
{
  cat(paste("rMSIproc Peak Matrix containing"), length(x$numPixels), "images:\n")
  for(i in 1:length(x$names))
  {
    cat(paste0("- Image ", i, ": ", x$names[i], " (", x$numPixels[i], " pixels)\n"))
  }
  cat(paste0("\n",sum(x$numPixels)," pixels in total.\n"))
  cat(paste0(length(x$mass), " mass channels ranging from m/z ", round(min(x$mass), digits = 2), " to ", round(max(x$mass), digits = 2), ".\n" ))
  cat(paste0(format(object.size(pks), units = "Mb"), " of used memory.\n\n"))
}

#' Displays a summary of a peak matrix.
#'
#' @param x rMSIproc peak matrix object.
#'
#' @export
print.rMSIprocPeakMatrix <- function(x)
{
  summary(x)
}


