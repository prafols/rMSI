#########################################################################
#     rMSI - R package for MSI data handling and visualization
#     Copyright (C) 2014 Pere Rafols Soler
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

#Image ramdisk API

#' Loads a part of a img data in RAM.
#'
#' Loads a part of a rMSI object in computer memory. The spectra to load is specified by its identifiers (Ids).
#'
#' @param Img the rMSI object where the data is stored.
#' @param Ids Identifiers of spectra to load.
#'
#' @return a matrix containing the loaded spectra.
#'
#' @export
#'
loadImgChunkFromIds<-function(Img, Ids)
{
  #Avoid duplicates
  Ids <- unique(Ids)
  
  #Prepare the imzML file
  ibd_file <- path.expand(file.path(Img$data$path, paste0(Img$data$imzML$file, ".ibd")))
  if(!file.exists(ibd_file))
  {
    #TODO rise this error properly to the GUI
    stop("imzML data is not available")
  }
      
  return (Cload_imzMLSpectra(Img, Ids - 1)) #Ids-1 to translate from R-style indexting to C indexing
}


#' Stores a data matrix to a part of the img data.
#'
#' Overwrites a part of rMSI object with the provided data matrix. The part of rMSI object that is overwrited is specified by the Ids vector.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Ids Identifiers of spectra to be overwrited.
#' @param dm a matrix containing the data to overwrite the rMSI object.
#'
#' @export
#'
saveImgChunkAtIds<-function(Img, Ids, dm)
{
  #TODO re-implement using new rMSIXBin format
  
  #Avoid duplicates
  Ids <- unique(Ids)
  
  cubeRows <- getCubeRowFromIds(Img, Ids)

  # sort dm for fast wrinting (cube by cube)
  if(nrow(dm) > 1)
  {
    sortedIds <- unlist(lapply(cubeRows, function(x){x$id}))
    sorted_rows <- rep(0, length(Ids))
    for( i in 1:length(Ids))
    {
      sorted_rows[i] <- which( sortedIds[i] == Ids )
    }
    dm <- dm[ sorted_rows, ]
  }
  
  # Store data
  istart <- 1
  for( i in 1:length(cubeRows))
  {
    istop <- istart + length(cubeRows[[i]]$row) - 1
    Img$data[[cubeRows[[i]]$cube]][cubeRows[[i]]$row, ] <- dm[istart:istop , ]
    istart<-istop + 1
  }
}

#' Loads a part of the img data in RAM.
#'
#' Loads a part of a rMSI object in computer memory. The spectra to load is specified by its pixel coordinates (Coords).
#'
#' @param Img the rMSI object where the data is stored.
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a matrix containing the loaded spectra.
#'
#' @export
#'
loadImgChunkFromCoords<-function(Img, Coords)
{
  return(loadImgChunkFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#' Stores a data matrix to a part of img data.
#'
#' Overwrites a part of rMSI object with the provided data matrix. The part of rMSI object that is overwrited is specified by the Coords vector.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#' @param dm a matrix containing the data to overwrite the rMSI object.
#'
#' @export
#'
saveImgChunkAtCoords<-function(Img, Coords, dm)
{
  saveImgChunkAtIds(Img, getIdsFromCoords(Img, Coords), dm)
}

#' Obtain the image Identifiers from a given set of images coords
#'
#' Calculate the pixel identifiers from the pixel coordinates.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Coords a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a vector of identifiers corresponding to given coords.
#'
#'
#' @export
#'
getIdsFromCoords<-function(Img, Coords)
{
  Zpos <-complex(real = Img$pos[,"x"], imaginary = Img$pos[,"y"]) #Convert positions to a complex numbers
  Ids<-unlist(sapply(Coords, function(x) { which(Zpos == x) }))
  return(Ids)
}

#' Obtain the image Identifiers from a given set of images motor coords
#'
#' Calculate the pixel identifiers from the pixel motor coordinates.
#'
#' @param Img the rMSI object where the data is stored.
#' @param MotorCoords a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a vector of identifiers corresponding to given coords.
#'
#'
#' @export
#'
getIdsFromMotorCoords<-function(Img, MotorCoords)
{
  Zpos <-complex(real = Img$posMotors[,"x"], imaginary = Img$posMotors[,"y"]) #Convert positions to a complex numbers
  Ids<-unlist(sapply(MotorCoords, function(x) { which(Zpos == x) }))
  return(Ids)
}

#' getCoordsFromIds Obtain images coords from a set of Ids.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Ids from which coords will be obtained.
#'
#' @return a matrix containing the coords.
#' @export
#'
getCoordsFromIds<-function(Img, Ids)
{
  return( Img$pos[Ids, ] )
}

#' Load a image slice from a specified datacube columns.
#'
#' Loads a slice of an rMSI object into RAM. The slice is determined by the columns specified in Cols parameter.
#' The returned value is not a MS image, in order to obtain a plotable image the function builRasterImageFromCols must be used.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Cols the columns indexes from which data will be taken.
#'
#' @return a data matrix containing the image slice.
#'
#' @export
#'
loadImageSliceFromCols<-function(Img, Cols)
{
  #TODO this method is the one which actually created images! but it seams it is not creating the image directely... but new format is... see how it is used and re-implement all!
  dm <- matrix(nrow = nrow(Img$pos), ncol = length(Cols))
  ptr<-1
  for( i in 1:length(Img$data))
  {
      dm[ ptr:(ptr + nrow(Img$data[[i]]) - 1), ]<-Img$data[[i]][,Cols]
      ptr <- ptr + nrow(Img$data[[i]])
  }

  return(dm)
}

#' Overwrite an image slice of an rMSI object.
#'
#' Overwrite a part of a ramdisk in rMSI object with the provided data matrix in the specified columns.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Cols  the columns indexes at which data will be overwriten.
#' @param dm the new data matrix overwrite current rMSI slected columns.
#'
#' @export
#'
saveImageSliceAtCols<-function(Img, Cols, dm)
{
  #TODO this methos is incompatible with new data format, re-implement!
  ptr<-1
  for( i in 1:length(Img$data))
  {
    Img$data[[i]][,Cols]<- dm[ ptr:(ptr + nrow(Img$data[[i]]) - 1), ]
    ptr <- ptr + nrow(Img$data[[i]])
  }
}

#' Obtain the rMSI datacube columns corresponding to a mass.
#'
#' Obtains the rMSI datacube columns that corresponds to a mass range defined by a central mass and a tolerance.
#' The central mass and tolerance are also recalculated to fit acurately in the dataset mass axis.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Mass the slected mass.
#' @param Tolerance a tolerance expressed in daltons around the slected mass.
#'
#' @return A list with: a Cols vector containing columns of the datacube, Mass and Tolerance recalculated to correctly fit the dataset.
#'
#' @export
#'
getImageColsFromMass<-function(Img, Mass, Tolerance)
{

  l1 <-which.min( abs( Img$mass - ( Mass - Tolerance ) ) )
  l2 <-which.min( abs( Img$mass - ( Mass + Tolerance ) ) )

  l1 <- l1[1] #Store only the first element
  l2 <- l2[length(l2)] #store only the last element
  Mass <- round(mean(c(Img$mass[l2] , Img$mass[l1])), digits = 4)
  Tolerance <- round(0.5*(Img$mass[l2] - Img$mass[l1]), digits = 4)

  return(list( Cols = l1:l2, Mass = Mass, Tolerance = Tolerance ))
}

#' Load a image slice from a specified mass and tolerance.
#'
#' Loads a slice of an rMSI object into RAM. The slice is determined by a mass and a tolerance.
#' The returned value is not a MS image, in order to obtain a plotable image the function builRasterImageFromMass must be used.
#'
#' @param Img the rMSI object where the data is stored.
#' @param Mass the slected mass.
#' @param Tolerance a tolerance expressed in daltons around the slected mass.
#'
#' @return a list with data matrix containing the image slice, the used mass and tolerance.
#'
#' @export
#'
loadImageSliceFromMass<-function(Img, Mass, Tolerance)
{
  location<-getImageColsFromMass(Img, Mass, Tolerance)
  return( list( data = loadImageSliceFromCols(Img,location$Cols), Mass = location$Mass, Tolerance = location$Tolerance ) )
}


#' Overwrite an image slice of an rMSI object.
#'
#' Overwrite a part of a ramdisk in rMSI object with the provided data matrix in the specified mass range.
#' The mass range is obtained by the provided central mass and the necessary columns to fit the full data matrix in rMSI.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Mass  the central mass at which data will be overwriten.
#' @param dm the new data matrix overwrite current rMSI slected columns.
#'
#' @export
#'
saveImageSliceAtMass<-function(Img, Mass, dm)
{
  
  #TODO deprecated method for the new data format, re-implement
  central_col <- getImageColsFromMass(Img, Mass, 0)$Cols
  sel_cols <- ((central_col - floor(0.5*ncol(dm))) : (central_col + floor(0.5*ncol(dm))))

  #Remove the most distant case if the resulting vector is to long
  if(length(sel_cols) > ncol(dm))
  {
    sel_cols <- sel_cols[ -which.max( abs( Mass - Img$mass[sel_cols] ) )]
  }

  #Trim sel_cols to data limits
  sel_cols <- c(sel_cols[ sel_cols >= 1 ], seq( from = ( 1 + sel_cols[length(sel_cols)]),  to = ( sel_cols[length(sel_cols)] + length( sel_cols[ sel_cols < 1  ] )), length.out =  length( sel_cols[ sel_cols < 1  ] )))
  sel_cols <- c(seq( from = ( sel_cols[1] - length( sel_cols[ sel_cols > length(Img$mass)  ] )),  to =  (sel_cols[ 1 ] - 1 ), length.out =  length( sel_cols[ sel_cols > length(Img$mass)  ] )), sel_cols[ sel_cols <= length(Img$mass)]  )
  sel_cols <- sel_cols[ sel_cols >= 1 & sel_cols <= length(Img$mass)]

  saveImageSliceAtCols(Img, sel_cols, dm)
}

#' Build a image slice from specified datacube columns.
#'
#' Builds a image from the selected columns in the rMSI object. The image is returned arranged in a matrix containing each pixel value.
#'
#' @param Img the rMSI object where the data is stored.
#' @param IonIndex the starting ion index to construc the image.
#' @param IonCount the number of ions used to build the image.
#' @param Normalization optionally a vector of the normalization factor for each pixel.
#'
#' @return a matrix with the same size as image size containing the pixel values.
#'
#' @export
#'
builRasterImageFromCols<-function( Img, IonIndex, IonCount, Normalization = NULL)
{
  #Set no normalization if is null
  if(is.null(Normalization))
  {
    Normalization <- rep(1, nrow(Img$pos))
  }

  #Get the ion image
  zplots <- Cload_rMSIXBinIonImage(Img, IonIndex, IonCount, Normalization, parallel::detectCores())
  
  return( zplots )
}

#' Build a image slice from specified datacube masses.
#'
#' Builds a image from the selected masses in the rMSI object. The image is returned arranged in a matrix containing each pixel value.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Mass the slected mass.
#' @param Tolerance  a tolerance expressed in daltons around the slected mass.
#' @param Normalization  optionally a vector of the normalization factor for each pixel.
#'
#' @return list with a matrix with the same size as image size containing the pixel values, used mass and tolerance.
#'
#' @export
#'
builRasterImageFromMass<-function( Img, Mass, Tolerance, Normalization = NULL)
{
  location<-getImageColsFromMass(Img, Mass, Tolerance)
  return( list( pixels = builRasterImageFromCols(Img, location$Cols[1], length(location$Cols), Normalization),  Mass = location$Mass, Tolerance = location$Tolerance ) )
}

#' Inserts a image at specified Cols of a rMSI object.
#'
#' A raster image provided as a matrix is inserted at given Cols with a gaussian shape.
#' The raster_matrix has nrows as X direction.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Cols the columns indexes from which data will be inserted
#' @param raster_matrix a raster image represented as a matrix with pixel values.
#'
#' @export
#'
insertRasterImageAtCols<-function( Img, Cols, raster_matrix)
{
  #TODO  re-implement for the new data format
  dm <- matrix(0, nrow = nrow(Img$pos), ncol = length(Cols) )

  for( i in 1:nrow(dm))
  {
    dm[i ,] <- raster_matrix[ Img$pos[i, "x"], Img$pos[i, "y"] ] * dnorm( Cols, mean = mean(Cols) )
  }

  saveImageSliceAtCols(Img, Cols, dm)
}

#'  Inserts a image at specified mass of a rMSI object.
#'
#' A raster image provided as a matrix is inserted at given Mass with a gaussian shape in the given Tolerance.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Mass the slected mass.
#' @param Tolerance  a tolerance expressed in daltons around the slected mass.
#' @param raster_matrix a raster image represented as a matrix with pixel values.
#'
#' @export
#'
insertRasterImageAtMass<-function( Img, Mass, Tolerance, raster_matrix)
{
  #TODO  re-implement for the new data format
  insertRasterImageAtCols(Img, getImageColsFromMass(Img, Mass, Tolerance)$Cols, raster_matrix)
}
