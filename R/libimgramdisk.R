#Image ramdisk API

#' Loads a part of a ff data img in RAM.
#'
#' Loads a part of a rMSI object in computer memory. The spectra to load is specified by its identifiers (Ids).
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Ids Identifiers of spectra to load.
#'
#' @return a matrix containing the loaded spectra.
#'
#' @export
#'
loadImgCunckFromIds<-function(Img, Ids)
{
  ind <- getCubeRowFromIds(Img, Ids)

  dM <-matrix(ncol = length(Img$mass), nrow = length(Ids))
  istart <- 1
  istop <- 0
  for(i in 1:length(ind))
  {
    istop <- istop + length(ind[[i]]$row)
    dM [istart:istop , ]<- Img$data[[ind[[i]]$cube]][ind[[i]]$row, ]
    istart<-istop + 1
  }

  return(dM)
}


#' Stores a data matrix to a part of ff data img.
#'
#' Overwrites a part of rMSI object with the provided data matrix. The part of rMSI object that is overwrited is specified by the Ids vector.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Ids Identifiers of spectra to be overwrited.
#' @param dm a matrix containing the data to overwrite the rMSI object.
#'
#' @export
#'
saveImgCunckAtIds<-function(Img, Ids, dm)
{
  ind <- getCubeRowFromIds(Img, Ids)

  istart <- 1
  istop <- 0
  for(i in 1:length(ind))
  {
    istop <- istop + length(ind[[i]]$row)
    Img$data[[ind[[i]]$cube]][ind[[i]]$row, ] <-  dm [istart:istop , ]
    istart<-istop + 1
  }
}

#' Loads a part of a ff data img in RAM.
#'
#' Loads a part of a rMSI object in computer memory. The spectra to load is specified by its pixel coordinates (Coords).
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a matrix containing the loaded spectra.
#'
#' @export
#'
loadImgCunckFromCoords<-function(Img, Coords)
{
  return(loadImgCunckFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#' Stores a data matrix to a part of ff data img.
#'
#' Overwrites a part of rMSI object with the provided data matrix. The part of rMSI object that is overwrited is specified by the Coords vector.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#' @param dm a matrix containing the data to overwrite the rMSI object.
#'
#' @export
#'
saveImgCunckAtCoords<-function(Img, Coords, dm)
{
  saveImgCunckAtIds(Img, getIdsFromCoords(Img, Coords), dm)
}

#' Loads a part of a ff data img in RAM.
#'
#' This function loads a specified ff datacubes.
#' It loads full cubes so it should be the most eficient and fast way of loading data to RAM.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Cube The datacube index to load.
#'
#' @return a matrix containing the loaded spectra.
#'
#' @export
#'
loadImgCunckFromCube<-function(Img, Cube)
{
  return(Img$data[[Cube]][,])
}

#' Save a data matrix to a whole cube in rMSI object
#'
#' @param Img the rMSI object where the data will be overwrited (ramdisk).
#' @param Cube the cube index that will be overwrited
#' @param dm the data matrix that will overwrite the selected cube
#'
#' @export
#'
saveImgCunckToCube<-function(Img, Cube, dm)
{
  Img$data[[Cube]][,] <- dm
}

#' Obtain the cube index and cube row of a given image coords.
#'
#' Calculates where is located the data provided by pixel coordinates. The cube index and the cube row are returned.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a list of two vectors of cube index and row index in ff data objects list corresponding to given coords spectra.
#'
#' @export
#'
getCubeRowFromCoords<-function(Img, Coords)
{
  return(getCubeRowFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#' Obtain the cube index and cube row of a given image Ids.
#'
#' Calculates where is located the data provided by pixel Identifiers. The cube index and the cube row are returned.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Ids Identifiers of spectra.
#'
#' @return a list of two vectors of cube index and row index in ff data objects list corresponding to given coords spectra.
#'
#' @export
#'
getCubeRowFromIds<-function(Img, Ids)
{
  max_nrow<-nrow(Img$data[[1]])
  icube<-(1+((Ids - 1) %/% max_nrow))
  irow<- (Ids - (icube -1) * max_nrow)

  cls<-list()
  cls[[1]]<- list( cube = icube[1], row = c(irow[1]))

  if(length(Ids) > 1)
  {
    for(i in 2:length(Ids))
    {
      if(icube[i] == icube[i-1] )
      {
        #Append row
        cls[[length(cls)]]$row <- c(cls[[length(cls)]]$row, irow[i])
      }
      else
      {
        #Add cube
        cls[[length(cls)+1]] <- list( cube = icube[i], row = c(irow[i]))
      }
    }
  }

  return(cls)
}

#' Obtain the image Identifiers from a given set of images coords
#'
#' Calculate the pixel identifiers from the pixel coordinates.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
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

#' Load a image slice from a specified datacube columns.
#'
#' Loads a slice of an rMSI object into RAM. The slice is determined by the columns specified in Cols parameter.
#' The returned value is not a MS image, in order to obtain a plotable image the function builRasterImageFromCols must be used.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Cols the columns indexes from which data will be taken.
#'
#' @return a data matrix containing the image slice.
#'
#' @export
#'
loadImageSliceFromCols<-function(Img, Cols)
{
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
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Cols  the columns indexes at which data will be overwriten.
#' @param dm the new data matrix overwrite current rMSI slected columns.
#'
#' @export
#'
saveImageSliceAtCols<-function(Img, Cols, dm)
{
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
#' @param Img the rMSI object where the data is stored (ramdisk).
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
#' @param Img the rMSI object where the data is stored (ramdisk).
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

#Private method to obtain pixels values from image slice data matrix dm
getPixelValuesFromImageSlice <- function(dm, method)
{
  pixel_values <- apply(dm, 1, method)
}

#' Build a image slice from specified datacube columns.
#'
#' Builds a image from the selected columns in the rMSI object. The image is returned arranged in a matrix containing each pixel value.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Cols the columns indexes from which data will be taken.
#' @param method the method used to calculate the pixel value. Can be max or mean.
#' @param Normalization optionally a vector of the normalization factor for each pixel.
#'
#' @return a matrix with the same size as image size containing the pixel values.
#'
#' @export
#'
builRasterImageFromCols<-function( Img, Cols, method = "max", Normalization = NULL)
{
  #Set the method
  if( method == "max")
  {
    fmethod <- max
  }
  else if(method == "mean")
  {
    fmethod <- mean
  }
  else
  {
    stop(paste("The specified method", method, "is invalid\n"))
  }

  #Set no normalization if is null
  if(is.null(Normalization))
  {
    Normalization <- rep(1, nrow(Img$pos))
  }

  dm <- loadImageSliceFromCols(Img, Cols)
  zplots<-matrix(0, nrow=Img$size["x"], ncol=Img$size["y"]) #Now I'm using a zero instead of NA to display a completely black background
  pixel_values<-getPixelValuesFromImageSlice(dm, fmethod)

  for( i in 1:nrow(Img$pos))
  {
    zplots[Img$pos[ i , "x" ], Img$pos[ i , "y" ]] <- pixel_values[i] * Normalization[i]
  }

  return( zplots )
}

#' Build a image slice from specified datacube masses.
#'
#' Builds a image from the selected masses in the rMSI object. The image is returned arranged in a matrix containing each pixel value.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Mass the slected mass.
#' @param Tolerance  a tolerance expressed in daltons around the slected mass.
#' @param method the method used to calculate the pixel value. Can be max or mean.
#' @param Normalization  optionally a vector of the normalization factor for each pixel.
#'
#' @return list with a matrix with the same size as image size containing the pixel values, used mass and tolerance.
#'
#' @export
#'
builRasterImageFromMass<-function( Img, Mass, Tolerance, method = "max", Normalization = NULL)
{
  location<-getImageColsFromMass(Img, Mass, Tolerance)
  return( list( pixels = builRasterImageFromCols(Img, location$Cols, method, Normalization),  Mass = location$Mass, Tolerance = location$Tolerance ) )
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
  insertRasterImageAtCols(Img, getImageColsFromMass(Img, Mass, Tolerance)$Cols, raster_matrix)
}

#Internal method to create an empty ramdisk of specified size
.CreateEmptyRamdisk<-function(num_of_columns, num_of_spectrums, ff_data_folder, max_ff_file_size_MB = 50)
{

  #Do not create to big files (> 50 MB by default)
  max_nrow <- floor((max_ff_file_size_MB*1024*1024)/(4*num_of_columns))
  max_nrow <- min(max_nrow, floor(.Machine$integer.max / (num_of_columns)))

  dataCube<-list()
  for(i in 1:ceiling(num_of_spectrums/ max_nrow))
  {
    nrows<-min(max_nrow, (num_of_spectrums - sum(unlist(lapply(dataCube, nrow)))))
    dataCube[[i]]<-ff::ff(vmode = "integer", dim = c(nrows, num_of_columns), filename = file.path(ff_data_folder, paste("ramdisk",i,".dat",sep = "")))
    names(dataCube)[i]<-paste("ramdisk",i,".dat",sep = "")
  }
  return(dataCube)
}
