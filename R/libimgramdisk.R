#Image ramdisk API

#' Loads a part of a ff data img in RAM.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Ids Identifiers of spectra to load.
#'
#' @return a matrix containing the loaded spectra.
#'
#' Loads a part of a rMSI object in computer memory. The spectra to load is specified by its identifiers (Ids).
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

#' Loads a part of a ff data img in RAM.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a matrix containing the loaded spectra.
#'
#' Loads a part of a rMSI object in computer memory. The spectra to load is specified by its pixel coordinates (Coords).
#'
#' @export
#'
loadImgCunckFromCoords<-function(Img, Coords)
{
  return(loadImgCunckFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#' Loads a part of a ff data img in RAM.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Cube The datacube index to load.
#'
#' @return a matrix containing the loaded spectra.
#'
#' This function loads a specified ff datacubes.
#' It loads full cubes so it should be the most eficient and fast way of loading data to RAM.
#'
#' @export
#'
loadImgCunckFromCube<-function(Img, Cube)
{
  return(Img$data[[Cube]][,])
}

#' Obtain the cube index and cube row of a given image coords.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Coords a coordinates vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a list of two vectors of cube index and row index in ff data objects list corresponding to given coords spectra.
#'
#' Calculates where is located the data provided by pixel coordinates. The cube index and the cube row are returned.
#'
#' @export
#'
getCubeRowFromCoords<-function(Img, Coords)
{
  return(getCubeRowFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#' Obtain the cube index and cube row of a given image Ids.
#'
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Ids Identifiers of spectra.
#'
#' @return a list of two vectors of cube index and row index in ff data objects list corresponding to given coords spectra.
#'
#' Calculates where is located the data provided by pixel Identifiers. The cube index and the cube row are returned.
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
#' @param Img the rMSI object where the data is stored (ramdisk).
#' @param Coords a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y.
#'
#' @return a vector of identifiers corresponding to given coords.
#'
#' Calculate the pixel identifiers from the pixel coordinates.
#'
#' @export
#'
getIdsFromCoords<-function(Img, Coords)
{
  Zpos <-complex(real = Img$pos[,"x"], imaginary = Img$pos[,"y"]) #Convert positions to a complex numbers
  Ids<-unlist(sapply(Coords, function(x) { which(Zpos == x) }))
  return(Ids)
}


###TODO functions to store data in FF object
