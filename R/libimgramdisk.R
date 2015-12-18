#Image ramdisk API

#Loads a part of a ff data img in RAM
# Img - the image object list where the data is stored (ramdisk)
# Ids - Identifiers of spectra to load
# returns - a matrix containing the loaded spectras
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

#Loads a part of a ff data img in RAM
# Img - the image object list where the data is stored (ramdisk)
# Coords - a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y
# returns - a matrix containing the loaded spectras
loadImgCunckFromCoords<-function(Img, Coords)
{
  return(loadImgCunckFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#Loads a part of a ff data img in RAM
#This function loads a specified ff datacubes.
#It loads full cubes so it should be the most eficient and fast way of loading data to RAM
# Img - the image object list where the data is stored (ramdisk)
# Cube - The datacube index of spectra to load
# returns - a matrix containing the loaded spectras
loadImgCunckFromCube<-function(Img, Cube)
{
  return(Img$data[[Cube]][,])
}

#Obtain the cube index from a given image coords
# Img - the image object list where the data is stored (ramdisk)
# Coords - a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y
# returns - a list of two vectors of cube index and row index in ff data objects list corresponding to given coords spectra
getCubeRowFromCoords<-function(Img, Coords)
{
  return(getCubeRowFromIds(Img, getIdsFromCoords(Img, Coords)))
}

#Obtain the cube index from a given spectra identifiers
# Img - the image object list where the data is stored (ramdisk)
# Ids - Identifiers of spectra
# returns - a list of two vectors of cube index and row index in ff data objects list corresponding to given coords spectra
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

#Obtain the image Identifiers from a given set of images coords
# Img - the image object list where the data is stored (ramdisk)
# Coords - a coords vector of spectra to load represented as complex numbers where real part corresponds to X and imaginary to Y
# returns - a vector of identifiers corresponding to given coords spectra
getIdsFromCoords<-function(Img, Coords)
{
  Zpos <-complex(real = Img$pos[,"x"], imaginary = Img$pos[,"y"]) #Convert positions to a complex numbers
  Ids<-unlist(sapply(Coords, function(x) { which(Zpos == x) }))
  return(Ids)
}


###TODO functions to store data in FF object
