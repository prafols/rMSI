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

#' Save a rMSI object to disk in a compressed .tar file.
#'
#' @param imgData a rMSI objecte with a image with a working ramdisk.
#' @param data_file Full path to hdd location to save image as a tar file.
#'
#' @export
SaveMsiData<-function(imgData, data_file)
{
  cat("Saving Image...\n")
  pt<-proc.time()

  setwd(dirname(data_file)) #Path of the resulting .tar file
  data_dir<-file.path(dirname(data_file), "ImgData")
  dir.create(data_dir)
  uuidObj<-imgData$uuid
  save(uuidObj, file = file.path(data_dir, "uuid.ImgR")) #Save UUID
  massObj<-imgData$mass
  save(massObj, file = file.path(data_dir, "mass.ImgR")) #Save mass axis
  sizeObj<-imgData$size
  save(sizeObj, file = file.path(data_dir, "size.ImgR")) #Save size Object
  posObj<-imgData$pos
  save(posObj, file = file.path(data_dir, "pos.ImgR")) #Save pos Object
  if(!is.null(imgData$posMotors))
  {
    posMotorsObj<-imgData$posMotors
    save(posMotorsObj, file = file.path(data_dir, "posMotors.ImgR")) #Save posMotors Object
  }
  meanSpcData<-imgData$mean
  save(meanSpcData, file = file.path(data_dir, "mean.SpcR")) #Save mean spectra
  resolutionObj<-imgData$pixel_size_um
  save(resolutionObj, file = file.path(data_dir, "pixel_size_um.ImgR")) #Save pixel size um Object
  if(!is.null(imgData$normalizations))
  {
    normalizationsObj <- imgData$normalizations
    save(normalizationsObj, file = file.path(data_dir, "normalizations.ImgR")) #Save normalizations Object
  }

  #Store also a vector of names of ff data in order to be able to restore it
  ffDataNames <- paste(names(imgData$data), "_ffzip", sep = "")
  save(ffDataNames, file = file.path(data_dir, "ffnames.ImgR")) #Save ff filenames Object

  ##New approach to avoid long and recursive path issues
  pb<-txtProgressBar(min = 0, max = length(ffDataNames), style = 3 )
  for(i in 1:length(ffDataNames))
  {
    setTxtProgressBar(pb, i)
    dm<-imgData$data[[i]][,]
    ffObj<-ff::ff(vmode = attr(attr(imgData$data[[i]],"physical"), "vmode"), dim = c(nrow(dm), ncol(dm)), filename =  file.path(data_dir, ffDataNames[i]))
    ffObj[,]<-dm
    ff::ffsave(ffObj , file =  file.path(data_dir, ffDataNames[i]))
    ff::delete(ffObj)
    rm(ffObj)
  }

  tar(tarfile =  file.path(getwd(), basename(data_file)), files = "ImgData")
  unlink(data_dir, recursive = T) #Remove intermediate data

  close(pb)
  pt<-proc.time() - pt
  cat(paste("Saving time:",round(pt["elapsed"], digits = 1),"seconds\n"))
}


#' Load rMSI data from a compressed tar.
#'
#' @param data_file The tar o imzML file containing the MS image in rMSI format or imzML.
#' @param restore_path Where the ramdisk will be created.
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param ff_overwrite Tell ff to overwrite or not current ramdisk files.
#' @param fun_label This is a callback function to update the progress bar dialog text.
#' @param close_signal function to be called if the loading process is aborted.
#' @param imzMLChecksum if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#'
#' @return an rMSI object pointing to ramdisk stored data
#'
#' Loads a rMSI data object from .tar compressed file or imzML format. It will be uncompressed at specified restore_path.
#' fun_progress can be NULL or a function with the following prototipe: fun_progress( currentState ). If NULL is used
#' a default command line progress bar is used.
#' This function will be called periodically to monitor the loading status. This is usefull to implement progressbars.
#' If ramdisk is already created befor calling this method the parameter ff_overwrite will control the loadin behaviour. If it is set to false (default)
#' The ramdisk will be kept and the imaged loaded imediatelly. Otherwise if is set to true, the while dataset will be reloaded from tar file.
#'
#' @export
LoadMsiData<-function(data_file, restore_path = file.path(dirname(data_file), paste("ramdisk",basename(data_file), sep = "_")) , fun_progress = NULL, ff_overwrite = F, fun_label = NULL, close_signal = NULL, imzMLChecksum = F)
{
  cat("Loading Image...\n")
  pt<-proc.time()

  #1- Check if the specified image ramdisk exists in the restore_path location
  datacube<-.FastLoad(restore_path)
  if(!is.null(datacube) && !ff_overwrite)
  {
    if(!is.null(fun_progress))
    {
      fun_progress(100)
    }

    pt<-proc.time() - pt
    cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))
    class(datacube) <- "rMSIObj"
    return(datacube)
  }

  fileExtension <- unlist(strsplit(basename(data_file), split = "\\."))
  fileExtension <- as.character(fileExtension[length(fileExtension)])
  if( fileExtension == "imzML")
  {
    return(import_imzML(data_file, ramdisk_path = restore_path, fun_progress = fun_progress, fun_text = fun_label, close_signal = close_signal, verifyChecksum = imzMLChecksum))
  }
  else if(fileExtension == "tar")
  {
    return(import_rMSItar(data_file,restore_path, fun_progress, fun_text = fun_label, close_signal = close_signal))
  }
  else
  {
    cat("The slected file is not valid.\n")
  }
  pt<-proc.time() - pt
  cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))
}

#' import_rMSItar.
#'
#' @param data_file The tar file containing the MS image in rMSI format.
#' @param restore_path Where the ramdisk will be created.
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param fun_text This is a callback function to update the label widget of loading data. See details for more information.
#' @param close_signal function to call if error.
#'
#'  Imports an rMSI data object in a tar data file
#'  It is recomanded to use rMSI::LoadMsiData directly instead of this function.
#'
#' @return  an rMSI data object.
#' @export
#'
import_rMSItar<-function(data_file, restore_path, fun_progress = NULL, fun_text = NULL, close_signal = NULL )
{
  setPbarValue<-function(progress)
  {
    setTxtProgressBar(pb, progress)
    return(T)
  }

  if(is.null(fun_progress))
  {
    pb<-txtProgressBar(min = 0, max = 100, style = 3 )
    fun_progress <- setPbarValue
  }
  else
  {
    pb<-NULL
  }

  #2 - Image is not preloaded so... the slow way
  if(!is.null(fun_text))
  {
    fun_text("Unpacking data...")
  }

  untar(tarfile = data_file, exdir = restore_path)
  img_path<-file.path(restore_path, "ImgData")


  if( file.exists(file.path(img_path, "uuid.ImgR") ))
  {
    load(file.path(img_path, "uuid.ImgR"))
  }
  else
  {
    #For compatibility with old datasets, just set no id
    uuidObj <- ""
  }
  
  load(file.path(img_path, "mass.ImgR"))
  load(file.path(img_path, "size.ImgR"))
  load(file.path(img_path, "pos.ImgR"))
  if( file.exists(file.path(img_path, "posMotors.ImgR")))
  {
    load(file.path(img_path, "posMotors.ImgR"))
  }
  else
  {
    posMotorsObj <- NULL
  }
  load(file.path(img_path, "mean.SpcR"))
  load(file.path(img_path, "ffnames.ImgR"))
  if(file.exists(file.path(img_path, "pixel_size_um.ImgR")))
  {
    load(file.path(img_path, "pixel_size_um.ImgR"))
  }
  else
  {
    print("Warning: Old image without resolution object. It is set to 9999 um by default!")
    resolutionObj <- 9999
  }
  if(file.exists(file.path(img_path, "normalizations.ImgR")))
  {
    load(file.path(img_path, "normalizations.ImgR"))
  }
  else
  {
    normalizationsObj <- NULL
  }

  if(!is.null(fun_text))
  {
    fun_text("Loading data in the ramdisk...")
  }
  spectra<-list()
  ppStep<-100/length(ffDataNames)
  pp<-0
  for(i in 1:length(ffDataNames))
  {
    ff::ffload(file.path(img_path, ffDataNames[i]), rootpath = restore_path, overwrite = T)
     #Get Hdd space back asap
    unlink(file.path(img_path, paste0(ffDataNames[i], ".RData")))
    unlink(file.path(img_path, paste0(ffDataNames[i], ".ffData")))
    spectra[[i]]<-ffObj
    names(spectra)[i]<-paste("ramdisk",i,".dat",sep = "")

    pp_ant<-pp
    pp<-pp+ppStep
    if(!is.null(fun_progress) && (round(pp) > round(pp_ant)) )
    {
      #Update progress bar
      if( !fun_progress(pp) )
      {
        return(NULL) #progress bar function must return true if the loading process is to be continued.
      }
    }
  }

  lapply(spectra, ff::open.ff)
  datacube<-list(name = basename(data_file), uuid = uuidObj, mass = massObj,  size = sizeObj,  pos = posObj, pixel_size_um = resolutionObj, mean = meanSpcData, data = spectra)
  if(!is.null(posMotorsObj))
  {
    datacube$posMotors <- posMotorsObj
  }
  if(!is.null(normalizationsObj))
  {
    datacube$normalizations <- normalizationsObj
  }

  unlink(img_path, recursive = T)

  save(datacube, file = file.path(restore_path, "datacube.RImg"))
  if(!is.null(pb))
  {
    close(pb)
  }
  if(!is.null(fun_text))
  {
    fun_text("Done")
  }
  
  class(datacube) <- "rMSIObj"
  return(datacube)
}

#Re-Loads a previously loaded image which still have the ff files on HDD
#The LoadMsiData whill produce a R objecte image file in the same directory of ramdisk in HDD
#This function will check for the R objecte image in the provided path if it exist and ff object is correct it will return with the image object
#In case it is not possible loading the image it will return NULL
.FastLoad<-function( restorePath )
{
  #1- Check if restorePath exists
  if( !dir.exists(restorePath) )
  {
    cat("\nNo ramdisk directory, it will be created\n")
    return(NULL)
  }

  #2- Check if fast loading R image object file exists
  if( !file.exists( file.path(restorePath, "datacube.RImg") ))
  {
    cat("\nNo datacube.RImg in ramdisk, new ramdisk will be created\n")
    return(NULL)
  }

  #3- Load the Image object
  load(file.path(restorePath, "datacube.RImg") )

  #4- Check the ramdisk and load it
  ramDiskExists <- is.element(T, unlist(lapply(datacube$data, function(x) { file.exists(attr(attr(x, "physical"), "filename"))  })))
  if(!ramDiskExists)
  {
    cat("\nCurrent ramdisk has been corrupted, it will be created\n")
    return(NULL)
  }
  lapply(datacube$data, ff::open.ff)

  cat("\nRamdisk has been sucessfully restored\n")
  return(datacube)
}


#'Remove an rMSI object ramdisk
#'
#' @param img an rMSI object.
#'
#' Removes a rMSI objecte ramdisk.
#'
#' @export
DeleteRamdisk<-function(img)
{
  lapply(img$data, function(x){ ff::close.ff(x) })
  ramdisk_path <- dirname(attr(attributes(img$data[[1]])$physical, "filename"))
  ramdisk_path_splited <- unlist(strsplit(ramdisk_path, "/"))
  if(.Platform$OS.type == "unix")
  {
    ramdisk_path <- "/"
  }
  else
  {
    ramdisk_path <- ""
  }
  ramdisk_path <- paste0(ramdisk_path,ramdisk_path_splited[1])
  if(length(ramdisk_path_splited) > 1)
  {
    for( i in 2:length(ramdisk_path_splited))
    {
      ramdisk_path<-file.path(ramdisk_path, ramdisk_path_splited[i])
      if( length(grep("ramdisk_", ramdisk_path_splited[i])) > 0 )
      {
        break;
      }
    }
  }
  unlink(ramdisk_path, recursive = T)
}

#' Create an empty rMSI object with defined mass axis and size.
#'
#' @param x_size the number of pixel in X direction.
#' @param y_size the number of pixel in Y direction.
#' @param mass_axis the mass axis.
#' @param pixel_resolution defined pixel size in um.
#' @param img_name the name for the image.
#' @param ramdisk_folder where ramdisk will be stored.
#' @param data_type a string determining data type used to store data.
#' @param uuid a string containing an universal unique identifier for the image. If it is not provided it will be created using a time code.
#'
#' Creates an empty rMSI object with the provided parameters. This method is usefull to implement importation of new data formats
#' and synthetic datasets to test and develop processing methods and tools.
#'
#' data_type possible values are:
#'  byte	  8 bit signed integer with NA.
#'  ubyte	  8 bit unsigned integer without NA.
#'  short   16 bit signed integer with NA.
#'  ushort  16 bit unsigned integer without NA.
#'  integer 32 bit signed integer with NA.
#'  single  32 bit float.
#'  double  64 bit float.
#'
#' @return the created rMSI object
#' @export
#'
CreateEmptyImage<-function(x_size, y_size, mass_axis, pixel_resolution, img_name = "New empty image", ramdisk_folder = getwd(), data_type = "integer", uuid = NULL)
{
  img<-list()
  img$name <- img_name
  img$mass <- mass_axis
  img$size <- c( x_size, y_size )
  names(img$size) <- c("x", "y")
  
  #Fill the UUID string
  if(is.null(uuid))
  {
    img$uuid <- format(Sys.time(), "%Y%m%d%H%M%S")
  }
  else
  {
    img$uuid <- uuid  
  }

  #Prepare the pos matrix
  img$pos <- matrix( ncol = 2, nrow = x_size*y_size )
  colnames(img$pos)<- c("x", "y")
  i <- 1
  for( xi in 1:x_size)
  {
    for( yi in 1:y_size)
    {
      img$pos[i,]<- c(xi, yi)
      i<-i+1
    }
  }

  img$pixel_size_um <-  pixel_resolution
  img$mean <- rep(0, length(mass_axis))

  #Prepare an empty datacube
  img$data<-.CreateEmptyRamdisk(length(mass_axis), nrow(img$pos), ramdisk_folder, vmode_type = data_type)

  class(img) <- "rMSIObj"
  return(img)
}

#' Create an empty rMSI object with defined mass axis and total number of pixels.
#'
#' @param num_of_pixels Total number of spectrums/pixels.
#' @param mass_axis the mass axis.
#' @param pixel_resolution defined pixel size in um.
#' @param img_name the name for the image.
#' @param ramdisk_folder where ramdisk will be stored.
#' @param data_type a string determining data type used to store data.
#' @param uuid a string containing an universal unique identifier for the image. If it is not provided it will be created using a time code.
#'
#' Creates an empty rMSI object with the provided parameters. This method is usefull to implement importation of new data formats
#' and synthetic datasets to test and develop processing methods and tools.
#' img$size is initialized with c(NA, NA) and the pos matrix with NA coords. Size and pos matrix must be filled by user.
#'
#' data_type possible values are:
#'  byte	  8 bit signed integer with NA.
#'  ubyte	  8 bit unsigned integer without NA.
#'  short   16 bit signed integer with NA.
#'  ushort  16 bit unsigned integer without NA.
#'  integer 32 bit signed integer with NA.
#'  single  32 bit float.
#'  double  64 bit float.
#'
#' @return the created rMSI object
#' @export
#'
CreateEmptyImage<-function(num_of_pixels, mass_axis, pixel_resolution, img_name = "New empty image", ramdisk_folder = getwd(), data_type = "integer", uuid = NULL)
{
  img<-list()
  img$name <- img_name
  img$mass <- mass_axis
  img$size <- c( NA, NA )
  names(img$size) <- c("x", "y")

  #Fill the UUID string
  if(is.null(uuid))
  {
    img$uuid <- format(Sys.time(), "%Y%m%d%H%M%S")
  }
  else
  {
    img$uuid <- uuid  
  }
  
  #Prepare the pos matrix
  img$pos <- matrix( NA, ncol = 2, nrow = num_of_pixels )
  colnames(img$pos)<- c("x", "y")

  img$pixel_size_um <-  pixel_resolution
  img$mean <- rep(0, length(mass_axis))

  #Prepare an empty datacube
  img$data<-.CreateEmptyRamdisk(length(mass_axis), nrow(img$pos), ramdisk_folder, vmode_type = data_type)
  class(img) <- "rMSIObj"
  return(img)
}

#' AverageSpectrum.
#' Computes the average spectrum of a whole rMSI image.
#'
#' @param img the rMSI image object.
#'
#' @return A vector with the average spectrum intensities. Masses are the same as the rMSI object.
#' @export
#'
AverageSpectrum <- function(img)
{
  cat("Calculating Average Spectrum...\n")
  if(!PackageChecker("rMSIproc","0.1"))
    {
    pbavg <- txtProgressBar(min = 0, max = length(img$data), style = 3)
    avgI <- rep(0, length(img$mass))
    for( i in 1:length(img$data))
    {
      setTxtProgressBar(pbavg, i)
      avgI <- avgI +  colSums(img$data[[i]][,])
    }
    avgI <- avgI/nrow(img$pos)
    close(pbavg)
    return(avgI)
  }
  return(rMSIproc::MTAverageSpectrum(img))
}

#' BaseSpectrum.
#' Computes the base spectrum of a whole rMSI image where the intensity value
#' for each mass bin is calculated as the maxium of all mass channels.
#'
#' @param img the rMSI image object.
#'
#' @return A vector with the base spectrum intensities. Masses are the same as the rMSI object.
#' @export
#'
BaseSpectrum <- function(img)
{
  cat("Calculating Base Spectrum...\n")
  pb <- txtProgressBar(min = 0, max = length(img$data), style = 3)
  maxSub <- rep(0, length(img$mass))
  for( i in 1:length(img$data))
  {
    setTxtProgressBar(pb, i)
    maxSub <- maxSub + apply(img$data[[i]][,], 2, max)
  }
  close(pb)
  return(maxSub)
}

#' SortIDsByAcquisition: Order the rMSI pixel IDs according the acqusition sequence (first acquired pixel is the first one).
#'
#' @param img a rMSI image.
#'
#' @return a vector of ID's ordered acording acquisiton.
#' @export
#'
SortIDsByAcquisition <- function(img)
{
  idxArray <- matrix( c(1:nrow(img$pos), img$pos[,"x"], img$pos[,"y"]), nrow = nrow(img$pos), ncol = 3, byrow = F )
  colnames(idxArray) <- c("id", "x", "y")
  idxArray <- idxArray[order(idxArray[,"y"], idxArray[,"x"]), ]
  return(idxArray[,"id"])
}


#' PlotClusterImage.
#' Plot a segmentation image with the user-given clusters.
#'
#' @param posMat a two columns matrix where first column ara the x coodrinates of values and second column the y coordinates.
#' @param clusters a vector with integer number according the cluster of each pixel.
#' @param rotate rotation to apply.
#' @param pixel_size_um the pixel resolution in um.
#'
#' @return a vector with the used color for each cluster sorted according clustering numering in assending order.
#' @export
#'
PlotClusterImage <- function( posMat, clusters,  rotate = 0,  pixel_size_um = 100 )
{
  img <- list()
  img$pos <- posMat
  colnames(img$pos) <- c("x", "y")
  img$size <- c(max(img$pos[ ,"x"]), max(img$pos[ ,"y"]))
  names(img$size) <- c("x", "y")
  img$pixel_size_um <- pixel_size_um

  #Prepare image matrix
  zplots<-matrix(0, nrow=img$size["x"], ncol=img$size["y"]) #Now I'm using a zero instead of NA to display a completely black background
  for( i in 1:nrow(img$pos))
  {
    zplots[img$pos[ i , "x" ], img$pos[ i , "y" ]] <- clusters[i]
  }

  #Create the raster
  my_raster <- raster::raster( nrow = ncol(zplots), ncol = nrow(zplots), xmn= 0, xmx= nrow(zplots), ymn= 0, ymx= ncol(zplots))
  raster::values(my_raster) <- as.vector(zplots)
  rm(zplots)

  #Put zplots matrix and some metadata in a list
  img_sgn <- list(raster = my_raster, mass = "", tolerance = 0, cal_resolution = img$pixel_size_um)
  rm(my_raster)

  raster_RGB<-.BuildSingleIonRGBImage(img_sgn, XResLevel = 3, light = 5)
  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotate, display_axes = F, display_syscoords = F)

  #Get the colors used for clusters as a plotable form (RGB code) using the same rMSI internal function as raster image
  colras <- raster::raster( nrow = 1, ncol = 1+length(unique(clusters)))
  raster::values(colras) <- sort(c(0, unique(clusters)))
  rgbColRas <- .ReMappingIntensity2HSV(colras, value_multiplier = 5)
  clusterColors <- c()
  Rchannel <- raster::values(rgbColRas$layer.1)
  Gchannel <- raster::values(rgbColRas$layer.2)
  Bchannel <- raster::values(rgbColRas$layer.3)
  for( i in 1:length(unique(clusters))) #I'm avoiding the fist values because is the zero used to draw the background
  {
    clusterColors <- c(clusterColors, rgb( Rchannel[i + 1], Gchannel[i + 1], Bchannel[i + 1], 255, maxColorValue = 255))
  }

  return(clusterColors)
}


#' PlotValues.
#'
#' Plot values in a image using the same methods as plotting an MS image.
#' The raster position of each value is definied in posMat.
#'
#' @param posMat a two columns matrix where first column ara the x coodrinates of values and second column the y coordinates.
#' @param values the values to plot.
#' @param rotate rotation to apply.
#' @param scale_title a text label for the color scale.
#' @param pixel_size_um the pixel resolution in um.
#'
#' @export
#'
PlotValues <- function(posMat, values, rotate = 0, scale_title = "", pixel_size_um = 100)
{
  fooImg <- list()
  fooImg$pos <- posMat
  colnames(fooImg$pos) <- c("x", "y")
  fooImg$size <- c(max(fooImg$pos[ ,"x"]), max(fooImg$pos[ ,"y"]))
  names(fooImg$size) <- c("x", "y")
  fooImg$pixel_size_um <- pixel_size_um
  PlotTICImage( fooImg, values, rotate, scale_title )
}

#' PlotTICImage.
#'
#' @param img an rMSI object.
#' @param TICs a vector of TIC values ordered acording pos array in img object.
#'
#' @export
#'
PlotTICImage <- function(img, TICs = NULL, rotate = 0, scale_title = "TIC")
{
  #Calculate TICs
  if(is.null(TICs))
  {
    pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
    setTxtProgressBar(pb, 0)
    TICs <- c()
    for( i in 1:length(img$data))
    {
      TICs <- c(TICs, rowSums(img$data[[i]][,]))
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  #Do not plot infinites
  TICs[which(is.infinite(TICs))] <- 0

  #Prepare image matrix
  zplots<-matrix(0, nrow=img$size["x"], ncol=img$size["y"]) #Now I'm using a zero instead of NA to display a completely black background
  for( i in 1:nrow(img$pos))
  {
    zplots[img$pos[ i , "x" ], img$pos[ i , "y" ]] <- TICs[i]
  }

  #Create the raster
  my_raster <- raster::raster( nrow = ncol(zplots), ncol = nrow(zplots), xmn= 0, xmx= nrow(zplots), ymn= 0, ymx= ncol(zplots))
  raster::values(my_raster) <- as.vector(zplots)
  rm(zplots)

  #Put zplots matrix and some metadata in a list
  img_sgn <- list(raster = my_raster, mass = scale_title, tolerance = 0, cal_resolution = img$pixel_size_um)
  rm(my_raster)

  raster_RGB<-.BuildSingleIonRGBImage(img_sgn, XResLevel = 3, light = 5)

  layout( matrix( 2:1, ncol = 2, nrow = 1, byrow = TRUE ), widths = c(7, rep(1, 2)) )
  .plotIntensityScale(img_sgn, light =  5)
  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotate, display_axes = F, display_syscoords = F)
}

#' ConvertrMSIimg2Bin.
#'
#' @param img a rMSI image object to be converted.
#' @param out_data_dir the resulting output data directory.
#'
#' @return nothing.
#' @export
#'
ConvertrMSIimg2Bin <- function( img, out_data_dir)
{

  dir.create(out_data_dir)

  #Export global m/z axis to txt
  write(img$mass, file.path(out_data_dir, "mz.txt" ), ncolumns = 1)

  #Export metadata as a plain ascii file:
  img_name_char <- paste("Original image name:\t", img$name, "\n", sep = "")
  img_size_charX <- paste("Image size in X:\t", img$size["x"], "\n", sep = "")
  img_size_charY <- paste("Image size in Y:\t", img$size["y"], "\n", sep = "")
  pixel_size_char <- paste("Pixel size (um):\t", img$pixel_size_um, "\n", sep = "")
  write( paste(img_name_char, img_size_charX, img_size_charY, pixel_size_char, sep =""), file= file.path(out_data_dir, "metadata.txt")  )

  #Export spectra intensity data
  spc_id <- 1
  subdir <- 1
  dir.create(file.path(out_data_dir, subdir))
  spc_in_subdir_count <- 0
  for( ic in 1:length(img$data))
  {
    cat(paste("Processing cube", ic, "of", length(img$data), "\n"))
    cube <- loadImgChunkFromCube(img, ic)

    for( i in 1:nrow(cube))
    {
      if(spc_in_subdir_count >= 1e4)
      {
        subdir <- subdir + 1
        dir.create(file.path(out_data_dir, subdir))
        spc_in_subdir_count <- 0
      }

      fname <- paste ("ID",spc_id, "_X", img$pos[spc_id, "x"] , "_Y", img$pos[spc_id, "y"], ".bin", sep = "")
      spc_id <- spc_id + 1
      spc <- cube[i, ]
      writeBin(spc, file.path(out_data_dir, subdir, fname ), size = 4, endian="little")
      spc_in_subdir_count <- spc_in_subdir_count + 1
    }

  }
}

#' ConvertBin2rMSIimg.
#'
#' @param in_data_dir data dir where the bin image is located.
#' @param out_img_tar_file if not NULL the imported image will be also stored as a tar file.
#'
#' @return the rMSI image object.
#' @export
#'
ConvertBin2rMSIimg <- function( in_data_dir, out_img_tar_file = NULL )
{
  #Get the metadata
  mz_axis <- as.numeric(read.table(file.path(in_data_dir, "mz.txt"))[,1])
  metadata <- read.table(file.path(in_data_dir, "metadata.txt"), sep = "\t", colClasses = "character")
  imgName <- as.character(metadata[1,2])
  xSize <- as.numeric(metadata[2,2])
  ySize <- as.numeric(metadata[3,2])
  resolution <- as.numeric(metadata[4,2])

  #List bin files
  bin_files<-list.files(path=in_data_dir, include.dirs = F, recursive = T, pattern = "*.bin", full.names = T)

  #Fill image fields
  img <-CreateEmptyImage( num_of_pixels = length(bin_files) , pixel_resolution = resolution, img_name =  imgName, mass_axis = mz_axis )
  img$size["x"] <- xSize
  img$size["y"] <- ySize

  #Prepare a dataframe with each bin file info
  print("Parsing bin file names...")
  pb<-txtProgressBar(min = 0, max = length(bin_files), style = 3 )
  pb_i <- 0
  ID <- c()

  for( bin_file in bin_files)
  {
    pb_i <- pb_i + 1
    setTxtProgressBar(pb, pb_i)
    pixel_fields <- strsplit(as.character(strsplit(basename(bin_file), split = "\\.")[[1]])[1], split = "_")[[1]]
    id <- as.numeric(strsplit(pixel_fields[1], split = "ID")[[1]])[2]
    X_cord <- as.numeric(strsplit(pixel_fields[2], split = "X")[[1]])[2]
    Y_cord <- as.numeric(strsplit(pixel_fields[3], split = "Y")[[1]])[2]
    img$pos[id, "x"] <- X_cord
    img$pos[id, "y"] <- Y_cord
    ID <- c(ID, id)
  }
  data_inf <- data.frame( ID, bin_files )
  data_inf <- data_inf[order(ID), ] #Sort by ID's (faster datacubes writing)
  close(pb)

  #Extract bin files
  print("Extracting spectra from bin files...")
  pb<-txtProgressBar(min = 0, max = nrow(data_inf), style = 3 )
  dm <- matrix(nrow = 100, ncol = length( mz_axis )) #Read HDD in chunks of 100 spectra
  partial_ids <- c()
  dm_irow <- 1
  for( i in 1:nrow(data_inf))
  {
    setTxtProgressBar(pb, i)
    dm[dm_irow, ] <- as.integer(readBin(as.character(data_inf[i, "bin_files"]), integer(),n=length(mz_axis) ,size = 4, signed=T, endian="little" ))
    partial_ids <- c(partial_ids, data_inf[i, "ID"])
    dm_irow <- dm_irow + 1

    if( dm_irow > nrow(dm) || i == nrow(data_inf))
    {
      saveImgChunkAtIds( img,  partial_ids, dm[1:length(partial_ids), ] )
      partial_ids <- c()
      dm_irow <- 1
    }
  }
  close(pb)

  print("Calculating average spectrum...")
  img$mean <- AverageSpectrum(img)

  if(!is.null(out_img_tar_file))
  {
    SaveMsiData(img, out_img_tar_file)
  }
  return(img)
}

#' plotParamAcqOrdered.
#'
#' @param img rMSI object from wich data must be ploted.
#' @param Param a vector of elements to plot.
#' @param yAxisLabel.
#'
#' Param will be ploted ordered according the order of pixels in MALDI acqusition.
#'
#' @export
#'
plotParamAcqOrdered <- function( img, Param, yAxisLabel = "Param" )
{
  idxArray <- matrix( c(1:nrow(img$pos), img$pos[,"x"], img$pos[,"y"]), nrow = nrow(img$pos), ncol = 3, byrow = F )
  colnames(idxArray) <- c("id", "x", "y")
  idxArray <- idxArray[order(idxArray[,"y"], idxArray[,"x"]), ]
  plot(Param[idxArray[,"id"]], type="l", col ="red", ylab = yAxisLabel, xlab = "Pixel" )
}



#' remap2ImageCoords.
#'
#' @param dataPos a pos matrix as it is in rMSI data object.
#'
#' This function should be only used to implement data importers from foreign formats.
#' This functions maps a MALDI motors coors space to a image coord space.
#' dataPos matrix is a two columns matrix where first column stores x positions and second y pixel positions.
#' a remapped dataPos matrix do not contain empty raster positions neighter offsets.
#'
#' @return the dataPos matrix remapped.
#' @export
#'
remap2ImageCoords <- function(dataPos)
{
  #1- Calc offsets and subtract it
  x_offset<-min(dataPos[,"x"])
  y_offset<-min(dataPos[,"y"])
  for(i in 1:nrow(dataPos))
  {
    dataPos[i, "x"] <- dataPos[i, "x"] - x_offset + 1
    dataPos[i, "y"] <- dataPos[i, "y"] - y_offset + 1
  }

  #2- Compute Motor coords range
  x_size<-max(dataPos[,"x"])
  y_size<-max(dataPos[,"y"])

  #3- Map MALDI motor coords to image cords (1-pixels steps)
  #It is important to map MALDI motors coords to image coords.
  #Otherwise, null extra pixels may be added leading to bad reconstruction
  px_map <- matrix( 0, nrow = x_size, ncol = y_size)
  for(i in 1:nrow(dataPos))
  {
    xi <- dataPos[i, "x"]
    yi <- dataPos[i, "y"]
    px_map[xi, yi]<- i
  }

  colNull <- which( base::colSums(px_map) == 0)
  rowNull <- which( base::rowSums(px_map) == 0)
  remap<-FALSE
  if( length(colNull) > 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) > 0 && length(rowNull) == 0 )
  {
    px_map_ <- px_map[ , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) == 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , ]
    remap<-TRUE
  }

  if(remap)
  {
    for(ix in 1:nrow(px_map_))
    {
      for(iy in 1:ncol(px_map_))
      {
        if(px_map_[ix, iy] > 0)
        {
          dataPos[px_map_[ix, iy], "x"] <- ix
          dataPos[px_map_[ix, iy], "y"] <- iy
        }
      }
    }
  }

  return(dataPos)
}

#' .controlled_loadAbort.
#'
#' @param text to be promt at R console
#' @param close_signal the function to call
#'
.controlled_loadAbort <- function(text, close_signal = NULL)
{
  if(!is.null(close_signal))
  {
    close_signal()
  }
  stop(text)
}

#' ParseBrukerXML.
#'
#' Reads a Bruker's xml file exported using fleximaging.
#' A list is returned where each element in the list is named according the ROI name.
#' Each element in the list consists in a data.frame with the pixels XY coordinates inside each ROI.
#'
#' @param xml_path the full path where XML file is stored.
#'
#' @return ROI pixel coordinates arranged in a named list.
#' @export
#' 
ParseBrukerXML <- function(xml_path)
{
  roilst <- CparseBrukerXML(path.expand(xml_path))
  
  if( !is.null(roilst$Error))
  {
    stop(roilst$Error)
  }
  
  return(roilst)
}

#' ReadBrukerRoiXML.
#' 
#' Reads a Bruker ROI XML file and matches it to an rMSI image object.
#' A list of rMSI ID's contained in each Bruker's ROI will be returned.
#'
#' @param img an rMSI MS image object.
#' @param xml_file a full path to a Bruker XML ROI file.
#'
#' @return a list containing lists of ID's for each ROI.
#' @export
#'
ReadBrukerRoiXML <- function(img, xml_file)
{
  if( is.null(img$posMotors))
  {
    stop("ERROR: image posMotros matrix not available.\nYou need to re-import MS data using a recent verison of rMSI.\n")
  }
  
  spectraListRois <- ParseBrukerXML(xml_file)
  lstRois <- list()
  imPosMat <- complex( real = img$posMotors[, "x"], imaginary = img$posMotors[, "y"])
  
  for( i in 1:length(spectraListRois))
  {
    cat(paste0("Parsing ROI ", i, " of ", length(spectraListRois), "\n"))
    lstRois[[i]] <- list(name = spectraListRois[[i]]$name, id = c())
    for( j in 1:nrow(spectraListRois[[i]]$pos))
    {
      #Extract original X Y Bruker Coords
      imCoord <- complex(real = spectraListRois[[i]]$pos$x[j], imaginary = spectraListRois[[i]]$pos$y[j])
      
      #Look for this Bruker coords in image pos matrix
      matchXY <- which(imPosMat == imCoord)
      if( length(matchXY) > 0)
      {
        lstRois[[i]]$id <- c( lstRois[[i]]$id, matchXY[1])
        if(  length(matchXY) > 1 )
        {
          cat(paste0("WARNING: roi ",spectraList$name, " coordinates x", Re(imCoord), " , y", Im(imCoord), " are duplicated.\n" ))
        }
      }
      else
      {
        cat(paste0("WARNING: roi ",spectraList$name, " coordinates x", Re(imCoord), " , y", Im(imCoord), " not found in image.\n" ))
      }
      
    }
  }
  
  return(lstRois)
}


#' CreateSubDataset.
#' 
#' Creates a new rMSI image object from sub-set of selected pixels by ID's.
#'
#' @param img the original rMSI object.
#' @param id a vector of ID's to retain in the sub data set.
#' @param ramdisk_path a full disk path where the new ramdisk will be stored.
#' @param new_mass a new mass axis if resampling must be used.
#'
#' @return the sub rMSI object.
#' @export
#'
CreateSubDataset <- function(img, id, ramdisk_path, new_mass = img$mass)
{
  cat("Creating the sub image...\n")

  if(!dir.exists(ramdisk_path))
  {
    dir.create(ramdisk_path, showWarnings = F, recursive = T)
  }
  
  #Resample data only if the supplied mass axis is different
  bResampleData <- !identical(img$mass, new_mass)
  
  subImg <- CreateEmptyImage(num_of_pixels = length(id), 
                             mass_axis = new_mass, 
                             pixel_resolution = img$pixel_size_um, 
                             img_name = paste0(img$name, "_sub"), 
                             ramdisk_folder = ramdisk_path, 
                             data_type = attr(attr(img$data[[1]], "physical"), "vmode"), 
                             uuid = img$uuid )
  
  subImg$pos <- remap2ImageCoords( img$pos[id, ] )
  
  if(!is.null(img$posMotors))
  {
    subImg$posMotors <- img$posMotors[id, ]
  }
  if( !is.null(img$normalizations))
  {
    subImg$normalizations <- img$normalizations
    for( i in 1:length(subImg$normalizations))
    {
      subImg$normalizations[[i]] <- subImg$normalizations[[i]][id]
    }
  }
  
  #Copy data on given id list
  pb <- txtProgressBar(min = 0, max = length(subImg$data), style = 3)
  istart <- 1
  for( i in 1:length(subImg$data))
  {
    setTxtProgressBar(pb, i)
    subID <- id[ istart:(istart + nrow(subImg$data[[i]]) - 1) ]
    istart <- istart +  nrow(subImg$data[[i]])
    dm <- loadImgChunkFromIds(img, subID)
    
    #Resampling...
    if(bResampleData)
    {
      dmSub <- matrix(nrow = nrow(dm), ncol = length(subImg$mass))
      for( irow in 1:nrow(dm))
      {
        dmSub[irow, ] <- approx(img$mass, dm[irow,], xout = subImg$mass, ties = "ordered", yleft = 0, yright = 0)$y
      }
    }
    else
    {
      dmSub <- dm
    }
    
    subImg$data[[i]][ , ] <- dmSub
  }
  close(pb)
  
  subImg$size <- c( max( subImg$pos[, "x"] ), max( subImg$pos[, "y"] ))
  names(subImg$size) <- c("x", "y")
  
  subImg$mean <- AverageSpectrum(subImg)
  return(subImg)
}

#' PackageChecker.
#' 
#' Cheks if there is the specified package in the correct version in the user library.
#'
#' @param PackageName the name of the package to be checked.
#' @param PackageVersion the minimum version of the package required.
#'
#' @return boolean.
#' @export
#'
PackageChecker <- function(PackageName,PackageVersion)
{
  info<-sessionInfo()
  
  if (length(info$otherPkgs)==0)  #Looking for other packages
  {
    return(FALSE)
  }
  
  for (a in 1:length(info$otherPkgs)) 
  {
    if ((names(info$otherPkgs[a]) == PackageName) & 
        (info$otherPkgs[[a]]$Version >= PackageVersion))
    {
      return(TRUE)
    }
  }
  return(FALSE)
}


