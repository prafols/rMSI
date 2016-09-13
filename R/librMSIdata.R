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
  massObj<-imgData$mass
  save(massObj, file = file.path(data_dir, "mass.ImgR")) #Save mass axis
  sizeObj<-imgData$size
  save(sizeObj, file = file.path(data_dir, "size.ImgR")) #Save size Object
  posObj<-imgData$pos
  save(posObj, file = file.path(data_dir, "pos.ImgR")) #Save pos Object
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
    ffObj<-ff::ff(vmode = "integer", dim = c(nrow(dm), ncol(dm)), filename =  file.path(data_dir, ffDataNames[i]))
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
#' @param data_file The tar file containing the MS image in rMSI format.
#' @param restore_path Where the ramdisk will be created.
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param ff_overwrite Tell ff to overwrite or not current ramdisk files.
#'
#' @return an rMSI object pointing to ramdisk stored data
#'
#' Loads a rMSI data object from .tar compressed file. It will be uncompressed at specified restore_path.
#' fun_progress can be NULL or a function with the following prototipe: fun_progress( currentState ). If NULL is used
#' a default command line progress bar is used.
#' This function will be called periodically to monitor the loading status. This is usefull to implement progressbars.
#' If ramdisk is already created befor calling this method the parameter ff_overwrite will control the loadin behaviour. If it is set to false (default)
#' The ramdisk will be kept and the imaged loaded imediatelly. Otherwise if is set to true, the while dataset will be reloaded from tar file.
#'
#' @export
LoadMsiData<-function(data_file, restore_path = file.path(dirname(data_file), paste("ramdisk",basename(data_file), sep = "_")) , fun_progress = NULL, ff_overwrite = F)
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
    return(datacube)
  }

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
  untar(tarfile = data_file)
  img_path<-file.path(getwd(), "ImgData")

  load(file.path(img_path, "mass.ImgR"))
  load(file.path(img_path, "size.ImgR"))
  load(file.path(img_path, "pos.ImgR"))
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


  spectra<-list()
  ppStep<-100/length(ffDataNames)
  pp<-0
  for(i in 1:length(ffDataNames))
  {
    ff::ffload(file.path(img_path, ffDataNames[i]), rootpath = restore_path, overwrite = T)
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
  datacube<-list(name = basename(data_file), mass = massObj,  size = sizeObj,  pos = posObj, pixel_size_um = resolutionObj, mean = meanSpcData, data = spectra)
  if(!is.null(normalizationsObj))
  {
    datacube$normalizations <- normalizationsObj
  }

  unlink("ImgData", recursive = T)

  save(datacube, file = file.path(restore_path, "datacube.RImg"))
  if(!is.null(pb))
  {
    close(pb)
  }
  pt<-proc.time() - pt
  cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))
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
  ramdisk_path <- dirname(attr(attributes(img$data[[1]])$physical, "filename"))
  ramdisk_path_splited <- unlist(strsplit(ramdisk_path, "/"))
  ramdisk_path <- ""
  for( i in 1:length(ramdisk_path_splited))
  {
    ramdisk_path<-file.path(ramdisk_path, ramdisk_path_splited[i])
    if( length(grep("ramdisk_", ramdisk_path_splited[i])) > 0 )
    {
      break;
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
#'
#' Creates an empty rMSI object with the provided parameters. This method is usefull to implement importation of new data formats
#' and synthetic datasets to test and develop processing methods and tools.
#'
#' @return the created rMSI object
#' @export
#'
CreateEmptyImage<-function(x_size, y_size, mass_axis, pixel_resolution, img_name = "New empty image", ramdisk_folder = getwd())
{
  img<-list()
  img$name <- img_name
  img$mass <- mass_axis
  img$size <- c( x_size, y_size )
  names(img$size) <- c("x", "y")

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
  img$data<-.CreateEmptyRamdisk(length(mass_axis), nrow(img$pos), ramdisk_folder)

  return(img)
}

#' Create an empty rMSI object with defined mass axis and total number of pixels.
#'
#' @param num_of_pixels Total number of spectrums/pixels.
#' @param mass_axis the mass axis.
#' @param pixel_resolution defined pixel size in um.
#' @param img_name the name for the image.
#' @param ramdisk_folder where ramdisk will be stored.
#'
#' Creates an empty rMSI object with the provided parameters. This method is usefull to implement importation of new data formats
#' and synthetic datasets to test and develop processing methods and tools.
#' img$size is initialized with c(NA, NA) and the pos matrix with NA coords. Size and pos matrix must be filled by user.
#'
#' @return the created rMSI object
#' @export
#'
CreateEmptyImage<-function(num_of_pixels, mass_axis, pixel_resolution, img_name = "New empty image", ramdisk_folder = getwd())
{
  img<-list()
  img$name <- img_name
  img$mass <- mass_axis
  img$size <- c( NA, NA )
  names(img$size) <- c("x", "y")

  #Prepare the pos matrix
  img$pos <- matrix( NA, ncol = 2, nrow = num_of_pixels )
  colnames(img$pos)<- c("x", "y")

  img$pixel_size_um <-  pixel_resolution
  img$mean <- rep(0, length(mass_axis))

  #Prepare an empty datacube
  img$data<-.CreateEmptyRamdisk(length(mass_axis), nrow(img$pos), ramdisk_folder)

  return(img)
}

#' Computes the average spectrum of a whole rMSI image.
#'
#' @param img the rMSI image object.
#'
#' @return A vector with the average spectrum intensities. Masses are the same as the rMSI object.
#' @export
#'
AverageSpectrum <- function(img)
{
  avgI<-apply(matrix(unlist(lapply(img$data, function(x){ ff::ffrowapply(colSums(x[i1:i2,,drop=FALSE]), X=x, RETURN = TRUE, CFUN = "csum", FF_RETURN = FALSE) })), nrow = length(img$data), byrow = T), 2, sum)
  avgI<-avgI/( sum(unlist(lapply(img$data, nrow))) )

  return(avgI)
}

#' NormalizeTIC: Calculates the TIC normalizatin of each pixel as the sum of all intensities.
#'
#' @param img the rMSI image object.
#' @param remove_empty_pixels boolean detailing if pixels detected to not contain data must be removed from normalization (smaller than mean-*sd).
#'
#' @return  a rMSI image containing the normalizations$TIC field or TICne if remove_empty_pixels is true.
#' @export
#'
NormalizeTIC <- function(img, remove_empty_pixels = FALSE)
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

  if(remove_empty_pixels)
  {
    minallowedTIC <- mean(TICs) - sd(TICs)
    TICs[which(TICs < minallowedTIC)] <- Inf #Remove pixels divinding them by infinite
    img <- AppendNormalizationCoefs(img, "TICne", TICs)
  }
  else
  {
    img <- AppendNormalizationCoefs(img, "TIC", TICs)
  }
  return(img)
}

#' NormalizeMAX: Calculates the MAX normalizatin of each pixel as the maximum of all intensities.
#'
#' @param img the rMSI image object.
#' @param remove_empty_pixels boolean detailing if pixels detected to not contain data must be removed from normalization (smaller than mean-*sd).
#'
#' @return  a rMSI image containing the normalizations$MAX field or MAXne if remove_empty_pixels is true..
#' @export
#'
NormalizeMAX <- function(img, remove_empty_pixels = FALSE)
{
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)
  MAXs <- c()
  for( i in 1:length(img$data))
  {
    MAXs <- c(MAXs, apply(img$data[[i]][,], 1, max))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  if(remove_empty_pixels)
  {
    minallowedMAX <- mean(MAXs) - sd(MAXs)
    MAXs[which(MAXs < minallowedMAX)] <- Inf #Remove pixels divinding them by infinite
    img <- AppendNormalizationCoefs(img, "MAXne", MAXs)
  }
  else
  {
    img <- AppendNormalizationCoefs(img, "MAX", MAXs)
  }
  return(img)
}

#' NormalizeByAcqusitionDegradation: Normalizes an rMSI image to compensate ionization source degradation.
#'
#' @param img  the rMSI object to normalize.
#' @param winSize the window size use for smoothing (0 to 1).
#'
#' @return a rMSI image containing the normalizations$AcqTic field.
#' @export
#'
NormalizeByAcqusitionDegradation <- function( img, winSize = 0.1 )
{
  #Calc all tics
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)
  TICs <- c()
  for( i in 1:length(img$data))
  {
    TICs <- c(TICs, rowSums(img$data[[i]][,]))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #Sort TICs according acquisition
  idxArray <- matrix( c(1:nrow(img$pos), img$pos[,"x"], img$pos[,"y"]), nrow = nrow(img$pos), ncol = 3, byrow = F )
  colnames(idxArray) <- c("id", "x", "y")
  idxArray <- idxArray[order(idxArray[,"y"], idxArray[,"x"]), ]
  TICs <- matrix( c( idxArray[,"id"], TICs[idxArray[,"id"]]), ncol = 2)
  colnames(TICs) <- c("id", "TIC")

  smTICs <- TICs
  winLen <- floor(0.5*nrow(TICs) * winSize)
  for( i in 1:nrow(TICs))
  {
    iStart <- i - winLen
    iStop <- i + winLen
    if( iStart < 1)
    {
      iStop <- iStop + (1 - iStart)
      iStart <- 1
    }
    if( iStop > nrow(TICs))
    {
      iStart <- iStart + (nrow(TICs) - iStop)
      iStop <- nrow(TICs)
    }

    discardLows <- mean(TICs[iStart:iStop, "TIC"]) - sd(TICs[iStart:iStop, "TIC"])
    dWind <- TICs[iStart:iStop, "TIC"]
    dWind <- dWind[which(dWind > discardLows)]
    smTICs[i, "TIC"] <- mean(dWind)
    #TODO pensar que fer si es NAN
  }

  #Order smTICs according Id's and append it to image normalization
  smTICs <- smTICs[order(smTICs[,"id"]), ]
  img <- AppendNormalizationCoefs(img, "AcqTic", smTICs[,"TIC"])
  return(img)
}

#' AppendNormalizationCoefs: Appends a new normalization coefinients to rMSI object.
#'
#' @param img the rMSI object to append normalization coeficients.
#' @param normName the given name of the normalization.
#' @param normCoefs the normalizations coeficients.
#'
#' @return  a rMSI image containing the new normalization field.
#' @export
#'
AppendNormalizationCoefs <- function(img, normName, normCoefs)
{
  if(is.null(img$normalizations))
  {
    img$normalizations <- list()
  }

  img$normalizations[[normName]] <- normCoefs

  return(img)
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
  .plotMassImageRGB(raster_RGB, cal_um2pixels = img$pixel_size_um, rotation = rotate, display_axes = F)
}

#' calMzAxis.
#'
#' @param avgSpc_mz  The mass axis to calibrate.
#' @param ref_mz a vector of reference masses (for exaple the theorical gold peaks).
#' @param target_mz manually slected masses to be fittet to ref_masses (must be the same length than ref_mz).
#' @param use_zoo if true zoo package is used for mass error interpolation, it is much more tolerant to large mass errors. Use it carefuly!.
#'
#' @return the calibrated mass axis.
#' @export
#'
calMzAxis <- function(avgSpc_mz, ref_mz, target_mz, use_zoo = F )
{
  if(use_zoo)
  {
    real_idx <- unlist( lapply(target_mz, function(x){ which.min(abs(x - avgSpc_mz)) }))
    new_mz <- rep(NA, length(avgSpc_mz))
    new_mz[real_idx] <- ref_mz
    error <- new_mz - avgSpc_mz
    error <- zoo::na.spline( error )
    return(avgSpc_mz + error)
  }
  else
  {
    MQrefPeaks <- MALDIquant::createMassPeaks(ref_mz, rep(1, length(ref_mz)))
    MQtargetPeaks <- MALDIquant::createMassPeaks(target_mz, rep(1, length(target_mz)))
    MQwarp <- MALDIquant::determineWarpingFunctions(list(MQtargetPeaks), reference = MQrefPeaks, tolerance = 0.005, method = "lowess")
    MQSpectra <- MALDIquant::createMassSpectrum(mass = avgSpc_mz, intensity = rep(0, length(avgSpc_mz)))
    MQcalibrated <- MALDIquant::warpMassSpectra(list(MQSpectra), MQwarp)
    return(MQcalibrated[[1]]@mass)
  }
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
    cube <- loadImgCunckFromCube(img, ic)

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
      saveImgCunckAtIds( img,  partial_ids, dm[1:length(partial_ids), ] )
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
