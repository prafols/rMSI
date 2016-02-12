###TODO work here... input ha de ser un datacube, un sol objecte complet per facilitar guardar iumatges de PCA i altres cuentos....
###TODO els metodes d'aki estic segur k poden anar a la api d'acces a objecte rMSI....

#Save image using custom format
#data_file - full path to hdd location to save image as a tar file
#imgData - the image to save in custom ff data format
#meanSpcData - Maldiquant object containing the average spectrum

#' Save a rMSI object to disk in a compressed .tar file.
#'
#' @param data_file Full path to hdd location to save image as a tar file.
#' @param imgData TODO.
#' @param meanSpcData TODO.
#' @param um2pixel TODO.
#'
#' @return
#' @export
#'
#' @examples
SaveMsiData<-function(data_file, imgData, meanSpcData, um2pixel)
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
  save(meanSpcData, file = file.path(data_dir, "mean.SpcR")) #Save mean spectra
  resolutionObj<-um2pixel
  save(resolutionObj, file = file.path(data_dir, "pixel_size_um.ImgR")) #Save pixel size um Object


  #Store also a vector of names of ff data in order to be able to restore it
  ffDataNames <- paste(names(imgData$data), "_ffzip", sep = "")
  save(ffDataNames, file = file.path(data_dir, "ffnames.ImgR")) #Save ff filenames Object

  #Save the ff object
  pb<-txtProgressBar(min = 0, max = length(ffDataNames), style = 3 )
  for(i in 1:length(ffDataNames))
  {
    setTxtProgressBar(pb, i)
    ffObj<-imgData$data[[i]] #The ff::ffsave can only handle ff objects directly, not in lists!
    ff::ffsave(ffObj , file =  file.path(data_dir, ffDataNames[i]))
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
#' @param fun_progress_event This is a callback function to update the progress of loading data. See details for more information.
#' @param ff_overwrite Tell ff to overwrite or not current ramdisk files.
#'
#' @return an rMSI object pointing to ramdisk stored data
#'
#' Loads a rMSI data object from .tar compressed file. It will be uncompressed at specified restore_path.
#' fun_progress_event can be NULL or a function with the following prototipe: fun_progress_event( currentState, TotalNumberOfStates).
#' This function will be called periodically to monitor the loading status. This is usefull to implement progressbars.
#' If ramdisk is already created befor calling this method the parameter ff_overwrite will control the loadin behaviour. If it is set to false (default)
#' The ramdisk will be kept and the imaged loaded imediatelly. Otherwise if is set to true, the while dataset will be reloaded from tar file.
#'
LoadMsiData<-function(data_file, restore_path, fun_progress_event = NULL, ff_overwrite = F)
{
  #1- Check if the specified image ramdisk exists in the restore_path location
  datacube<-.FastLoad(restore_path)
  if(!is.null(datacube) && !ff_overwrite)
  {
    if(!is.null(fun_progress_event))
    {
      fun_progress_event(100)
    }
    return(datacube)
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
    print("Warining: Old image without resolution object. It is set to 9999 um by default!")
    resolutionObj <- 9999
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
    if(!is.null(fun_progress_event) && (round(pp) > round(pp_ant)) )
    {
      #Update progress bar
      if( !fun_progress_event(pp) )
      {
        return(NULL) #progress ar function must return true if the loading process is to be continued.
      }
    }
  }

  lapply(spectra, ff::open.ff)
  datacube<-list(name = basename(data_file), mass = massObj,  size = sizeObj,  pos = posObj, pixel_size_um = resolutionObj, mean = meanSpcData, data = spectra)
  unlink("ImgData", recursive = T)

  save(datacube, file = file.path(restore_path, "datacube.RImg"))
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
