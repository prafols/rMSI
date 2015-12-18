#MALDI Imaging pack 2 R data file

#data_processing_function must be a function with a MALDIquant object input and MALDIquant object output
MSI2Rdata<-function(raw_data_full_path, resolution_um, spot_selection_xml_file_full_path, output_data_filename, data_processing_function = NULL, ...)
{
  cat("Importing data to R session...\n")
  pt<-proc.time()
  ff_folder<-file.path(raw_data_full_path, "ffdata")
  dir.create(ff_folder)
  raw<-importBrukerXmassImg(raw_data_folder = raw_data_full_path, xml_file = spot_selection_xml_file_full_path, ff_data_folder = ff_folder)
  pt<-proc.time() - pt
  cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))

  #OLD CODE I will not convert to MALDIquant the whole dataset anymore
  ##cat("Converting data to MALDIquant object...\n")
  ##imgQuant<-ConvertMyImg2MassSpectrum(imgOrig)
  ##gc()

  cat("Calculating average spectrum...\n")
  pt<-proc.time()

  avgI<-apply(matrix(unlist(lapply(raw$data, function(x){ ff::ffrowapply(colSums(x[i1:i2,,drop=FALSE]), X=x, RETURN = TRUE, CFUN = "csum", FF_RETURN = FALSE) })), nrow = length(raw$data), byrow = T), 2, sum)
  avgI<-avgI/( sum(unlist(lapply(raw$data, nrow))) )

  ###TODO test if new splited approach works an remove the old one
  ###avgI<-(ff::ffrowapply(colSums(raw$data[i1:i2,,drop=FALSE]), X=raw$data, RETURN = TRUE, CFUN = "csum", FF_RETURN = FALSE))/nrow(raw$data)

  meanSpc<-createMassSpectrum(mass =  raw$mass, intensity = avgI)
  meanSpc<-smoothIntensity(meanSpc, halfWindowSize=2)
  pt<-proc.time() - pt
  cat(paste("Average spectrun time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  gc()

  ##TODO commented old code must be re.implemented with the new ff approach
  #Custom data processing
#   if(!is.null(data_processing_function))
#   {
#     startTime<-proc.time()
#     cat("Data Processing ", length(imgQuant) , "spectras...\n")
#     imgQuant<-data_processing_function(imgQuant, ...)
#     elapsedTime<-proc.time() - startTime
#     cat("Processing done in ", elapsedTime["elapsed"], " seconds\n" )
#
#     if(output_data_filename != "")
#     {
#       cat("Converting to imgOrig\n")
#       imgOrig<-ConvertMassSpectrum2MyImg(imgQuant)
#     }
#     else
#     {
#       return(imgQuant)
#     }
#   }
#   gc()

  if(output_data_filename != "")
  {
    #TODO commented old code
#     cat("Creating JPEG compressed imgage slices...\n")
#     img_jpeg<-CompressMyImg2JpegSlices(imgOrig)
#     save(img_jpeg, file = paste(output_data_filename, "jpeg", sep = "_" ), compress = T)

    cat("Packaging R data objects...\n")
    pt<-proc.time()
    SaveMsiData(output_data_filename, raw, meanSpc, resolution_um)
    pt<-proc.time() - pt
    cat(paste("Data saving time:", round(pt["elapsed"], digits = 1),"seconds\n"))

    cat("Done. Data is saved in:\n",  output_data_filename, "\n")
    gc()
  }

  #Remove ff file created on HDD
  lapply(raw$data, ff::delete)
  rm(raw)
  unlink(ff_folder, recursive = T) #Remove folder containing ffdata
}


#Save image using custom format
#data_file - full path to hdd location to save image as a tar file
#imgData - the image to save in custom ff data format
#meanSpcData - Maldiquant object containing the average spectrum
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

#Load image using custom format
#fun_progress_event is a callback function to update the progress of loading data
#Prototipe must be: fun_progress_event( currentState, TotalNumberOfStates)
LoadMsiData<-function(data_file, restore_path, fun_progress_event = NULL, ff_overwrite = T)
{
  #1- Check if the specified image ramdisk exists in the restore_path location
  datacube<-FastLoad(restore_path)
  if(!is.null(datacube))
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
    ff::ffload(file.path(img_path, ffDataNames[i]), rootpath = restore_path, overwrite = ff_overwrite)
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
FastLoad<-function( restorePath )
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
