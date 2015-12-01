#MSI Methods


#Load the ff image displaying a progress bar
#img_full_path - Full path to image .tar file
#Returns - The image object
LoadImageWithProgressBar <- function( img_full_path )
{
  setPbarValue<-function(progress)
  {
    setTxtProgressBar(pb, progress)
  }

  cat("Loading Image...\n")
  pt<-proc.time()
  pb<-txtProgressBar(min = 0, max = 100, style = 3 )
  RawData <- LoadMsiData(data_file = img_full_path, restore_path = file.path(dirname(img_full_path), paste("ramdisk",basename(img_full_path), sep = "_")), fun_progress_event = setPbarValue)
  close(pb)
  pt<-proc.time() - pt
  cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))

  return(RawData)
}

#Remove The ramdisk
RemoveImageRamdisk<-function(img)
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


#Camputes the warping function to align the average spectra to a given reference masses
#img - A image in ff data format
#ref - vector of reference masses
#snr - Minimum signal to noise ratio to keep a peak
#Tolerance - double, maximal deviation of a peak position (mass) to be considered as identical
#Method - Method used to determine the warping function. Options are: "lowess", "linear", "quadratic", "cubic"
#PlotResult - If true, opens a spectrawidget windows to check if the warping function is working properly
#             Original average spectrum is ploted in blue the calibrated spetrum is ploted in red
#Returns - The new mass axis
ObtainCalibrationFunction<-function(img, ref, snr = 10, Tolerance = 0.01, Method = "lowess", PlotResult = T)
{
  mean_pks<-detectPeaks(img$mean, SNR=snr, method= "SuperSmoother")
  Peak_refMass<-createMassPeaks(ref, rep(1, length(ref)))
  warp<-determineWarpingFunctions(mean_pks, Peak_refMass,  tolerance = Tolerance, method= Method)
  mean_cal <- warpMassSpectra(list(img$mean), warp)[[1]]

  if(PlotResult)
  {
    source('~/R/SCRIPTS/spectraWidget.R')
    plotSpc<-plotSpectra()
    plotSpc$SetRefMass(ref)
    plotSpc$AddSpectra(img$mean@mass, img$mean@intensity, col = "blue")
    plotSpc$AddSpectra(mean_cal@mass, mean_cal@intensity, col = "red")
    plotSpc$.ZoomResetClicked()
  }
  return(mean_cal@mass)
}


#Apply a new mass axis to a complete image
#img - A image in ff data format
#new_mass - The new mass axis
#output_fname - full path to store the output image
#Return - Nothing! New image is stored in hdd at spesified path
CalibrateImage<-function(img, new_mass, output_fname)
{
  #Copy the img objet
  calImg <- img

  #Replace mass axis
  calImg$mass <- new_mass
  calImg$mean@mass <- new_mass

  #Store
  SaveMsiData(output_fname, calImg, calImg$mean)
}
