#MSI Methods


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
