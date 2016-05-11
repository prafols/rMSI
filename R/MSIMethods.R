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


#' Apply a new mass axis to a complete image.
#' A GUI to calibrate the mean spectrum will be shown.
#'
#' @param A image in rMSI data format .
#' @param full path to store the output image.
#'
#' @export
#'
CalibrateImage<-function(img, output_fname)
{
  #Copy the img objet
  calImg <- img

  #Replace mass axis using GUI calibration
  if( class( calImg$mean) == "MassSpectrum")
  {
    #Old mean MALDIquant handling
    mIntensity <- calImg$mean@intensity
  }
  else
  {
    mIntensity <- calImg$mean
  }
  new_mass <- CalibrationWindow( calImg$mass, mIntensity )

  if(is.null(new_mass))
  {
    cat("Calibration process aborted by user\n")
    return()
  }

  #Ask for confirmation
  resp <- ""
  while(resp != "y" && resp != "n" && resp != "Y" && resp != "N")
  {
    resp <- readline(prompt = "Proced with calibration of the whole MS image? (this may thake some time) [y, n]:")
    if(resp != "y" && resp != "n" && resp != "Y" && resp != "N")
    {
      cat("Invalid response, valid responses are: y, n, Y and N. Try again.\n")
    }
  }

  if( resp == "n" || resp =="N")
  {
    cat("Calibration process aborted by user\n")
    return()
  }

  calImg$mass <- new_mass

  if( class( calImg$mean) == "MassSpectrum")
  {
    calImg$mean@mass <- new_mass
  }

  #Store
  SaveMsiData(calImg, output_fname)
}
