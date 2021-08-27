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


#' CalibrationWindow.
#'
#' @param mass The mass vector of spectrum to calibrate.
#' @param intensity The intensity vector of spectrum to calibrate.
#' @param CalibrationSpan the span of the loess method for calibration.
#'
#' @return a list containing the calibration model and a data frame with calibration data. NULL is returned if calibration is aborted.
#'
CalibrationWindow<-function( mass, intensity, peak_win_size = 20, win_title = "", CalibrationSpan = 0.75)
{

  options(guiToolkit="RGtk2") # Forca que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)
  
  ## Get the environment for this
  ## instance of the function.
  this <- environment()
  
  ##Class data members
  spectraWidget <- NULL
  dMass <- mass
  dIntensity <- intensity
  refMz <- NULL #a vector of user selected reference masses
  targetMz <- NULL #a vector of user selected target masses to calibrate
  refNames <- NULL #a vector of user selected reference names to calibrate
  PeakWindow <- peak_win_size
  rm(mass)
  rm(intensity)
  Span <- CalibrationSpan

  Tbl_ColNames <- list(name = "Name", ref = "Ref. m/z", sel = "Sel. m/z", errRAW = "Error RAW [m/z]", ppmRAW = "Error RAW [ppm]", errCAL = "Error Cal. [m/z]", ppmCAL = "Error Cal. [ppm]", active = "Active")
  CalModel <- NULL #Place holder for the calibration model
  TableCalData <- NULL

  #Create an mzTable
  CreateMzTable <- function( ref_names, ref_mz )
  {
    sel_mz <- rep(NaN, length(ref_mz))
    mzTable <- data.frame( ref_names, ref_mz, sel_mz, (ref_mz - sel_mz), (1e6*(ref_mz - sel_mz)/ref_mz), rep(NaN, length(ref_mz)), rep(NaN, length(ref_mz)) , rep(F, length(ref_mz)))
    colnames(mzTable) <- as.vector(unlist(this$Tbl_ColNames))
    return(mzTable)
  }
  
  #Enable calibration button only if are enough mz points to calibrate
  SetCalibrateButtonActiveState <- function()
  {
    if(length(which( this$Table_Ctl[ , this$Tbl_ColNames$active] ) ) >= 2)
    {
      gWidgets2::enabled(Btn_Calibrate) <-T
    }
    else
    {
      gWidgets2::enabled(Btn_Calibrate) <-F
    }
  }
  
  #A click on mass spectra widgets drives here. 
  SpectrumClicked <- function( channel, mass, tol )
  {
    tolppm <- 1e6*tol/mass
    tolppm <- max(c(500, tolppm))
    tol <- tolppm*mass/1e6
  
    idLo <- which.min( abs( (mass - tol) - this$dMass  ) )
    idHi <- which.min( abs( (mass + tol) - this$dMass  ) )
    subMass <- this$dMass[idLo:idHi]
    subIntensity <- this$dIntensity[idLo:idHi]
    pks <- DetectPeaks(subMass, subIntensity, SNR = 0.1, WinSize = this$PeakWindow)
    if( length(pks$intensity) == 0 )
    {
      return()
    }
    
    nearestPeak <- which.min(abs(mass - pks$mass))
    PeakMz <- pks$mass[nearestPeak]
    
    iRow <- this$Table_Ctl$get_selected()
    this$Table_Ctl[iRow, this$Tbl_ColNames$sel] <- PeakMz
    this$Table_Ctl[iRow, this$Tbl_ColNames$errRAW] <- this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz
    this$Table_Ctl[iRow, this$Tbl_ColNames$ppmRAW] <- 1e6*(this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz)/this$Table_Ctl[iRow, this$Tbl_ColNames$ref]
    this$Table_Ctl[iRow, this$Tbl_ColNames$active] <- T
    gWidgets2::svalue(this$Btn_Active) <-  this$Table_Ctl[iRow, this$Tbl_ColNames$active]
    gWidgets2::enabled(this$Btn_Active) <-T
    this$SetCalibrateButtonActiveState()
    
    this$spectraWidget$SetSelectedMassTol(1, PeakMz, 0)
  }
  
  #Load an ASCII with reference masses
  LoadRefMzAscii <- function (...)
  {
    fname <- gWidgets2::gfile("Select reference m/z file", type = "open", multi = F)
    if( length(fname) > 0)
    {
      dataM <- read.table(fname)
      ref_names <- as.vector(dataM[,1])
      ref_mz <- as.vector(dataM[,2])
      this$Table_Ctl$set_items( CreateMzTable( ref_names, ref_mz ))
      gWidgets2::size(this$Table_Ctl) <- list( width = -1, height = -1,  column.widths = rep(110, length(this$Tbl_ColNames)))
      this$spectraWidget$SetRefMass(ref_mz)
      gWidgets2::enabled(this$Table_Ctl) <- T
      gWidgets2::enabled(this$Btn_Active) <- T
      gWidgets2::enabled(this$Btn_AutoAssign) <- T
      this$SetCalibrateButtonActiveState()
      this$spectraWidget$SetActiveTool("Red")
    }
  }
  
  #Asign all ref. masses automatically
  AutoMzAssign <- function(...)
  {
    pks <- DetectPeaks( this$dMass,  this$dIntensity, SNR = 0.5, WinSize = this$PeakWindow) 
    if( length(pks$intensity) == 0 )
    {
      return()
    }

    for( iRow in 1:nrow( this$Table_Ctl))
    {
      nearestPeak <- which.min(abs(this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - pks$mass))
      PeakMz <- pks$mass[nearestPeak]
      this$Table_Ctl[iRow, this$Tbl_ColNames$sel] <- PeakMz
      this$Table_Ctl[iRow, this$Tbl_ColNames$errRAW] <- this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz
      this$Table_Ctl[iRow, this$Tbl_ColNames$ppmRAW] <- 1e6*(this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz)/this$Table_Ctl[iRow, this$Tbl_ColNames$ref]
      this$Table_Ctl[iRow, this$Tbl_ColNames$active] <- T
    }
    
    this$SetCalibrateButtonActiveState()
  }
  
  #Mz Table row selected
  RowSelected <- function ( ... )
  {
    #Zoom in spectrum
    iRow <- this$Table_Ctl$get_selected()
    SelMz <- this$Table_Ctl[ iRow, this$Tbl_ColNames$ref]
    if( SelMz > max(this$dMass) || SelMz < min(this$dMass))
    {
      cat(paste("Selected mass (", SelMz, ") is out of range: ", min(this$dMass), " - " , max(this$dMass),"\n", sep = ""))
      return()
    }
    this$spectraWidget$ZoomMzRange( SelMz - 0.007*SelMz, SelMz + 0.007*SelMz )
    RGtk2::gtkButtonSetLabel( gWidgets2::getToolkitWidget(this$Btn_Active) , paste("m/z",sprintf("%.2f",SelMz), "active") )
    gWidgets2::svalue(this$Btn_Active) <-  this$Table_Ctl[iRow, this$Tbl_ColNames$active]
    
    if(is.nan( this$Table_Ctl[ iRow, this$Tbl_ColNames$sel] ))
    {
      gWidgets2::enabled(Btn_Active) <-F
    }
    else
    {
      gWidgets2::enabled(Btn_Active) <-T
      this$spectraWidget$SetSelectedMassTol(1, this$Table_Ctl[ iRow, this$Tbl_ColNames$sel], 0)
    }
  }
  
  #Button to enable/disable mz
  BtnActiveChanged <- function( ... )
  {
    iRow <- this$Table_Ctl$get_selected()
    if( length(iRow) > 0)
    {
      this$Table_Ctl[iRow, this$Tbl_ColNames$active] <- gWidgets2::svalue(this$Btn_Active)
      this$SetCalibrateButtonActiveState()
    }
  }
  
  #Calibration Button clicked
  BtnCalibrate <- function( ... )
  {
    valid_rows <- which( this$Table_Ctl[, this$Tbl_ColNames$active] )
    
    this$CalModel <- calcMassModel(this$Table_Ctl[valid_rows, this$Tbl_ColNames$ref],
                                   this$Table_Ctl[valid_rows, this$Tbl_ColNames$sel],
                                   tolower(gWidgets2::svalue(this$Rad_Method)),
                                   CalSpan = gWidgets2::svalue(this$Spin_Span) )
    
    dMassCalibrated <- applyMassCalibration(this$CalModel, this$dMass)
    
    
    #Apply the calibration model to all the reference peaks to get an after-calibration error prediction
    tablePeaksCalMass <- applyMassCalibration(this$CalModel, this$Table_Ctl[, this$Tbl_ColNames$sel])
    this$Table_Ctl[, this$Tbl_ColNames$errCAL] <- this$Table_Ctl[, this$Tbl_ColNames$ref] - tablePeaksCalMass
    this$Table_Ctl[, this$Tbl_ColNames$ppmCAL] <- 1e6*(this$Table_Ctl[, this$Tbl_ColNames$errCAL])/this$Table_Ctl[, this$Tbl_ColNames$ref]
    
    #Copy the calibration data table to allow returning it
    this$TableCalData <- this$Table_Ctl[,]
    
    gWidgets2::enabled(Chk_ShowCal) <-T
    
    this$spectraWidget$AddSpectra(  dMassCalibrated, dIntensity, col = "darkgreen", name = "cal")
    gWidgets2::svalue(this$Chk_ShowRaw) <- F
    gWidgets2::svalue(this$Chk_ShowCal) <- T
    gWidgets2::enabled(Btn_Confirm) <- T
    this$spectraWidget$ZoomResetClicked()
    this$spectraWidget$SetActiveTool("Zoom")
  }
  
  #Clicked on checkbox
  ChkShowRawSpc <- function( ... )
  {
    this$spectraWidget$SetSpectrumEnabled("raw", gWidgets2::svalue(this$Chk_ShowRaw))
  }
  
  #Clicked on checkbox
  ChkShowCalSpc <- function( ... )
  {
    this$spectraWidget$SetSpectrumEnabled("cal", gWidgets2::svalue(this$Chk_ShowCal))
  }
  
  #Radio button to select method changed
  MethodRadioChanged <- function( ... )
  {
    gWidgets2::enabled(this$Grp_CalSpan) <- (gWidgets2::svalue(this$Rad_Method) == "Loess")
  }
  
  #Validate calibration and quit
  ValidateCalAndQuit <- function( ... )
  {
    gWidgets2::dispose(this$this$window)
  }
  
  #GUI builder
  this$window <- gWidgets2::gwindow ( paste("Spectrum Calibration -",win_title ), visible = F )
    Grp_Top <- gWidgets2::gpanedgroup(horizontal = F, container = window)
  
  Grp_Bot <- gWidgets2::ggroup(horizontal = T, container = Grp_Top)
  Table_Ctl <- gWidgets2::gtable(CreateMzTable( c(NaN), c(NaN) ), multiple = F, container = Grp_Bot, chosen.col = 2) #First empty MZ table
  
  Table_Ctl$set_editable(F)
  gWidgets2::size(this$Table_Ctl) <- list( width = -1, height = -1,  column.widths = rep(110, length(Tbl_ColNames)))
  gWidgets2::enabled(Table_Ctl) <-F
  
  Grp_Btn <- gWidgets2::ggroup(horizontal = F, container = Grp_Bot)
  Btn_LoadRef <- gWidgets2::gbutton("Load Ref m/z", container = Grp_Btn, handler = this$LoadRefMzAscii)
  Btn_AutoAssign <- gWidgets2::gbutton("Auto assign", container = Grp_Btn, handler = this$AutoMzAssign)
  Btn_Active <- gWidgets2::gcheckbox("m/z active", checked = F, use.togglebutton = T, handler = this$BtnActiveChanged, container = Grp_Btn)
  
  Frm_Method <- gWidgets2::gframe("Calibration Method", container = Grp_Btn)
  Grp_CalMethod <- gWidgets2::ggroup(horizontal = F, container = Frm_Method)
  Rad_Method <- gWidgets2::gradio( items = c("Linear", "Loess"), selected = 2, horizontal = F, handler = this$MethodRadioChanged, container = Grp_CalMethod)
  Grp_CalSpan <- gWidgets2::ggroup(horizontal = T, container = Grp_CalMethod)
  lblSpan <- gWidgets2::glabel("Span:", container = Grp_CalSpan)
  Spin_Span <- gWidgets2::gspinbutton(from = 0.1, to = 2, by = 0.05, value = Span, digits = 2, container = Grp_CalSpan)
  
  Btn_Calibrate <- gWidgets2::gbutton("Calibrate", container = Grp_Btn, handler = this$BtnCalibrate)
  gWidgets2::enabled(Btn_Active) <-F
  gWidgets2::enabled(Btn_AutoAssign) <-F
  gWidgets2::enabled(Btn_Calibrate) <-F
  
  Frm_PlotCtl <- gWidgets2::gframe("Plot control", container = Grp_Btn)
  Grp_PlotCtl <- gWidgets2::ggroup(horizontal = F, container = Frm_PlotCtl)
  Chk_ShowRaw <- gWidgets2::gcheckbox("RAW", checked = T, handler = this$ChkShowRawSpc, container = Grp_PlotCtl)
  Chk_ShowCal <- gWidgets2::gcheckbox("CAL", checked = F, handler = this$ChkShowCalSpc, container = Grp_PlotCtl)
  .setCheckBoxText(Chk_ShowRaw, "RAW spectrum", background = NULL, foreground = "darkblue", font_size = NULL, font_weight = "heavy")
  .setCheckBoxText(Chk_ShowCal, "CAL spectrum", background = NULL, foreground = "darkgreen", font_size = NULL, font_weight = "heavy")
  
  Btn_Confirm <- gWidgets2::gbutton("Validate CAL & quit", container = Grp_Btn, handler = this$ValidateCalAndQuit)
  gWidgets2::enabled(Btn_Confirm) <-F
  
  gWidgets2::enabled(Chk_ShowCal) <-F
  
  spectraFrame<-gWidgets2::gframe("", container = Grp_Top,  fill = T, expand = T, spacing = 5 )
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, top_window_widget = window, clicFuntion = this$SpectrumClicked, showOpenFileButton = F,  display_sel_red = T, display_sel_spins = F)
  
  gWidgets2::size( window )<- c(1024, 740)
  gWidgets2::size( spectraWidget )<- c(-1, 380)
  gWidgets2::size( Table_Ctl )<- c(-1, 200)
  
  visible(window)<-TRUE
  
  spectraWidget$AddSpectra(  dMass, dIntensity, col = "blue", name = "raw")
  spectraWidget$ZoomResetClicked()
  window$widget$present()
  Table_Ctl$add_handler_selection_changed(this$RowSelected)
  
  ## Set the name for the class
  class(this) <- append(class(this),"CalWindow")
  gc()
  
  #Do not return until this window is disposed...
  while(gWidgets2::isExtant(this$window ))
  {
    Sys.sleep(0.1)
  }
  
  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(list( model = this$CalModel,
               data = this$TableCalData
               ) )

}

#' calcMassModel.
#' 
#' Calculate the mass calibration model from a list of target masses and a list of reference masses. 
#' Two methods are available: linear and loess. The linear model is more robust to calibration artifacts and do not exhibit oscillations.
#' But, the loess model can provide a higher mass accuracy when the span parameter is tuned properly. When using the loess model, try to settle with a span value as low as possible.
#' In general, it is advised to use linear model for dataset acquired with high mass accuracy instruments (FTICR, Orbitrap) and keep the loess model for lower mass accuracy instruments (TOF).
#'
#' @param ref_mz a vector of reference masses (for exaple the theorical gold peaks).
#' @param target_mz manually slected masses to be fittet to ref_masses (must be the same length than ref_mz).
#' @param method a string with the method used for interpolation, valid methods are loess and linear
#' @param CalSpan the span for the loess method (ignored in linear mode).
#'
#' @return the calibration model.
#'
#' @export
#'
#' @seealso \code{\link{applyMassCalibration}} for a complete mass calibration example.
#' 
calcMassModel <- function(ref_mz, target_mz, method = "loess", CalSpan = 0.75 )
{
  a <- data.frame( targetMass = target_mz, refMass =  ref_mz)
  if(method == "linear")
  {
    fitmodel <- lm(refMass~targetMass, a)
  }
  else if (method == "loess")
  {
    fitmodel <- loess(refMass~targetMass, a, control = loess.control(surface = "direct"), span = CalSpan)
  }
  else
  {
    stop(paste("The specified method:", method, "is not a valid method\n"))
  }
  
  massModel <- list(method = method, model = fitmodel)
  class(massModel) <- "MassModel"
  
  return(massModel)
}


#' applyMassCalibration.
#' 
#' Calibrates a mass axis using an already calculated mass model. 
#'
#' @param massModel a massModel calculated using the calcMassModel() method.
#' @param massRAW a numeric vector containing the raw mass axist to be calibrated. 
#'
#' @return a numeric vector with the calibrated mass axis.
#' 
#' @export
#'
#' @examples
#' # Set  five target mass peaks as example
#' target_peaks <- c(110, 200, 500, 700, 750) 
#' 
#' # Create a synthetic spectrum containing the five example peaks
#' raw_mass <- seq(from = 100, to = 800, by = 0.05)
#' raw_intensity <- 0.4*dnorm(raw_mass, mean = target_peaks[1], sd = 0.1) +
#'   0.5*dnorm(raw_mass, mean = target_peaks[2], sd = 0.1) +
#'   0.8*dnorm(raw_mass, mean = target_peaks[3], sd = 0.1) +
#'   0.5*dnorm(raw_mass, mean = target_peaks[4], sd = 0.1) +
#'   0.3*dnorm(raw_mass, mean = target_peaks[5], sd = 0.1)
#' 
#' # Set the reference masses for the  five target peaks (theoretical mass values)
#' ref_peaks <- c(110.1, 200.05, 499.95, 700.1, 750.15)
#' 
#' # Calculate the calibration model using the linear and loess methods
#' linModel <- rMSI::calcMassModel(ref_mz = ref_peaks, target_mz = target_peaks, method = "linear")
#' loessModel <- rMSI::calcMassModel(ref_mz = ref_peaks, target_mz = target_peaks, method = "loess", CalSpan = 20)
#' 
#' # Apply the calibration models to get a calibrated spectrum
#' mass_linCal <- rMSI::applyMassCalibration(linModel, raw_mass)
#' mass_loessCal <- rMSI::applyMassCalibration(loessModel, raw_mass)
#' 
#' # Plot the synthetic raw spectrum overlaied with its calibrations
#' rMSI::plotSpectra( raw_mass, raw_intensity, col = "black", ref_mass = ref_peaks) #raw with mass references marked
#' rMSI::plotSpectra( mass_linCal, raw_intensity, col = "red") #linear calibration
#' rMSI::plotSpectra( mass_loessCal, raw_intensity, col = "blue") #loess calibration
#' 
applyMassCalibration <- function(massModel, massRAW)
{
  if(class(massModel) != "MassModel")
  {
    stop("Invalid calibration model. Use the calcMassModel() function to obtain a valid model.\n")
  }
  
  return( predict(massModel$model, newdata = data.frame( targetMass = massRAW)) ) 
}


#' applyMassCalibrationImage.
#' 
#' Apply a mass calibration model to a complete image. The new mass axis will be calculated using the supplied calibration model (massModel).
#' Then, the original imzML file will be overwritten with the calibrated mass axis. In case of data in continuous mode the calibration process is very fast since a 
#' single mass axis is recalculated and overwritten. For imzML in processed mode the calibration process may take some time since a new mass axis must be recalculated
#' for each pixel in the data.
#'
#' @param img A image in rMSI data format .
#' @param massModel a massModel calculated using the calcMassModel() method.
#' 
#' @return a rMSI img object with the mass calibration applied.
#' 
#' @export
#'
#' @seealso \code{\link{applyMassCalibration}} and \code{\link{calcMassModel}} 
#' 
applyMassCalibrationImage<-function(img, massModel)
{
  #Check massModel
  if(class(massModel) !=  "MassModel")
  {
    stop("Invalid calibration model. Use the calcMassModel() function to obtain a valid model.\n")
  }
  
  ibdFile <- path.expand(file.path(img$data$path, paste0(img$data$imzML$file, ".ibd")) )
  
  #Calculate the calibrated common mass axis
  img$mass <- applyMassCalibration(massModel, img$mass)
  
  if( img$data$imzML$continuous_mode )
  {
    #Data in continuous mode, so only the common mass axis must be overwritten
    CimzMLBinWriteModifyMass(ibdFname = ibdFile, 
                             NPixels = nrow(img$pos), 
                             mz_dataTypeString =  img$data$imzML$mz_dataType, 
                             int_dataTypeString = img$data$imzML$int_dataType, 
                             continuous = img$data$imzML$continuous_mode, 
                             mzNew = img$mass,
                             mzOffset = img$data$imzML$run$mzOffset[1]
                             )
  }
  else
  {
    #Data in contiuous mode, so each pixel mass axis must be calibrated
    cat("Calibrating the imzML file...\n")
    pb <- txtProgressBar(min = 0, max = nrow(img$data$imzML$run), initial = 0, style = 3)
    for(i in 1:nrow(img$data$imzML$run))
    {
      setTxtProgressBar(pb, i)
      #Read imzML mass axis to the mass variable
      mass <- CimzMLBinReadMass(ibdFname = ibdFile, 
                                NPixels = nrow(img$pos), 
                                N = img$data$imzML$run$mzLength[i], 
                                offset = img$data$imzML$run$mzOffset[i], 
                                dataTypeString = img$data$imzML$mz_dataType, 
                                continuous = img$data$imzML$continuous_mode
                                )
      
      #Calibrate
      mass <- applyMassCalibration(massModel, mass)
      
      #Overwrite the imzML mass axis
      CimzMLBinWriteModifyMass(ibdFname = ibdFile, 
                               NPixels = nrow(img$pos), 
                               mz_dataTypeString =  img$data$imzML$mz_dataType, 
                               int_dataTypeString = img$data$imzML$int_dataType, 
                               continuous = img$data$imzML$continuous_mode, 
                               mzNew = mass,
                               mzOffset = img$data$imzML$run$mzOffset[i]
                               )
    }
    close(pb)
  }
  
  return(img)
}
