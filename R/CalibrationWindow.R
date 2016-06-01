#' CalibrationWindow.
#'
#' @param mass The mass vector of spectrum to calibrate.
#' @param intensity The intensity vector of spectrum to calibrate.
#' @param use_zoo if the zoo package interpolation must be used
#'
#' @return a the calibrated mass axis.
#' @export
#'
CalibrationWindow<-function( mass, intensity, use_zoo = F)
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  spectraWidget <- NULL
  dMass <- mass
  dIntensity <- intensity
  rm(mass)
  rm(intensity)

  Tbl_ColNames <- list(name = "Name", ref = "Ref. m/z", sel = "Sel. m/z", err = "Error [m/z]", ppm = "Error [ppm]", active = "Active")
  dMassCalibrated <- NULL

  #Create an mzTable
  CreateMzTable <- function( ref_names, ref_mz )
  {
    sel_mz <- rep(NaN, length(ref_mz))
    mzTable <- data.frame( ref_names, ref_mz, sel_mz, (ref_mz - sel_mz), (1e6*(ref_mz - sel_mz)/ref_mz), rep(F, length(ref_mz)))
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

  #A click on mass spectra widgets drives here. From here image recostruction will be called for various widgets
  SpectrumClicked <- function( channel, mass, tol )
  {
    idLo <- which.min( abs( (mass - tol) - this$dMass  ) )
    idHi <- which.min( abs( (mass + tol) - this$dMass  ) )
    sel_int <- this$dIntensity
    sel_int[ -(idLo:idHi) ] <- 0
    idPk <- which.max(sel_int )
    PeakMz <- dMass[idPk]

    iRow <- this$Table_Ctl$get_selected()
    this$Table_Ctl[iRow, this$Tbl_ColNames$sel] <- PeakMz
    this$Table_Ctl[iRow, this$Tbl_ColNames$err] <- this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz
    this$Table_Ctl[iRow, this$Tbl_ColNames$ppm] <- 1e6*(this$Table_Ctl[iRow, this$Tbl_ColNames$ref] - PeakMz)/this$Table_Ctl[iRow, this$Tbl_ColNames$ref]
    this$Table_Ctl[iRow, this$Tbl_ColNames$active] <- T
    gWidgets2::svalue(this$Btn_Active) <-  this$Table_Ctl[iRow, this$Tbl_ColNames$active]
    gWidgets2::enabled(Btn_Active) <-T
    this$SetCalibrateButtonActiveState()

    return(list(selMz = PeakMz, selTol = 0))
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
      gWidgets2::enabled(Table_Ctl) <- T
      gWidgets2::enabled(Btn_Active) <- T
      this$SetCalibrateButtonActiveState()
      this$spectraWidget$SetActiveTool("Red")
     }
  }

  #Mz Table row selected
  RowSelected <- function ( ... )
  {
    #Zoom in spectrum
    iRow <- this$Table_Ctl$get_selected()
    SelMz <- this$Table_Ctl[ iRow, this$Tbl_ColNames$ref]
    this$spectraWidget$ZoomMzRange( SelMz - 0.005*SelMz, SelMz + 0.005*SelMz )
    RGtk2::gtkButtonSetLabel( gWidgets2::getToolkitWidget(this$Btn_Active) , paste("m/z",sprintf("%.2f",SelMz), "active") )
    gWidgets2::svalue(this$Btn_Active) <-  this$Table_Ctl[iRow, this$Tbl_ColNames$active]

    if(is.nan( this$Table_Ctl[ iRow, this$Tbl_ColNames$sel] ))
    {
      gWidgets2::enabled(Btn_Active) <-F
    }
    else
    {
      gWidgets2::enabled(Btn_Active) <-T
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
    refMz <- this$Table_Ctl[valid_rows, this$Tbl_ColNames$ref]
    targetMz <- this$Table_Ctl[valid_rows, this$Tbl_ColNames$sel]
    this$dMassCalibrated <- calMzAxis(this$dMass, refMz, targetMz, use_zoo)
    gWidgets2::enabled(Chk_ShowCal) <-T

    this$spectraWidget$AddSpectra(  this$dMassCalibrated, dIntensity, col = "darkgreen", name = "cal")
    gWidgets2::svalue(this$Chk_ShowRaw) <- F
    gWidgets2::svalue(this$Chk_ShowCal) <- T
    gWidgets2::enabled(Btn_Confirm) <- T
    this$spectraWidget$ZoomResetClicked()
    this$spectraWidget$SetActiveTool("Zoom")
  }

  #Clicked on checkbox
  ChkShowRawSpc <- function( ... )
  {
    spectraWidget$SetSpectrumEnabled("raw", gWidgets2::svalue(this$Chk_ShowRaw))
  }

  #Clicked on checkbox
  ChkShowCalSpc <- function( ... )
  {
    spectraWidget$SetSpectrumEnabled("cal", gWidgets2::svalue(this$Chk_ShowCal))
  }

  #Validate calibration and quit
  ValidateCalAndQuit <- function( ... )
  {
    gWidgets2::dispose(this$this$window)
  }

  #GUI builder
  this$window <- gWidgets2::gwindow ( "Spectrum Calibration" , visible = F )
  Grp_Top <- gWidgets2::gpanedgroup(horizontal = F, container = window)

  Grp_Bot <- gWidgets2::ggroup(horizontal = T, container = Grp_Top)
  Table_Ctl <- gWidgets2::gtable(CreateMzTable( c(NaN), c(NaN) ), multiple = F, container = Grp_Bot, chosen.col = 2) #First empty MZ table

  Table_Ctl$set_editable(F)
  gWidgets2::size(this$Table_Ctl) <- list( width = -1, height = -1,  column.widths = rep(110, length(Tbl_ColNames)))
  gWidgets2::enabled(Table_Ctl) <-F

  Grp_Btn <- gWidgets2::ggroup(horizontal = F, container = Grp_Bot)
  Btn_LoadRef <- gWidgets2::gbutton("Load Ref m/z", container = Grp_Btn, handler = this$LoadRefMzAscii)
  Btn_Active <- gWidgets2::gcheckbox("m/z active", checked = F, use.togglebutton = T, handler = this$BtnActiveChanged, container = Grp_Btn)
  Btn_Calibrate <- gWidgets2::gbutton("Calibrate", container = Grp_Btn, handler = this$BtnCalibrate)
  gWidgets2::enabled(Btn_Active) <-F
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
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, top_window_widget = window, clicFuntion = this$SpectrumClicked, showOpenFileButton = F,  display_sel_red = T)

  gWidgets2::size( window )<- c(800, 600)
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

  #Return the mass axis... it is NULL if no calibration was done
  return(this$dMassCalibrated)
}
