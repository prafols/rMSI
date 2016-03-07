#MALDI Image reconstruction by selected ion Top Windows
##########################################################

#' Open a MS image from Hdd directly allowing the user to chose it from a file dialog.
#'
#' A file choser dialog is presented to allow the selection of a .tar file containing an rMSI data object.
#' If a ramdisk has been previously created it will be used to speed up the loading process.
#' The image will be presented using the GUI and it will be returned as rMSI object.
#'
#' @return the rMSI object containing the MS image.
#'
#' @export
OpenMSI<-function()
{
  fname<-gfile("Select an MSI file to open", type="open", multi = F, filter =  c("tar"="tar"), initial.dir = path.expand("~/"))
  if(length(fname) == 0)
  {
    return ()
  }

  #Create a progress bar
  mPBar<-.ProgressBarDialog("Loading data please wait...")

  #Preloading 2 speedup
  raw<-LoadMsiData(data_file = fname, restore_path = file.path(dirname(fname), paste("ramdisk",basename(fname), sep = "_")), fun_progress = mPBar$setValue)

  if(is.null(raw))
  {
    #Process Aborted By User, return
    unlink( file.path(dirname(fname), paste("ramdisk",basename(fname), sep = "_")), recursive = T)
    return()
  }

  #Close the progressBar when data is loaded
  mPBar$close()

  #Test for errors...
  if(exists("raw"))
  {
    MSIWindow(img = raw)
    return(raw)
  }
  else
  {
    gmessage("Error: The selected file is not valid.", icon = "error", title = "Load error")
  }
}


#' Open the GUI to explore a MS image
#'
#' @param img a rMSI data object
#'
#'  Open the GUI to explore the MS image provided as parameter.
#'
#' @export
MSIWindow<-function(img)
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  spectraWidget <- 0
  msiWidget <- 0 #TODO currently only one but must be a vector to allow multiple images loaded

  #A click on mass spectra widgets drives here. From here image recostruction will be called for various widgets
  SpectrumClicked <- function( channel, mass, tol )
  {
    ##TODO addapt for various widgets
    ret<-this$msiWidget$ImgBuildFun(channel, mass, tol)
    this$msiWidget$PlotMassImageRGB()
    return(ret)
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  AddSpectra <- function ( mass_axis, intensity_list, color_list)
  {
    ##TODO ill need to rework that for multiple image display
    this$spectraWidget$ClearSpectra()

    for( i in 1:length(intensity_list))
    {
      this$spectraWidget$AddSpectra(mass_axis, intensity_list[[i]], col = color_list[[i]])
    }
  }

  #GUI builder
  window <- gWidgets2::gwindow ( "MSI Reconstruction" , visible = F )
  Grp_Top <- gWidgets2::gpanedgroup(horizontal = F, container = window)
  msiWidget <- .MSImagePlotWidget(in_img = img , parent_widget = Grp_Top, AddSpectra_function = this$AddSpectra)
  spectraFrame<-gWidgets2::gframe("Average Spectra", container = Grp_Top,  fill = T, expand = T, spacing = 5 )
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, top_window_widget = window, clicFuntion = this$SpectrumClicked, showOpenFileButton = F)

  visible(window)<-TRUE

  ##TODO revisar aixo quan tinguis multiples plots com ho faras???
  if( class( img$mean) == "MassSpectrum")
  {
    #Old mean MALDIquant handling
    spectraWidget$AddSpectra(  img$mass, img$mean@intensity, col = "red")
  }
  else
  {
    spectraWidget$AddSpectra(  img$mass, img$mean, col = "red")
  }
  spectraWidget$ZoomResetClicked()

  #Plot a initial image which is the maximum peak in mean spectrum with a tolereance of 100 ppm of the mass range
  if( class( img$mean) == "MassSpectrum")
  {
    SpectrumClicked( channel = 1, mass = img$mass[which.max(img$mean@intensity)], tol = 100/1e6*(max(img$mass)-min(img$mass)))
  }
  else
  {
    SpectrumClicked( channel = 1, mass = img$mass[which.max(img$mean)], tol = 100/1e6*(max(img$mass)-min(img$mass)))
  }


  window$widget$present() ##Tot i forzzar aki un raise-up del top widget en RStudio no deixa!

  ## Set the name for the class
  class(this) <- append(class(this),"MsiWindows")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
}

#Global method to set checkbox text whit colors, font size and weight
.setCheckBoxText <- function(checkbox, text, background = NULL, foreground = NULL, font_size = NULL, font_weight = NULL)
{
  pango_str <- "<span"
  if(!is.null(foreground))
  {
    pango_str <- paste(pango_str, " foreground=\"", foreground,"\"", sep = "" )
  }
  if(!is.null(background))
  {
    pango_str <- paste(pango_str, " background=\"", background,"\"", sep = "" )
  }
  if(!is.null(font_size))
  {
    pango_str <- paste(pango_str, " size=\"",font_size, "\"", sep = "" )
  }
  if(!is.null(font_weight))
  {
    pango_str <- paste(pango_str, " weight=\"", font_weight, "\"", sep = "" )
  }
  pango_str <- paste(pango_str, ">", text, "</span>", sep = "" )

  RGtk2::gtkLabelSetMarkup(gWidgets2::getToolkitWidget(checkbox)[[1]],pango_str)
}
