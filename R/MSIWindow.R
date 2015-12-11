#MALDI Image reconstruction by selected ion Top Windows
##########################################################

#Open img from Hdd directly, the easy way...
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
  raw<-LoadMsiData(data_file = fname, restore_path = file.path(dirname(fname), paste("ramdisk",basename(fname), sep = "_")), fun_progress_event = mPBar$setValue)

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

MSIWindow <- function( img )
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
  AddSpectra <- function (...)
  {
    this$spectraWidget$AddSpectra(...)
  }

  #GUI builder
  window <- gwindow ( "MSI Reconstruction" , visible = F )
  Grp_Top <- gpanedgroup(horizontal = F, container = window)
  msiWidget <- .MSImagePlotWidget(in_img = img , parent_widget = Grp_Top, AddSpectra_function = this$AddSpectra)
  spectraFrame<-gframe("Average Spectra", container = Grp_Top,  fill = T, expand = T, spacing = 5 )
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, clicFuntion = this$SpectrumClicked, showOpenFileButton = F)

  visible(window)<-TRUE

  ##TODO revisar aixo quan tinguis multiples plots com ho faras???
  spectraWidget$AddSpectra(  img$mass, img$mean@intensity, col = "red")
  spectraWidget$ZoomResetClicked()


  window$widget$present() ##Tot i forzzar aki un raise-up del top widget en RStudio no deixa!

  ## Set the name for the class
  class(this) <- append(class(this),"MsiWindows")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)
}
