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
  #Load data using GUI
  MSI_obj<-LoadTwoMsImages()

  #Display the MSIWindow
  if( !is.null(MSI_obj$img1_obj) && is.null(MSI_obj$img2_obj) )
  {
    MSIWindow(img1 = MSI_obj$img1_obj)
  }
  if( is.null(MSI_obj$img1_obj) && !is.null(MSI_obj$img2_obj) )
  {
    MSIWindow(img1 = MSI_obj$img2_obj)
  }
  if( !is.null(MSI_obj$img1_obj) && !is.null(MSI_obj$img2_obj) )
  {
    MSIWindow(img1 = MSI_obj$img1_obj, img2 = MSI_obj$img2_obj)
  }

  return(list(img1 = MSI_obj$img1_obj, img2 = MSI_obj$img2_obj ))
}


#' Open the GUI to explore a MS image
#'
#' @param img a rMSI data object
#'
#'  Open the GUI to explore the MS image. A MS image can be provided as a parameter.
#'  Up to two images can be displayed at once using multiple arguments.
#'
#' @export
MSIWindow<-function(img1, img2 = NULL)
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  spectraWidget <- NULL
  msiWidget1 <- NULL
  msiWidget2 <- NULL

  #A click on mass spectra widgets drives here. From here image recostruction will be called for various widgets
  SpectrumClicked <- function( channel, mass, tol )
  {
    ret1<-this$msiWidget1$ImgBuildFun(channel, mass, tol)
    this$msiWidget1$PlotMassImageRGB()
    if(!is.null(this$msiWidget2))
    {
      ret2<-this$msiWidget2$ImgBuildFun(channel, mass, tol)
      this$msiWidget2$PlotMassImageRGB()
    }
    else
    {
      ret2 <- ret1
    }
    return(list(selMz = mean(c(ret1$selMz, ret2$selMz)), selTol = mean(c(ret1$selTol, ret2$selTol))))
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  AddSpectra <- function ( mass_axis, intensity_list, color_list, id_list, calling_image_string, normalizations,  mz_min, mz_max )
  {
    for( i in 1:length(intensity_list))
    {
      this$spectraWidget$AddSpectra(mass_axis, intensity_list[[i]]/normalizations[i], col = color_list[[i]], name = paste(calling_image_string,"_ID",as.character(id_list[[i]]), sep = ""))
    }

    if( !is.null(mz_min) && !is.null(mz_max))
    {
      this$spectraWidget$SetPlottedMassRange(mz_min, mz_max)
    }
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  ClearSpectraPlot <- function(id_list, calling_image_string)
  {
    for( i in 1:length(id_list))
    {
      this$spectraWidget$RmSpectra(paste(calling_image_string,"_ID",id_list[i], sep = ""))
    }
  }

  #A connector between spectraWidget and msiWidgets because them can not be joined directly
  GetPlotedSpectraInfo <- function( strImgName )
  {
    SpcData <-this$spectraWidget$GetSpectraInfo()
    myImgIds <- c ()
    myImgColors <- c()
    for( i in 1:nrow(SpcData))
    {
      name_id <- unlist(strsplit(as.character(SpcData[i, "names"]), "_ID"))
      if(length(name_id) == 2)
      {
        if( name_id[1] == strImgName )
        {
          myImgIds <- c(myImgIds, as.numeric(name_id[2]) )
          myImgColors <- c(myImgColors, as.character(SpcData[i, "colors"]))
        }
      }
    }

    #Get plotted mass range
    mz_range <- this$spectraWidget$GetPlottedMassRange()

    return(list(ID=myImgIds, color =  myImgColors, mz_min = mz_range$mz_min, mz_max = mz_range$mz_max))
  }

  #GUI builder
  window <- gWidgets2::gwindow ( "MSI Reconstruction" , visible = F )
  Grp_Top <- gWidgets2::gpanedgroup(horizontal = F, container = window)
  Grp_Ims <- gWidgets2::gpanedgroup(horizontal = T, container = Grp_Top)
  msiWidget1 <- .MSImagePlotWidget(in_img = img1 , parent_widget = Grp_Ims, AddSpectra_function = this$AddSpectra, GetSpectraInfo_function = this$GetPlotedSpectraInfo, ClearSpectraPlot_function = this$ClearSpectraPlot, meanSpectrumColor = "red", widget_name = "imgLeft")
  if( !is.null(img2))
  {
    msiWidget2 <- .MSImagePlotWidget(in_img = img2 , parent_widget = Grp_Ims, AddSpectra_function = this$AddSpectra, GetSpectraInfo_function = this$GetPlotedSpectraInfo, ClearSpectraPlot_function = this$ClearSpectraPlot, meanSpectrumColor = "blue", widget_name = "imgRight")
  }
  spectraFrame<-gWidgets2::gframe("Average Spectra", container = Grp_Top,  fill = T, expand = T, spacing = 5 )
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, top_window_widget = window, clicFuntion = this$SpectrumClicked, showOpenFileButton = F,  display_sel_red = T, display_sel_green = T, display_sel_blue = T)

  visible(window)<-TRUE

  if( class( img1$mean) == "MassSpectrum")
  {
    #Old mean MALDIquant handling
    spectraWidget$AddSpectra(  img1$mass, img1$mean@intensity, col = "red", name = "imgLeft_ID0")
  }
  else
  {
    spectraWidget$AddSpectra(  img1$mass, img1$mean, col = "red", name = "imgLeft_ID0")
  }

  if(!is.null(img2))
  {
    if( class( img2$mean) == "MassSpectrum")
    {
      #Old mean MALDIquant handling
      spectraWidget$AddSpectra(  img2$mass, img2$mean@intensity, col = "blue", name = "imgRight_ID0")
    }
    else
    {
      spectraWidget$AddSpectra(  img2$mass, img2$mean, col = "blue", name = "imgRight_ID0")
    }
  }

  spectraWidget$ZoomResetClicked()

  #Plot a initial image which is the maximum peak in mean spectrum with a tolereance of 100 ppm of the mass range
  if( class( img1$mean) == "MassSpectrum")
  {
    SpectrumClicked( channel = 1, mass = img1$mass[which.max(img1$mean@intensity)], tol = 100/1e6*(max(img1$mass)-min(img1$mass)))
  }
  else
  {
    SpectrumClicked( channel = 1, mass = img1$mass[which.max(img1$mean)], tol = 100/1e6*(max(img1$mass)-min(img1$mass)))
  }

  window$widget$present()
  ###RGtk2::gtkWindowMaximize(gWidgets2::getToolkitWidget(window)) #Start maximized, currently disabled to addres windows redraw issue

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
