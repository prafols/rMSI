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

###A GUI to display spectra in an interactive way

#' Plot Mass Spectra in a interactive way.
#'
#' @param mass the mass axis of spectrum
#' @param intensity a intensity vector of spectrum
#' @param peaks_mass a vector of peak masses to be labeled on spectrum
#' @param peaks_intensity a vector of peak intensities to be labeled on spectrum
#' @param ref_mass a vector of reference masses to be represented as vertical dashed lines
#' @param col the color for the spectrum
#'
#' Simply call this function each time you want to plot a mass spectrum and it will be overlayed with the current spectra.
#' Peaks can be specified in the same call and it will be labeled over the spectrum.
#'
#' @export
#'
plotSpectra<-function( mass = NULL, intensity = NULL, peaks_mass = NULL, peaks_intensity = NULL, ref_mass = NULL, col = "" )
{
  if( !exists( x = ".SpectraWidget", mode = "environment") )
  {
    ###TODO: Add a left-list showing all loaded spectra. The list will also allow to disable one of them, deleting, coloring...
    ###TODO: Yes the list will be added here because MSIWindows has its own list. So I don't want to modify spectraWidget
    #Spectra plot windows does not exists, create it
    window_spectra <- gWidgets2::gwindow ( "Spectra Plot" , visible = F , width = 750, height = 440)
    .SpectraWidget<<-.SpectraPlotWidget( parent_widget = window_spectra )
    gWidgets2::addHandlerDestroy( obj = window_spectra, handler = .plotSpectraWindow_Disposed )

    #window_spectra$widget$present()
    gWidgets2::visible(window_spectra) <- T

  }

  if( !is.null(mass) && !is.null(intensity) )
  {
    .SpectraWidget$AddSpectra( mass_data = mass, intensity_data = intensity,
                             mass_peaks = peaks_mass, intensity_peaks = peaks_intensity, col = col )
  }

  if( !is.null(ref_mass) )
  {
    .SpectraWidget$SetRefMass( ref_mass  )
  }
}

.plotSpectraWindow_Disposed <- function (evt, ...)
{
  rm(.SpectraWidget, envir = .GlobalEnv)
  gc()
}

.SpectraPlotWidget <- function( parent_widget=gwindow ( "Default SpectraPlotWidget" , visible = FALSE ), top_window_widget = NULL,  clicFuntion = NULL, showOpenFileButton = T,
                                display_sel_red = F, display_sel_green = F, display_sel_blue = F, max_spectra_limit = 50)
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members

  parent <- parent_widget
  rm(parent_widget)

  if(is.null(top_window_widget))
  {
    top_window <- parent
  }
  else
  {
    top_window <- top_window_widget
  }
  rm(top_window_widget)
  clicFun <- clicFuntion
  rm(clicFuntion)
  AnteriorClickMz <- 0
  plot_device <- 0
  radio_buttons <- 0
  lbl_mz_coords <- 0
  lbl_in_coords <- 0
  LABEL_LENGTH <- 10
  MouseWheelFunction <- 0
  spectra_data <- list() #Data structure controlled by functions
  ref_mass <- NULL
  mz_lim <- c(0, 1)
  in_lim <- c(0, 1)
  data_mass_range <- c() #c(min(mass_data),max(mass_data)),
  PointerCoords <- c(0,0)
  SelIon_mz_R <- NA
  SelIon_tol_R <- NA
  SelIon_mz_G <- NA
  SelIon_tol_G <- NA
  SelIon_mz_B <- NA
  SelIon_tol_B <- NA
  CurrentSelTool <- "Zoom" #Stores the curren state of the sel tool, can be: Zoom, Red, Green and Blue
  MAX_SPECTRA_LIMIT <- max_spectra_limit #Maximum number of spectra that can be added
  ReDraw <- F #Signal when spectra must be redraw
  MAX_MASS_SEL_RANGE <- 4 #Max range of masses to allow selection (in Da)
  ReDrawRedImg <- F
  ReDrawGreenImg <- F
  ReDrawBlueImg <- F

  #Stop gtimer if widget is distroyed
  Widget_Disposed <- function (evt, ...)
  {
    #cat("Stopping spectraWidget draw timer\n")
    this$redrawTimer$stop_timer()
  }

  #Create spectrum object
  #A spectrum is created with a defined mass and intensity. This two params are not modifiable.
  #The original color must be also specified, but this param can be changed after.
  #A spectrum can also provide peaks. This can be added/modified latter.
  CreateSpectrumObj <- function (mass, intensity, color, mass_peaks = NULL, intensity_peaks = NULL, active = T )
  {
    return(list(mass = mass, intensity = intensity, color = color, mass_peaks = mass_peaks, intensity_peaks = intensity_peaks, enabled = active ))
  }

  #Set spectrum color
  SetSpectrumColors <- function(spcObj, color)
  {
    spcObj$color <- color
    return(spcObj)
  }

  #Set spectrum peaks
  SetSpectrumPeaks <- function(spcObj, mass_peaks, intensity_peaks )
  {
    spcObj$mass_peaks <- mass_peaks
    spcObj$intensity_peaks <- intensity_peaks
    return(spcObj)
  }

  #Set new color to the internal data list of spectra
  ChangeSpectrumColor <- function(name, color )
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    this$spectra_data[[as.character(name)]]<-this$SetSpectrumColors(this$spectra_data[[as.character(name)]], color)
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Set new peaks to the internal data list of spectra
  ChangeSpectrumPeaks <- function(name, peaks_mass, peaks_intensity )
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    this$spectra_data[[as.character(name)]]<-this$SetSpectrumPeaks(this$spectra_data[[as.character(name)]], peaks_mass, peaks_intensity)
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Set enabled state of a spectrum, if enabled then it is visble
  SetSpectrumEnabled <- function(name, enabled)
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    this$spectra_data[[as.character(name)]]$enabled <- enabled
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Return the spectrum enabled state from a given name
  GetSpectrumEnabled <- function(name)
  {
    if( is.null(this$spectra_data[[as.character(name)]]$mass))
    {
      cat(paste("No spectrum with the provided name:", name, "\n"))
      return()
    }

    return(his$spectra_data[[as.character(name)]]$enabled)
  }

  #Return all names of the current spectra list
  GetSpectraInfo <- function()
  {
    spcNames <- names(this$spectra_data)
    spcColors <- c()
    for( i in 1:length(spcNames))
    {
      spcColors <- c(spcColors, this$spectra_data[[as.character(spcNames[i])]]$color)
    }
    return(data.frame( names = spcNames, colors = spcColors))
  }

  #Return the current plotted mass range
  GetPlottedMassRange <- function()
  {
    return(list( mz_min = this$mz_lim[1], mz_max = this$mz_lim[2]))
  }

  #Set a mew plotted mass range
  SetPlottedMassRange <- function(mz_min, mz_max)
  {
    if(length(this$spectra_data) == 0) return()
    this$mz_lim[1] <- mz_min
    this$mz_lim[2] <- mz_max
    this$AutoZoomIntensity()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Add spectrum data
  AddSpectra <- function( mass_data, intensity_data, mass_peaks = NULL, intensity_peaks = NULL, col = "", name = "", add_enabled = T )
  {
    if( length(this$spectra_data) >= this$MAX_SPECTRA_LIMIT)
    {
      gWidgets2::gmessage(paste("The limit number of spectra (",this$MAX_SPECTRA_LIMIT , ") has been reached. Remove some spectrum to add a new one.", sep =""), icon = "error")
      return()
    }

    if(col == "") #Set a default color if color is not provided
    {
      col<-as.character(sample(rainbow(100), 1))
    }
    if(name == "" ) #Provide a unique name if name is not specified
    {
      while(name == "")
      {
        name <- as.character(sample(1:1e3, 1))
        if( length(which( names(this$spectra_data) == name)))
        {
          name <- ""
        }
      }
    }
    this$spectra_data[[as.character(name)]]<-this$CreateSpectrumObj(mass_data, intensity_data, col, mass_peaks, intensity_peaks, add_enabled)

    this$data_mass_range <- c(min(sapply(this$spectra_data, function(x){ return( x$mass[1] ) })),max(sapply(this$spectra_data, function(x){ return( x$mass[length(x$mass)] ) })))

    #Set Spin_massSel range properly
    gWidgets2::blockHandlers(this$Spin_massSel)
    RGtk2::gtkSpinButtonSetRange( gWidgets2::getToolkitWidget( this$Spin_massSel),  this$data_mass_range[1],  this$data_mass_range[2] )
    gWidgets2::unblockHandlers(this$Spin_massSel)

    if(length(this$spectra_data) == 1) #First spectrum added, so zoom properly
    {
      this$mz_lim <- this$data_mass_range
    }
    this$AutoZoomIntensity()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Set ref mass data, ref masses will be ploted as vertical dashed lines
  SetRefMass <- function( mass_data )
  {
    this$ref_mass <- mass_data
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Remove spectrum data
  RmSpectra <- function( name )
  {
    this$spectra_data<- this$spectra_data[-(which(names(this$spectra_data) == as.character(name)))]

    #Recompute mz_limits, data_mass range to scale acording the removed spectra
    if(length(this$spectra_data) > 0)
    {
      this$mz_lim <- c(min(sapply(this$spectra_data, function(x){ return( x$mass[1] ) })),max(sapply(this$spectra_data, function(x){ return( x$mass[length(x$mass)] ) })))
      this$in_lim <- c(0 ,max(sapply(this$spectra_data, function(x){ return( max(x$intensity)) }))*1.1)
      this$data_mass_range <- this$mz_lim

      #Set Spin_massSel range properly
      gWidgets2::blockHandlers(this$Spin_massSel)
      RGtk2::gtkSpinButtonSetRange( gWidgets2::getToolkitWidget( this$Spin_massSel),  this$data_mass_range[1],  this$data_mass_range[2] )
      gWidgets2::unblockHandlers(this$Spin_massSel)
    }
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Clear all spectrum data
  ClearSpectra <- function( )
  {
    this$spectra_data<-list()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Redraw MS image on parent widget in a given channel
  ReDrawParentMSI <- function()
  {
    if(this$ReDrawRedImg) #Red
    {
      this$clicFun(1, this$SelIon_mz_R, this$SelIon_tol_R)
    }
    if(this$ReDrawGreenImg) #Green
    {
      this$clicFun(2, this$SelIon_mz_G, this$SelIon_tol_G)
    }
    if(this$ReDrawBlueImg) #Blue
    {
      this$clicFun(3, this$SelIon_mz_B, this$SelIon_tol_B)
    }

    mz_sel_spin <- NULL
    mz_tol_spin <- NULL
    if(this$CurrentSelTool == "Red" )
    {
      mz_sel_spin <- this$SelIon_mz_R
      mz_tol_spin <- this$SelIon_tol_R
    }

    if(this$CurrentSelTool == "Green" )
    {
      mz_sel_spin <- this$SelIon_mz_G
      mz_tol_spin <- this$SelIon_tol_G
    }

    if(this$CurrentSelTool == "Blue" )
    {
      mz_sel_spin <- this$SelIon_mz_B
      mz_tol_spin <- this$SelIon_tol_B
    }

    #Set Spin_massSel range properly
    if( !is.null(mz_sel_spin) && !is.null(mz_tol_spin))
    {
      gWidgets2::blockHandlers(this$Spin_massSel)
      gWidgets2::blockHandlers(this$Spin_TolSel)
      gWidgets2::svalue(this$Spin_massSel) <- mz_sel_spin
      gWidgets2::svalue(this$Spin_TolSel) <- mz_tol_spin
      RGtk2::gtkSpinButtonSetIncrements( gWidgets2::getToolkitWidget( this$Spin_massSel), mz_tol_spin, mz_tol_spin)
      gWidgets2::unblockHandlers(this$Spin_massSel)
      gWidgets2::unblockHandlers(this$Spin_TolSel)
    }

    this$ReDrawRedImg <- F
    this$ReDrawGreenImg <- F
    this$ReDrawBlueImg <- F
  }

  #Redraw ggraph with interpolation using a timer
  ReDrawByTimer <- function( data )
  {
    if( this$ReDraw )
    {
      #Redraw parent MS image if necessari
      if( !is.null( this$clicFun ))
      {
          ReDrawParentMSI()
      }

      #Reduce spectra data size to speed-up ploting
      mass_range <- this$mz_lim
      npoints <- 2*gWidgets2::size(this$plot_device)["width"]

      #Visible before plot() forces the target divice for ploting
      visible(this$plot_device)<-TRUE
      par(mar = c(3.1, 5.1, 0.5, 0.5), cex = 0.7, xaxs = "i", yaxs = "i")

      #Init Plot
      #i_axt<-pretty(this$in_lim[1]:this$in_lim[2], n = 5)
      i_axt<-seq(from = this$in_lim[1], to = this$in_lim[2], length.out = 5)
      plot(x=0, xlim = this$mz_lim, ylim = this$in_lim, type = "n", xlab = "", ylab ="", yaxt ="n")
      axis(2, at = i_axt, labels = sprintf("%.2e",i_axt), las = 1)

      #Draw ref masses as vertical lines
      if(!is.null(this$ref_mass))
      {
        abline(v = this$ref_mass, col = "grey", lty = 2)
      }

      #Plot selection range Red
      if(!is.na(this$SelIon_mz_R) && !is.na(this$SelIon_tol_R))
      {
        rect(xleft = this$SelIon_mz_R - this$SelIon_tol_R, xright = this$SelIon_mz_R + this$SelIon_tol_R, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "lightsalmon", border = "red3")
      }

      #Plot selection range Green
      if(!is.na(this$SelIon_mz_G) && !is.na(this$SelIon_tol_G))
      {
        rect(xleft = this$SelIon_mz_G - this$SelIon_tol_G, xright = this$SelIon_mz_G + this$SelIon_tol_G, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "lightgreen", border = "green3")
      }

      #Plot selection range Blue
      if(!is.na(this$SelIon_mz_B) && !is.na(this$SelIon_tol_B))
      {
        rect(xleft = this$SelIon_mz_B - this$SelIon_tol_B, xright = this$SelIon_mz_B + this$SelIon_tol_B, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "lightblue", border = "blue3")
      }

      if(length(this$spectra_data) > 0)
      {
        for(li in 1:length(this$spectra_data))
        {
          if( this$spectra_data[[li]]$enabled)
          {
            #Work with a copy which is much more fastter than accesing a pointer trought R.oo package
            mz_copy<-this$spectra_data[[li]]$mass
            int_copy <- this$spectra_data[[li]]$intensity

            imin<-which(mz_copy <=  mass_range[1], arr.ind = T)
            imax<-which(mz_copy >=  mass_range[2], arr.ind = T)

            #Limiting to real mass range
            if(length(imin) == 0)
            {
              imin<-1
            }
            if(length(imax) == 0)
            {
              imax<-length(mz_copy)
            }
            imin<-imin[length(imin)]
            imax<-imax[1]

            #Interpolation
            if(imax - imin +1 > npoints)
            {
              new_axes<-approx(imin:imax, mz_copy[imin:imax], n = npoints) #mz axes interpolation
              mzData<-new_axes$y
              dst.vect<-(new_axes$x[-1] - new_axes$x[-length(new_axes$x)])/2
              dst.vect<-c(0, dst.vect, 0) #Add sides
              inData<-c()

              for(i in 1:length(mzData))
              {
                inMax<-max(int_copy[ceiling(new_axes$x[i] - dst.vect[i]):floor(new_axes$x[i] + dst.vect[i+1])])
                inData<-c(inData, inMax)
              }
            }
            else
            {
              mzData<-mz_copy[imin:imax]
              inData<-int_copy[imin:imax]
            }

            lines(x = mzData, y = inData, col= this$spectra_data[[li]]$color)
            if(length(mzData) <= 200)
            {
              points(x = mzData, y = inData, col= this$spectra_data[[li]]$color, pch = 20)
            }

            #Plot labels
            if( !is.null(this$spectra_data[[li]]$mass_peaks) && !is.null(this$spectra_data[[li]]$intensity_peaks) )
            {
              #cat(paste("Mz peaks:", this$spectra_data[[li]]$mass_peaks, "int peaks:", this$spectra_data[[li]]$intensity_peaks, "\n"))
              pk_lbl <- sprintf("%.4f", this$spectra_data[[li]]$mass_peaks) #4 decimals
              text(x = this$spectra_data[[li]]$mass_peaks, y = this$spectra_data[[li]]$intensity_peaks, labels = pk_lbl, pos = 3, cex = 0.8, offset = 1.1)
              un_lbl<-sapply(pk_lbl, function(x) { paste(rep("_", nchar(x)), collapse = "") })
              text(x = this$spectra_data[[li]]$mass_peaks, y = this$spectra_data[[li]]$intensity_peaks, labels = un_lbl, pos = 3, cex = 0.8, offset = 1.1)
              text(x = this$spectra_data[[li]]$mass_peaks, y = this$spectra_data[[li]]$intensity_peaks, labels = rep("|", length(pk_lbl)), pos = 3, cex = 0.8)
            }
          }
        }
      }
    }
    this$SetStateAccordingSelTool()
    this$ReDraw <- F #Reset the redraw signaling
  }

  #OpenTXT
  OpenTXT <- function( evt, ... )
  {
    fname<-file.choose()
    spect<-read.table(fname, header = F, sep = " ")
    this$AddSpectra( spect[,1], spect[,2], name = basename(fname))
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Reset Zoom button clicked
  ZoomResetClicked <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    this$mz_lim <- this$data_mass_range
    this$in_lim <- c(0 ,max(sapply(this$spectra_data, function(x){ return( max(x$intensity)) }))*1.1)
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Auto Zoomin Mz axis
  ZoomMzClicked <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    this$mz_lim <- this$data_mass_range
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Auto Zomming Intensity Axis
  ZoomInClicked <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    this$AutoZoomIntensity()
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Intensity auto-zoom
  AutoZoomIntensity <- function( )
  {
    if(length(this$spectra_data) == 0) return()
    this$in_lim <- c(0 ,max(unlist(lapply(this$spectra_data, function(x){ max(x$intensity[ which( x$mass >= this$mz_lim[1] & x$mass <= this$mz_lim[2], arr.ind = T ) ]) })))*1.1)
  }

  #Mz zoom in a defined range with autoscaling intensity
  ZoomMzRange <- function(mzLow, mzHigh)
  {
    #Mz zooming
    mzLow<-max(mzLow, this$data_mass_range[1])
    mzHigh<-min(mzHigh, this$data_mass_range[2])
    this$mz_lim<- c(mzLow, mzHigh)

    #Intensity zoom
    this$AutoZoomIntensity()

    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  #Grab Mouse Selection Changes on plot
  OnSelection <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()
    if(this$CurrentSelTool == "Zoom")
    {
      if(abs(evt$x[1] - evt$x[2]) > 1)
      {
        #Mz zooming
        this$ZoomMzRange(min(evt$x), max(evt$x))
      }
    }
    else
    {
      top_left <- min(evt$x)
      top_right <- max(evt$x)
      mz_tol <-(top_right - top_left)/2
      mz_sel <- top_left + mz_tol
      mz_tol<-round(mz_tol, digits = 2)
      mz_sel<-round(mz_sel, digits = 2)
      mz_tol<-max(mz_tol, 0)

      #Use the tolerance spin if the spectrum was just clicked
      if(mz_tol == 0)
      {
        mz_tol <- gWidgets2::svalue(this$Spin_TolSel)
      }

      #Limit selection to 5 Da to avoid selecting large parts of spectra and filling RAM
      if(mz_tol*2 > this$MAX_MASS_SEL_RANGE )
      {
        cat(paste("Ion selection in a range of", mz_tol*2, "Da has been aborted. To large data sector.\n"))
        return(TRUE)
      }

      if(this$CurrentSelTool == "Red")
      {
        this$SelIon_mz_R <- mz_sel
        this$SelIon_tol_R <- mz_tol
        this$ReDrawRedImg <- T
      }
      else if ( this$CurrentSelTool == "Green")
      {
        this$SelIon_mz_G <- mz_sel
        this$SelIon_tol_G <- mz_tol
        this$ReDrawGreenImg <- T
      }
      else if ( this$CurrentSelTool == "Blue")
      {
        this$SelIon_mz_B <- mz_sel
        this$SelIon_tol_B <- mz_tol
        this$ReDrawBlueImg <- T
      }
      #Plot selection in spectra
      this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
    }

    return(TRUE)
  }

  #Grab mouse cursor on plot
  OnMouseMotion <- function( evt, ... )
  {
    if(length(this$spectra_data) == 0) return()

    #Update pointer
    this$PointerCoords<-c(evt$x, evt$y)

    #Update Labels
    if(evt$x >= this$mz_lim[1] &&
       evt$x <= this$mz_lim[2] &&
       evt$y >= this$in_lim[1] &&
       evt$y <= this$in_lim[2])
    {
      mz_txt<-sprintf(paste("%-",this$LABEL_LENGTH, ".4f", sep = ""), this$PointerCoords[1])
      in_txt<-sprintf(paste("%-",this$LABEL_LENGTH, ".2e", sep = ""), this$PointerCoords[2])
    }
    else
    {
      in_txt<-mz_txt<-paste(rep(" ", this$LABEL_LENGTH), collapse = "")
    }
    this$lbl_mz_coords$set_value(mz_txt)
    this$lbl_in_coords$set_value(in_txt)
    return(TRUE)
  }

  #Zoom on spectra plot handler
  ScrollEventOnSpectra <- function( obj, evt, ... )
  {
    if(length(this$spectra_data) == 0) return()

    dir<- as.double(evt$direction)
    dir<-dir*(-2) + 1 # 1 is up, -1 is down
    mz_min<-this$data_mass_range[1]
    mz_max<-this$data_mass_range[2]

    if(this$MouseWheelFunction == 0)
    {
      #Mz scrolling
      range<-abs(this$mz_lim[2] - this$mz_lim[1])
      mz_lim<-this$mz_lim - dir*range*0.05
      range<-abs(mz_lim[2] - mz_lim[1])
      if(mz_lim[1] < mz_min)
      {
        #Clip to min
        mz_lim<-c(mz_min, mz_min + range)
      }
      if(mz_lim[2] > mz_max)
      {
        #Clip to max
        mz_lim<-c(mz_max - range, mz_max)
      }
      if(range > abs(mz_max - mz_min))
      {
        #out of range! clip in full spectra
        mz_lim<-c(mz_min, mz_max)
      }
      this$mz_lim<- mz_lim
    }

    else if(this$MouseWheelFunction == 1)
    {
      #Mz zooming
      pointer.x<-this$PointerCoords[1]
      top_left <- pointer.x -  abs(0.1 + dir)*(pointer.x - this$mz_lim[1])
      top_right <- pointer.x + abs(0.1 + dir)*(this$mz_lim[2] - pointer.x)
      top_left<-max(top_left, mz_min)
      top_right<-min(top_right, mz_max)
      this$mz_lim<- c(top_left, top_right)
    }
    else
    {
      #Intensity scaling
      range<-abs(this$in_lim[2] - this$in_lim[1])
      range<-range - dir*range*0.05
      this$in_lim<- c(0, range)
    }

    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt

    return(TRUE) #The scroll event requires this return
  }

  #Key Press handler
  OnKeyPress <- function( evt, ... )
  {
    if(evt[[4]][["keyval"]] == 65507)
    {
      if(this$MouseWheelFunction == 0)
      {
        this$MouseWheelFunction<-1
      }
    }
    else if(evt[[4]][["keyval"]] == 65505)
    {
      if(this$MouseWheelFunction == 0)
      {
        this$MouseWheelFunction<-2
      }
    }

    return(FALSE) #The key event requires this return to allow keyboard continue working for other widgets as spin buttons
  }

  #Key Release handler
  OnKeyRelease <- function( evt, ... )
  {
    if(evt[[4]][["keyval"]] == 65507)
    {
      if(this$MouseWheelFunction == 1)
      {
        this$MouseWheelFunction<-0
      }
    }
    else if(evt[[4]][["keyval"]] == 65505)
    {
      if(this$MouseWheelFunction == 2)
      {
        this$MouseWheelFunction<-0
      }
    }
    return(FALSE) #The key event requires this return to allow keyboard continue working for other widgets as spin buttons
  }

  #Windows lost focuts, used to restore zoom status
  OnLostFocus <- function( evt, ... )
  {
    this$MouseWheelFunction<-0
    return(TRUE) #This event requires this return
  }

  #Clear all spectra button click
  ClearSpectraClicked <- function( evt, ... )
  {
    this$ClearSpectra()
  }

  #Zoom tool has been selected
  ZoomToolSel <- function( ... )
  {
    if( gWidgets2::svalue(this$Btn_ZoomTool))
    {
      this$CurrentSelTool <- "Zoom"
      if(!is.null(this$Btn_SelRedTool))
      {
        gWidgets2::svalue(this$Btn_SelRedTool) <- F
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        gWidgets2::svalue(this$Btn_SelGreenTool) <- F
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        gWidgets2::svalue(this$Btn_SelBlueTool) <- F
      }

      #Set ion manually selection visibility
      gWidgets2::visible(Lbl_massSel) <- F
      gWidgets2::visible(Spin_massSel) <- F
      gWidgets2::visible(Lbl_TolSel) <- F
      gWidgets2::visible(Spin_TolSel) <- F
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      if(!is.null(this$Btn_SelRedTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelRedTool) | bTest
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelGreenTool) | bTest
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelBlueTool) | bTest
      }
      if( bTest == F)
      {
        gWidgets2::svalue(this$Btn_ZoomTool) <- T
      }
    }
    this$SetStateAccordingSelTool()
  }

  #Red tool selected
  RedToolSel <- function( ... )
  {
    if( gWidgets2::svalue(this$Btn_SelRedTool))
    {
      this$CurrentSelTool <- "Red"
      gWidgets2::svalue(this$Btn_ZoomTool) <- F
      if(!is.null(this$Btn_SelGreenTool))
      {
        gWidgets2::svalue(this$Btn_SelGreenTool) <- F
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        gWidgets2::svalue(this$Btn_SelBlueTool) <- F
      }

      #Set Spin_massSel range properly
      gWidgets2::blockHandlers(this$Spin_massSel)
      gWidgets2::blockHandlers(this$Spin_TolSel)
      gWidgets2::svalue(this$Spin_massSel) <- this$SelIon_mz_R
      gWidgets2::svalue(this$Spin_TolSel) <- this$SelIon_tol_R
      RGtk2::gtkSpinButtonSetIncrements( gWidgets2::getToolkitWidget( this$Spin_massSel), this$SelIon_tol_R, this$SelIon_tol_R)
      gWidgets2::unblockHandlers(this$Spin_massSel)
      gWidgets2::unblockHandlers(this$Spin_TolSel)

      #Set ion manually selection visibility
      gWidgets2::visible(Lbl_massSel) <- T
      gWidgets2::visible(Spin_massSel) <- T
      gWidgets2::visible(Lbl_TolSel) <- T
      gWidgets2::visible(Spin_TolSel) <- T
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      bTest <- gWidgets2::svalue(this$Btn_ZoomTool) | bTest
      if(!is.null(this$Btn_SelGreenTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelGreenTool) | bTest
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelBlueTool) | bTest
      }
      if( bTest == F)
      {
        gWidgets2::svalue(this$Btn_SelRedTool) <- T
      }
    }
    this$SetStateAccordingSelTool()
  }

  #Green tool selected
  GreenToolSel <- function( ... )
  {
    if( gWidgets2::svalue(this$Btn_SelGreenTool))
    {
      this$CurrentSelTool <- "Green"
      gWidgets2::svalue(this$Btn_ZoomTool) <- F
      if(!is.null(this$Btn_SelRedTool))
      {
        gWidgets2::svalue(this$Btn_SelRedTool) <- F
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        gWidgets2::svalue(this$Btn_SelBlueTool) <- F
      }

      #Set Spin_massSel range properly
      gWidgets2::blockHandlers(this$Spin_massSel)
      gWidgets2::blockHandlers(this$Spin_TolSel)
      gWidgets2::svalue(this$Spin_massSel) <- this$SelIon_mz_G
      gWidgets2::svalue(this$Spin_TolSel) <- this$SelIon_tol_G
      RGtk2::gtkSpinButtonSetIncrements( gWidgets2::getToolkitWidget( this$Spin_massSel), this$SelIon_tol_G, this$SelIon_tol_G)
      gWidgets2::unblockHandlers(this$Spin_massSel)
      gWidgets2::unblockHandlers(this$Spin_TolSel)

      #Set ion manually selection visibility
      gWidgets2::visible(Lbl_massSel) <- T
      gWidgets2::visible(Spin_massSel) <- T
      gWidgets2::visible(Lbl_TolSel) <- T
      gWidgets2::visible(Spin_TolSel) <- T
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      bTest <- gWidgets2::svalue(this$Btn_ZoomTool) | bTest
      if(!is.null(this$Btn_SelRedTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelRedTool) | bTest
      }
      if(!is.null(this$Btn_SelBlueTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelBlueTool) | bTest
      }
      if( bTest == F)
      {
        gWidgets2::svalue(this$Btn_SelGreenTool) <- T
      }
    }
    this$SetStateAccordingSelTool()
  }

  #Blue tool selected
  BlueToolSel <- function( ... )
  {
    if( gWidgets2::svalue(this$Btn_SelBlueTool))
    {
      this$CurrentSelTool <- "Blue"
      gWidgets2::svalue(this$Btn_ZoomTool) <- F
      if(!is.null(this$Btn_SelRedTool))
      {
        gWidgets2::svalue(this$Btn_SelRedTool) <- F
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        gWidgets2::svalue(this$Btn_SelGreenTool) <- F
      }

      #Set Spin_massSel range properly
      gWidgets2::blockHandlers(this$Spin_massSel)
      gWidgets2::blockHandlers(this$Spin_TolSel)
      gWidgets2::svalue(this$Spin_massSel) <- this$SelIon_mz_B
      gWidgets2::svalue(this$Spin_TolSel) <- this$SelIon_tol_B
      RGtk2::gtkSpinButtonSetIncrements( gWidgets2::getToolkitWidget( this$Spin_massSel), this$SelIon_tol_B, this$SelIon_tol_B)
      gWidgets2::unblockHandlers(this$Spin_massSel)
      gWidgets2::unblockHandlers(this$Spin_TolSel)

      #Set ion manually selection visibility
      gWidgets2::visible(Lbl_massSel) <- T
      gWidgets2::visible(Spin_massSel) <- T
      gWidgets2::visible(Lbl_TolSel) <- T
      gWidgets2::visible(Spin_TolSel) <- T
    }
    else #Check if all ara false and avoid such situation
    {
      bTest <- F
      bTest <- gWidgets2::svalue(this$Btn_ZoomTool) | bTest
      if(!is.null(this$Btn_SelRedTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelRedTool) | bTest
      }
      if(!is.null(this$Btn_SelGreenTool))
      {
        bTest <- gWidgets2::svalue(this$Btn_SelGreenTool) | bTest
      }
      if( bTest == F)
      {
        gWidgets2::svalue(this$Btn_SelBlueTool) <- T
      }
    }
    this$SetStateAccordingSelTool()
  }

  SetStateAccordingSelTool <- function ()
  {
    if (this$CurrentSelTool == "Zoom")
    {
      c=RGtk2::gdkCursorNew("GDK_SIZING")
    }
    else
    {
      c=RGtk2::gdkCursorNew("GDK_BASED_ARROW_DOWN")
    }
    gWidgets2::getToolkitWidget(this$plot_device)$getWindow()$setCursor(c)
  }

  #Set the tool to use externally, can be: Zoom, Red, Green or Blue
  SetActiveTool <- function (tool)
  {
    if( tool == "Zoom")
    {
      gWidgets2::svalue(this$Btn_ZoomTool) <- T
    }
    else if( tool == "Red")
    {
      if(!is.null(this$Btn_SelRedTool))
      {
        gWidgets2::svalue(this$Btn_SelRedTool) <- T
      }
    }
    else if( tool == "Green")
    {
      if(!is.null(this$Btn_SelGreenTool))
      {
        gWidgets2::svalue(this$Btn_SelGreenTool) <- T
      }
    }
    else if( tool == "Blue")
    {
      if(!is.null(this$Btn_SelBlueTool))
      {
        gWidgets2::svalue(this$Btn_SelBlueTool) <- T
      }
    }
  }

  MassSelSpinChanged <- function(...)
  {
    if(this$CurrentSelTool == "Red")
    {
      this$SelIon_mz_R <- gWidgets2::svalue(this$Spin_massSel)
      this$ReDrawRedImg <- T
    }
    else if ( this$CurrentSelTool == "Green")
    {
      this$SelIon_mz_G <- gWidgets2::svalue(this$Spin_massSel)
      this$ReDrawGreenImg <- T
    }
    else if ( this$CurrentSelTool == "Blue")
    {
      this$SelIon_mz_B <- gWidgets2::svalue(this$Spin_massSel)
      this$ReDrawBlueImg <- T
    }

    #Set zoom range to see the selected ion centered
    currentRange <- this$mz_lim[2] - this$mz_lim[1]
    this$ZoomMzRange( gWidgets2::svalue(this$Spin_massSel) - (currentRange/2), gWidgets2::svalue(this$Spin_massSel) + (currentRange/2))

    #Plot selection in spectra
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  TolSelSpinChanged <- function(...)
  {
    if(this$CurrentSelTool == "Red")
    {
      this$SelIon_tol_R <-  gWidgets2::svalue(this$Spin_TolSel)
      this$ReDrawRedImg <- T
    }
    else if ( this$CurrentSelTool == "Green")
    {
      this$SelIon_tol_G <-  gWidgets2::svalue(this$Spin_TolSel)
      this$ReDrawGreenImg <- T
    }
    else if ( this$CurrentSelTool == "Blue")
    {
      this$SelIon_tol_B <-  gWidgets2::svalue(this$Spin_TolSel)
      this$ReDrawBlueImg <- T
    }
    #Plot selection in spectra
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

  SetSelectedMassTol <- function (channel, mass, tol)
  {
    if(channel == 1)
    {
      this$SelIon_mz_R <- mass
      this$SelIon_tol_R <- tol
    }
    if(channel == 2)
    {
      this$SelIon_mz_G <- mass
      this$SelIon_tol_G <- tol
    }
    if(channel == 3)
    {
      this$SelIon_mz_B <- mass
      this$SelIon_tol_B <- tol
    }
    this$ReDraw <- T #Signal a redraw request, redraw will be performed on the next timer interrupt
  }

    #Build GUI
  Grp_Top <- gWidgets2::ggroup(horizontal = F, container = this$parent, fill = T, expand = T)
  Grp_Buttons<- gWidgets2::ggroup(horizontal = T, container = Grp_Top, fill = F, expand = F)
  Btn_reset_zoom<- gWidgets2::gbutton(text = "Reset Zoom", handler = this$ZoomResetClicked, action = this, container = Grp_Buttons)
  Btn_auto_zoom_mz<- gWidgets2::gbutton(text = "Auto m/z", handler = this$ZoomMzClicked, action = this, container = Grp_Buttons)
  Btn_auto_zoom_in<- gWidgets2::gbutton(text = "Auto Intensity", handler = this$ZoomInClicked, action = this, container = Grp_Buttons)
  Btn_RemoveSpectra<- gWidgets2::gbutton(text = "Clear all", handler = this$ClearSpectraClicked, action = this, container = Grp_Buttons)
  Btn_ZoomTool <- gWidgets2::gcheckbox("Zoom", checked = T, handler = this$ZoomToolSel, container = Grp_Buttons, use.togglebutton = F)
  gWidgets2::visible(Btn_ZoomTool) <- display_sel_red | display_sel_green | display_sel_blue #Only dispaly the zoom tool selector if at leas one sel. ion is visible

  if( display_sel_red )
  {
    Btn_SelRedTool <- gWidgets2::gcheckbox("SelRed", checked = F, handler = this$RedToolSel, container = Grp_Buttons, use.togglebutton = F)
    .setCheckBoxText(Btn_SelRedTool, "Sel.Red", background = NULL, foreground = "darkred", font_size = NULL, font_weight = "heavy")
  }
  else
  {
    Btn_SelRedTool <- NULL
  }

  if( display_sel_green )
  {
    Btn_SelGreenTool <- gWidgets2::gcheckbox("GreenRed", checked = F, handler = this$GreenToolSel, container = Grp_Buttons, use.togglebutton = F)
    .setCheckBoxText(Btn_SelGreenTool, "Sel.Green", background = NULL, foreground = "darkgreen", font_size = NULL, font_weight = "heavy")
  }
  else
  {
    Btn_SelGreenTool <- NULL
  }

  if( display_sel_blue )
  {
    Btn_SelBlueTool <- gWidgets2::gcheckbox("SelBlue", checked = F, handler = this$BlueToolSel, container = Grp_Buttons, use.togglebutton = F)
    .setCheckBoxText(Btn_SelBlueTool, "Sel.Blue", background = NULL, foreground = "darkblue", font_size = NULL, font_weight = "heavy")
  }
  else
  {
    Btn_SelBlueTool <- NULL
  }

  #Dislay mass and tolerance spin boxes
  Lbl_massSel <- gWidgets2::glabel("m/z:", container = Grp_Buttons)
  Spin_massSel <- gWidgets2::gspinbutton( from = 0, to = 1, value = 0, digits = 4, by = 0.1, container =  Grp_Buttons, handler = this$MassSelSpinChanged)
  Lbl_TolSel <- gWidgets2::glabel("+/-", container = Grp_Buttons)
  Spin_TolSel <- gWidgets2::gspinbutton( from = 0, to = MAX_MASS_SEL_RANGE, value = 0.1, digits = 4, by = 0.01, container =  Grp_Buttons, handler = this$TolSelSpinChanged)
  gWidgets2::visible(Lbl_massSel) <- F
  gWidgets2::visible(Spin_massSel) <- F
  gWidgets2::visible(Lbl_TolSel) <- F
  gWidgets2::visible(Spin_TolSel) <- F

  gWidgets2::addSpring(Grp_Buttons)
  if(showOpenFileButton)
  {
    Btn_file_open<-gWidgets2::gbutton(text = "Open spectra TXT", handler = this$OpenTXT, action = this, container = Grp_Buttons)
  }
  rm(showOpenFileButton)

  this$plot_device <- gWidgets2::ggraphics(  )
  gWidgets2::size( this$plot_device )<- c(-1, 170)
  gWidgets2::add(obj = Grp_Top, child = this$plot_device, expand = T, fill = T)

  Grp_BottomLabel<-gWidgets2::ggroup(horizontal = T, container = Grp_Top)
  lbl_help_info<-gWidgets2::glabel(" Mouse wheel -> m/z scroll\n Ctrl + Mouse wheel -> m/z zooming\n Shift + Mouse Wheel -> intensity scaling", container = Grp_BottomLabel)
  gWidgets2::addSpring(Grp_BottomLabel)
  lbl_mz<-gWidgets2::glabel("m/z:", container = Grp_BottomLabel)
  this$lbl_mz_coords<-gWidgets2::glabel(paste(rep(" ", this$LABEL_LENGTH), collapse = ""), container = Grp_BottomLabel)
  lbl_int<-gWidgets2::glabel("Intensity:", container = Grp_BottomLabel)
  this$lbl_in_coords<-gWidgets2::glabel(paste(rep(" ", this$LABEL_LENGTH), collapse = ""), container = Grp_BottomLabel)

  font(this$lbl_mz_coords)<-list(family = "monospace", weight = "light", size = 8)
  font(this$lbl_in_coords)<-list(family = "monospace", weight = "light", size = 8)
  font(lbl_mz)<-list(family = "monospace", weight = "light", size = 8)
  font(lbl_int)<-list(family = "monospace", weight = "light", size = 8)
  font(lbl_help_info)<-list(family = "monospace", weight = "light", size = 8)

  #Signal handlers
  scroll_event_id <- RGtk2::gSignalConnect( gWidgets2::getToolkitWidget(this$plot_device),  signal = "scroll-event", f = this$ScrollEventOnSpectra, data = this )
  gWidgets2::addHandler(this$top_window, signal = "key-press-event", handler = this$OnKeyPress, action = this)
  gWidgets2::addHandler(this$top_window, signal = "key-release-event", handler = this$OnKeyRelease, action = this)
  gWidgets2::addHandler(this$top_window, signal = "focus-out-event", handler = this$OnLostFocus, action = this)
  gWidgets2::addHandlerMouseMotion(this$plot_device, handler = this$OnMouseMotion)
  gWidgets2::addHandlerSelectionChanged(this$plot_device, handler = this$OnSelection, action = this)
  gWidgets2::addHandlerDestroy( obj = this$top_window, handler = this$Widget_Disposed ) #Connect to widget dispose to stop the draw timer

  #Start the redraw timer
  redrawTimer <- gWidgets2::gtimer(10, this$ReDrawByTimer)

  ## Set the name for the class
  class(this) <- append(class(this),"SpectraPlotWidget")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
