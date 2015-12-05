###A GUI to display spectra in an interactive way


#Plot Sigles Spectra "User Friendly"
#Simply call this function each time you want to plot a mass spectrum and it will be overlayed
plotSpectra<-function( mass = NULL, intensity = NULL, peaks_mass = 0, peaks_intensity = 0, col = "" )
{
  require(gWidgets2)
  require(gWidgets2RGtk2)
  bFirstRun <- F
  if( !exists( x = ".SpectraWidget", mode = "environment") )
  {
    ###TODO: Add a left-list showing all loaded spectra. The list will also allow to disable one of them, deleting, coloring...
    ###TODO: Yes the list will be added here because MSIWindows has its own list. So I don't want to modify spectraWidget

    #Spectra plot windows does not exists, create it
    window_spectra <- gwindow ( "Spectra Plot" , visible = F , width = 750, height = 440) #using window just as a test! Finally it should be a Widget!
    .SpectraWidget<<-.SpectraPlotWidget( parent_widget = window_spectra )
    addHandlerDestroy( obj = window_spectra, handler = .plotSpectraWindow_Disposed )
    window_spectra$widget$present() ##Tot i forzzar aki un raise-up del top widget en RStudio no deixa!
    bFirstRun <- T
  }

  if( !is.null(mass) && !is.null(intensity) )
  {
    .SpectraWidget$AddSpectra( mass_data = mass, intensity_data = intensity,
                             mass_peaks = peaks_mass, intensity_peaks = peaks_intensity, col = col )
    if(bFirstRun)
    {
      .SpectraWidget$ZoomResetClicked()
    }
  }
}

.plotSpectraWindow_Disposed <- function (evt, ...)
{
  rm(.SpectraWidget, envir = .GlobalEnv)
  gc()
}

.SpectraPlotWidget <- function( parent_widget=gwindow ( "Default SpectraPlotWidget" , visible = FALSE ), clicFuntion = NULL, showOpenFileButton = T)
{
  require(cairoDevice)
  require(gWidgets2)
  require(gWidgets2RGtk2)

  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  parent <- parent_widget
  rm(parent_widget)
  clicFun <- clicFuntion
  rm(clicFuntion)
  AnteriorClickMz <- 0
  plot_device <- 0
  radio_buttons <- 0
  lbl_mz_coords <- 0
  lbl_in_coords <- 0
  LABEL_LENGTH <- 10
  MouseWheelFunction <- 0
  spectra_mass <- list() #mz_list,
  spectra_intensity <- list() #in_list,
  peak_mass <- list()
  peak_intensity <- list()
  ref_mass <- NULL
  colour <- c() #c(colour),
  mz_lim <- c() #c(min(mass_data),max(mass_data)),
  in_lim <- c() #c(0 ,max(intensity_data)*1.05),
  data_mass_range <- c() #c(min(mass_data),max(mass_data)),
  PointerCoords <- c(0,0)
  SelIon_mz <- NA
  SelIon_tol <- NA

  #Add spectrum data================================================================================
  AddSpectra <- function( mass_data, intensity_data, mass_peaks = 0, intensity_peaks = 0, col = "" )
  {
    this$spectra_mass[[length(this$spectra_mass) + 1]]<-mass_data
    this$spectra_intensity[[length(this$spectra_intensity) + 1]]<-intensity_data

    if(col == "")
    {
      this$colour<-c(this$colour,sample(colors(), 1))
    }
    else
    {
      this$colour<-c(this$colour, col)
    }

    #Recompute mz_limits, data_mass range to scale acording the added spectra
    #this$mz_lim <- c(min(sapply(this$spectra_mass, min)),max(sapply(this$spectra_mass, max)))
    #this$in_lim <- c(0 ,max(sapply(this$spectra_intensity, max))*1.1)
    #this$data_mass_range <- this$mz_lim

    this$AutoZoomIntensity()
    this$data_mass_range <- c(min(sapply(this$spectra_mass, min)),max(sapply(this$spectra_mass, max)))

    #Add peaklist
    if(mass_peaks != 0 && intensity_peaks != 0)
    {
      this$peak_mass[[length(this$peak_mass) + 1]]<-mass_peaks
      this$peak_intensity[[length(this$peak_intensity) + 1]]<-intensity_peaks
    }

    this$ReDraw()
  }

  #Set ref mass data, ref masses will be ploted as vertical dashed lines=============================
  SetRefMass <- function( mass_data )
  {
    this$ref_mass <- mass_data
    this$ReDraw()
  }

  #Remove spectrum data==============================================================================
  RmSpectra <- function( index )
  {
    this$spectra_mass[[index]]<-NULL
    this$spectra_intensity[[index]]<-NULL
    this$colour<-this$colour[-index]
    this$peak_mass[[index]]<-NULL
    this$peak_intensity[[index]]<-NULL

    #Recompute mz_limits, data_mass range to scale acording the removed spectra
    if(length(this$spectra_mass) > 0)
    {
      this$mz_lim <- c(min(sapply(this$spectra_mass, min)),max(sapply(this$spectra_mass, max)))
      this$in_lim <- c(0 ,max(sapply(this$spectra_intensity, max))*1.1)
      this$data_mass_range <- this$mz_lim
      this$ReDraw()
    }
  }

  #Clear all spectrum data===========================================================================
  ClearSpectra <- function( )
  {
    this$spectra_mass<-list()
    this$spectra_intensity<-list()
    this$colour<-this$colour<-c()
    this$peak_mass<-list()
    this$peak_intensity<-list()
  }

  #Redraw ggraph with interpolation, event args must be passed======================================
  ReDraw <- function( ) # TODO hi tenia ... com a param
  {
    if(length(this$spectra_mass) == 0) return()

    #Reduce spectra data size to speed-up ploting
    mass_range <- this$mz_lim
    npoints <- 2*size(this$plot_device)["width"]

    #Visible before plot() forces the target divice for ploting
    visible(this$plot_device)<-TRUE
    par(mar = c(3.1, 3.1, 0.5, 0.5), cex = 0.7, xaxs = "i", yaxs = "i")

    #Init Plot
    plot(x=0, xlim = this$mz_lim, ylim = this$in_lim, type = "n", xlab = "", ylab ="")

    #Draw ref masses as vertical lines
    if(!is.null(this$ref_mass))
    {
      abline(v = this$ref_mass, col = "grey", lty = 2)
    }

    #Plot selection range
    if(!is.na(this$SelIon_mz) && !is.na(this$SelIon_tol))
    {
      rect(xleft = this$SelIon_mz - this$SelIon_tol, xright = this$SelIon_mz + this$SelIon_tol, ybottom = this$in_lim[1], ytop = this$in_lim[2]*0.99, col = "snow2", border = "snow3")
    }

    for(li in 1:length(this$spectra_mass))
    {

      #Work with a copy which is much more fastter than accesing a pointer trought R.oo package
      mz_copy<-this$spectra_mass[[li]]
      int_copy <- this$spectra_intensity[[li]]

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

      #     if(li == 1)
      #     {
      #       plot(x = mzData, y = inData, xlim = this$mz_lim, ylim = this$in_lim, col= this$colour[[li]], type = "l", xlab = "", ylab ="")
      #     }
      #     else
      #     {
      lines(x = mzData, y = inData, col= this$colour[[li]])
      #    }
    }

    #Plot labels
    if(length(this$peak_mass) > 0)
    {
      for(i in 1:length(this$peak_mass))
      {
        pk_lbl <- sprintf("%.4f", this$peak_mass[[i]]) #4 decimals
        text(x = this$peak_mass[[i]], y = this$peak_intensity[[i]], labels = pk_lbl, pos = 3, cex = 0.8, offset = 1.1)
        un_lbl<-sapply(pk_lbl, function(x) { paste(rep("_", nchar(x)), collapse = "") })
        text(x = this$peak_mass[[i]], y = this$peak_intensity[[i]], labels = un_lbl, pos = 3, cex = 0.8, offset = 1.1)
        text(x = this$peak_mass[[i]], y = this$peak_intensity[[i]], labels = rep("|", length(pk_lbl)), pos = 3, cex = 0.8)
      }
    }

    #Update cursor
    this$CheckBox_Changed( ) #TODO hi tenia ... com a param
  }

  #OpenTXT==========================================================================================
  OpenTXT <- function( evt, ... )
  {
    fname<-file.choose()
    spect<-read.table(fname, header = F, sep = " ")
    this$AddSpectra( spect[,1], spect[,2])
    this$ReDraw()
  }

  #Reset Zoom button clicked========================================================================
  ZoomResetClicked <- function( evt, ... )
  {
    if(length(this$spectra_mass) == 0) return()
    this$mz_lim = this$data_mass_range
    this$in_lim <- c(0 ,max(sapply(this$spectra_intensity, max))*1.1)
    this$ReDraw()
  }

  #Auto Zoomin Mz axis==============================================================================
  ZoomMzClicked <- function( evt, ... )
  {
    if(length(this$spectra_mass) == 0) return()
    this$mz_lim = this$data_mass_range
    this$ReDraw()
  }

  #Auto Zomming Intensity Axis======================================================================
  ZoomInClicked <- function( evt, ... )
  {
    if(length(this$spectra_mass) == 0) return()
    this$AutoZoomIntensity()
    this$ReDraw()
  }

  #Auto Change cursor==============================================================================
  CheckBox_Changed <- function( evt, ... )
  {
    if (svalue(this$radio_buttons) == "Zoom")
    {
      c=gdkCursorNew("GDK_SIZING")
    }
    else
    {
      c=gdkCursorNew("GDK_BASED_ARROW_DOWN")
    }
    getToolkitWidget(this$plot_device)$getWindow()$setCursor(c)
  }

  #Intensity auto-zoom==============================================================================
  AutoZoomIntensity <- function( )
  {
    if(length(this$spectra_mass) == 0) return(TRUE)
    max_in<-lapply(this$spectra_mass, function(x){which(x >= this$mz_lim[1] & x <= this$mz_lim[2], arr.ind = T)})
    my_max<-0
    for(i in 1:length(max_in))
    {
      my_max<-max(my_max, this$spectra_intensity[[i]][max_in[[i]]])
    }

    this$in_lim<-c(0, my_max*1.1)
  }

  #Force a plot redraw when resizeing===============================================================
  OnResize <- function( evt, ... )
  {
    this$ReDraw()
    return(TRUE)
  }

  #Grab Mouse Selection Changes on plot=============================================================
  OnSelection <- function( evt, ... )
  {
    if(length(this$spectra_mass) == 0) return(TRUE)
    if(svalue(this$radio_buttons) == "Zoom")
    {
      if(abs(evt$x[1] - evt$x[2]) > 1)
      {
        #Mz zooming
        mz_min<-this$data_mass_range[1]
        mz_max<-this$data_mass_range[2]
        top_left <- min(evt$x)
        top_right <- max(evt$x)
        top_left<-max(top_left, mz_min)
        top_right<-min(top_right, mz_max)
        this$mz_lim<- c(top_left, top_right)

        #Intensity zoom
        this$AutoZoomIntensity()

        this$ReDraw()
      }
    }
    else if(svalue(this$radio_buttons) == "Select")
    {
      top_left <- min(evt$x)
      top_right <- max(evt$x)
      mz_tol <-(top_right - top_left)/2
      mz_sel <- top_left + mz_tol
      mz_tol<-round(mz_tol, digits = 2)
      mz_sel<-round(mz_sel, digits = 2)
      mz_tol<-max(mz_tol, 0)

      if(!is.null(this$clicFun))
      {
        ret <- this$clicFun(mz_sel, mz_tol)
        mz_sel <- ret$selMz
        mz_tol <- ret$selTol
      }
      #Plot selection in spectra
      this$SelIon_mz = mz_sel
      this$SelIon_tol = mz_tol
      this$ReDraw()
    }

    return(TRUE)
  }

  #Grab mouse cursor on plot=======================================================================
  OnMouseMotion <- function( evt, ... )
  {
    if(length(this$spectra_mass) == 0) return(TRUE)

    #Update pointer
    this$PointerCoords<-c(evt$x, evt$y)

    #Update Labels
    if(evt$x >= this$mz_lim[1] &&
       evt$x <= this$mz_lim[2] &&
       evt$y >= this$in_lim[1] &&
       evt$y <= this$in_lim[2])
    {
      mz_txt<-sprintf(paste("%-",this$LABEL_LENGTH, ".2f", sep = ""), this$PointerCoords[1])
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

  #Zoom on spectra plot handler====================================================================
  ScrollEventOnSpectra <- function( evt, ... )
  {
    if(length(this$spectra_mass) == 0) return(TRUE)

    dir<- as.double(evt[[4]][["direction"]]) # 0 is up, 1 is down
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

    this$ReDraw()
    return(TRUE) #The scroll event requires this return
  }

  #Key Press handler===============================================================================
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

    return(TRUE) #The key event requires this return
  }

  #Key Release handler=============================================================================
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
    return(TRUE) #The key event requires this return
  }

  #Windows lost focuts, used to restore zoom status================================================
  OnLostFocus <- function( evt, ... )
  {
    this$MouseWheelFunction<-0
    return(TRUE) #This event requires this return
  }

  #Build GUI======================================================================================
  Grp_Top <- gWidgets2::ggroup(horizontal = F, container = this$parent, fill = T, expand = T)
  Grp_Buttons<- gWidgets2::ggroup(horizontal = T, container = Grp_Top)
  Btn_reset_zoom<- gWidgets2::gbutton(text = "Reset Zoom", handler = this$ZoomResetClicked, action = this, container = Grp_Buttons)
  Btn_auto_zoom_mz<- gWidgets2::gbutton(text = "Auto Zoom m/z", handler = this$ZoomMzClicked, action = this, container = Grp_Buttons)
  Btn_auto_zoom_in<- gWidgets2::gbutton(text = "Auto Zoom Intensity", handler = this$ZoomInClicked, action = this, container = Grp_Buttons)
  this$radio_buttons <- gWidgets2::gradio( c("Zoom","Select"), selected = 1, horizontal = T, handler = this$CheckBox_Changed, action = this, container = Grp_Buttons)

  gWidgets2::addSpring(Grp_Buttons)
  if(showOpenFileButton)
  {
    Btn_file_open<-gWidgets2::gbutton(text = "Open spectra TXT", handler = this$OpenTXT, action = this, container = Grp_Buttons)
  }
  rm(showOpenFileButton)

  this$plot_device <- gWidgets2::ggraphics()
  gWidgets2::add(obj = Grp_Top, child = this$plot_device, expand = T, fill = T)

  Grp_BottomLabel<-ggroup(horizontal = T, container = Grp_Top)
  lbl_help_info<-glabel(" Mouse wheel -> m/z scroll\n Ctrl + Mouse wheel -> m/z zooming\n Shift + Mouse Wheel -> intensity scaling", container = Grp_BottomLabel)
  gWidgets2::addSpring(Grp_BottomLabel)
  lbl_mz<-glabel("m/z:", container = Grp_BottomLabel)
  this$lbl_mz_coords<-glabel(paste(rep(" ", this$LABEL_LENGTH), collapse = ""), container = Grp_BottomLabel)
  lbl_int<-glabel("Intensity:", container = Grp_BottomLabel)
  this$lbl_in_coords<-glabel(paste(rep(" ", this$LABEL_LENGTH), collapse = ""), container = Grp_BottomLabel)

  font(this$lbl_mz_coords)<-list(family = "monospace", weight = "light", size = 8)
  font(this$lbl_in_coords)<-list(family = "monospace", weight = "light", size = 8)
  font(lbl_mz)<-list(family = "monospace", weight = "light", size = 8)
  font(lbl_int)<-list(family = "monospace", weight = "light", size = 8)
  font(lbl_help_info)<-list(family = "monospace", weight = "light", size = 8)

  #Very first spectra plot...
  this$ReDraw()

  #Signal handlers
  gWidgets2::addHandler(this$plot_device, signal = "scroll-event", handler = this$ScrollEventOnSpectra, action = this)
  gWidgets2::addHandler(this$parent, signal = "key-press-event", handler = this$OnKeyPress, action = this)
  gWidgets2::addHandler(this$parent, signal = "key-release-event", handler = this$OnKeyRelease, action = this)
  gWidgets2::addHandler(this$parent, signal = "focus-out-event", handler = this$OnLostFocus, action = this)
  gWidgets2::addHandlerMouseMotion(this$plot_device, handler = this$OnMouseMotion)
  gWidgets2::addHandlerSelectionChanged(this$plot_device, handler = this$OnSelection, action = this)
  gWidgets2::addHandler(this$plot_device, signal = "size-allocate", handler = this$OnResize, action = this)

  ## Set the name for the class
  class(this) <- append(class(this),"SpectraPlotWidget")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
