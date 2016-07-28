###A GUI to display MS images from some ions/features

.MSImagePlotWidget <- function( in_img, parent_widget=gwindow ( "Default MSImagePlotWidget" , visible = FALSE ), AddSpectra_function = NULL, GetSpectraInfo_function = NULL, ClearSpectraPlot_function = NULL, meanSpectrumColor = "red", widget_name = "")
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  img <- in_img
  rm(in_img)
  parent <- parent_widget
  rm(parent_widget)
  mz_tolerance <- c()
  mz_selected <- c()
  imaging_dev <- 0
  Spin_Smooth <- 0
  Spin_Xmin <- 0
  Spin_Xmax <- 0
  Spin_Ymin <- 0
  Spin_Ymax <- 0
  Rotation <- 0
  Tbl_spotList <- 0
  Scale_light <- 0
  plotting_raster <- NULL #Current plotted MS image object
  ROI <- NULL #Current used ROI on image, NULL means no ROI
  ZOOM_win <- NULL #Current zoom windows
  IntLimit_ROI <- NULL #Intensity limiting ROI
  IntLimits <- NULL #A vector of intensity limits for each channel
  SidePanel_position <- 0 #Storing side panel position of spectra list in order to resotre it when it is hide/show
  AddSpectra_ptr <- AddSpectra_function #Pointer to a AddSpectra method of a spectraWidget to be able of ploting directly
  rm(AddSpectra_function)
  GetSpectraInfo_ptr <- GetSpectraInfo_function #Pointer to a GetSpectra method of a spectraWidget to be able to get current ploted spectra directly
  rm(GetSpectraInfo_function)
  ClearSpectraPlot_ptr <- ClearSpectraPlot_function #Pointer to a RmSpectra method of a spectraWidget to be able to clear ploted spectra directly
  rm(ClearSpectraPlot_function)
  myName <- widget_name
  NormalizationCoefs <- rep(1, nrow(img$pos))

  #Current image RGB layers
  Rlayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  Glayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )
  Blayer_raster <- .InitRGBEmptyRaster( img$size["x"], img$size["y"] )

  #==================================================================================================
  ImgBuildFun <- function(channel, mass, tol )
  {
    this$mz_tolerance[channel] <- tol
    this$mz_selected[channel] <- mass
    img_new<-.buildImageByPeak(this$img, mass.peak = mass, tolerance = tol, NormCoefs = this$NormalizationCoefs)

    #Apply intensity limitation directly to the raster object
    if( !is.null(this$IntLimit_ROI))
    {
      this$IntLimits[channel] <- max(raster::as.matrix(img_new$raster)[ (this$IntLimit_ROI[3]:this$IntLimit_ROI[4]), (this$IntLimit_ROI[1]:this$IntLimit_ROI[2])])
      raster::values(img_new$raster)[ raster::values(img_new$raster) > this$IntLimits[channel] ] <- this$IntLimits[channel]
    }
    else
    {
      this$IntLimits[channel] <- NA
    }

    mz_str <- this$mz_selected[channel]
    if(mz_str < 1000)
    {
      mz_str <- paste( round(mz_str, digits = 3), "Da" )
    }
    else
    {
      mz_str <- paste( round(mz_str/1000, digits = 3), "kDa" )
    }

    if( channel == 1)
    {
      this$Rlayer_raster<-img_new
      svalue(this$Lbl_RedMz) <- mz_str
    }
    else if (channel == 2)
    {
      this$Glayer_raster<-img_new
      svalue(this$Lbl_GreenMz) <- mz_str
    }
    else if (channel == 3)
    {
      this$Blayer_raster<-img_new
      svalue(this$Lbl_BlueMz) <- mz_str
    }

    #Return del buildImage
    return(list(selMz = img_new$mass, selTol = img_new$tolerance))
  }

  #==================================================================================================
  RedrawMSImage <-function()
  {
    visible(this$imaging_dev)<-TRUE
    if(is.null(this$ZOOM_win))
    {
      .plotMassImageRGB (this$plotting_raster, cal_um2pixels = this$img$pixel_size_um,  rotation = this$Rotation, display_axes = F,
                         roi_rectangle =  this$ROI, zoom = F)
    }
    else
    {
      .plotMassImageRGB (this$plotting_raster, cal_um2pixels = this$img$pixel_size_um,  rotation = this$Rotation, display_axes = F,
                         roi_rectangle =  this$ZOOM_win, zoom = T)
    }
  }

  #==================================================================================================
  PlotMassImageRGB <- function()
  {
    ch_count <- 0
    if( gWidgets2::svalue(this$Btn_RedEnable))
    {
      .setCheckBoxText(this$Btn_RedEnable, " ON ", background = "red", foreground = "white", font_size = "large", font_weight = "heavy")
      red_layer <- this$Rlayer_raster
      unique_layer <- red_layer
      ch_count <- ch_count + 1
    }
    else
    {
      .setCheckBoxText(this$Btn_RedEnable, " ON ", background = "darkred", foreground = "grey", font_size = "large", font_weight = "heavy")
      red_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }
    if( gWidgets2::svalue(this$Btn_GreenEnable))
    {
      .setCheckBoxText(this$Btn_GreenEnable, " ON ", background = "green", foreground = "white", font_size = "large", font_weight = "heavy")
      green_layer <- this$Glayer_raster
      unique_layer <- green_layer
      ch_count <- ch_count + 1
    }
    else
    {
      .setCheckBoxText(this$Btn_GreenEnable, " ON ", background = "darkgreen", foreground = "grey", font_size = "large", font_weight = "heavy")
      green_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }
    if( gWidgets2::svalue(this$Btn_BlueEnable))
    {
      .setCheckBoxText(this$Btn_BlueEnable, " ON ", background = "blue", foreground = "white", font_size = "large", font_weight = "heavy")
      blue_layer <- this$Blayer_raster
      unique_layer <- blue_layer
      ch_count <- ch_count + 1
    }
    else
    {
      .setCheckBoxText(this$Btn_BlueEnable, " ON ", background = "darkblue", foreground = "grey", font_size = "large", font_weight = "heavy")
      blue_layer <-.InitRGBEmptyRaster( this$img$size["x"], this$img$size["y"] )
    }

    if(ch_count < 1)
    {
      print("No selected data to plot image")
      return()
    }

    inter_level<-switch(svalue(this$Combo_Xres), x1 = 1, x2 = 2, x3 = 3, x4 = 4, x5 = 5)
    if(ch_count == 1)
    {
      this$plotting_raster<-.BuildSingleIonRGBImage( unique_layer,   XResLevel = inter_level, light =  svalue(this$Scale_light) )
    }
    else
    {
      this$plotting_raster<-.BuildRGBImage( imgR = red_layer, imgG = green_layer, imgB = blue_layer, XResLevel = inter_level, light =  svalue(this$Scale_light) )
    }

    this$RedrawMSImage()

    gWidgets2::visible(this$scaleRed_dev)<-TRUE
    if(ch_count == 1)
    {
      .plotIntensityScale(red_layer, light = svalue(this$Scale_light))
    }
    else
    {
      .plotIntensityScale(red_layer, "R", light = svalue(this$Scale_light) )
    }

    gWidgets2::visible(this$scaleGreen_dev)<-TRUE
    if(ch_count == 1)
    {
      .plotIntensityScale(green_layer, light = svalue(this$Scale_light))
    }
    else
    {
      .plotIntensityScale(green_layer, "G", light = svalue(this$Scale_light))
    }


    gWidgets2::visible(this$scaleBlue_dev)<-TRUE
    if(ch_count == 1)
    {
      .plotIntensityScale(blue_layer, light = svalue(this$Scale_light))
    }
    else
    {
      .plotIntensityScale(blue_layer, "B", light = svalue(this$Scale_light))
    }
  }

  #==================================================================================================
  BtnClearSpotList <- function( mass, tol, ... )
  {
    this$Tbl_spotList$set_items(data.frame(this$Tbl_spotList$get_items())[1,])
    gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0),
                               gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]],
                               background = 3 )

    render<-gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]]
    render$set( font = "bold")
    gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0), render)
    gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 3)$set(visible = F)
  }

  #==================================================================================================
  SpectraListSelChange <- function( ... )
  {
    selected <- gWidgets2::svalue(this$Tbl_spotList)
    df<-data.frame(this$Tbl_spotList$get_items())  #get data frame...
    color_list<-c()
    id_list<-c()

    for( i in selected)
    {
      selDf<-df[ df$ID == i,] #Get the correct data row
      color_list <- c(color_list, as.character(selDf$Colour))
      id_list <- c( id_list, i)
    }

    #Add spectra to plot
    plotSpectra( id_list, color_list)
  }

  #Plot spectra using parent plot function including normalization ==================================================
  plotSpectra <- function( id, colors)
  {
    intensity_list <- list()
    color_list <- list()
    id_list <- list()
    norm_vect <- c()

    id0_pos <- which(id == 0)
    if(length(id0_pos > 0))
    {
      #There is a mean spectrum
      if(class(this$img$mean) == "MassSpectrum")
      {
        #Addap to old data mean spectrum using MALDIquant object
        intensity_list[[length(intensity_list) + 1]] <- this$img$mean@intensity
      }
      else
      {
        intensity_list[[length(intensity_list) + 1]] <- this$img$mean
      }
      color_list[[length(color_list) + 1 ]] <- colors[id0_pos]
      norm_vect <- c(norm_vect, mean(this$NormalizationCoefs))
      id_list[[length(id_list) + 1]] <-  0
    }

    #Data
    id_fil_pos <- which(id > 0)
    if(length(id_fil_pos) > 0)
    {
      intensity_data <- loadImgCunckFromIds(this$img, id[id_fil_pos])
      for( i in 1:length(id_fil_pos))
      {
        intensity_list[[length(intensity_list) + 1]] <- intensity_data[i, ]
        color_list[[length(color_list) + 1 ]] <- colors[id_fil_pos[i]]
        norm_vect <- c(norm_vect, this$NormalizationCoefs[id[id_fil_pos[i]]])
        id_list[[length(id_list) + 1]] <-  id[id_fil_pos[i]]
      }
    }

    if(!is.null(this$AddSpectra_ptr))
    {
      this$AddSpectra_ptr(this$img$mass, intensity_list, color_list, id_list, this$myName, norm_vect)
    }
  }

  #==================================================================================================
  BtnExportSpotList <- function( ... )
  {
    indexes<-as.vector(data.frame(this$Tbl_spotList$get_items())$ID)
    indexes<-indexes[indexes > 0] #Remove zero!
    store_paths<-gWidgets2::gfile("Save current spots to txt files", type="selectdir", multi = F, initial.dir = path.expand("~/"))
    if(length(store_paths) == 0)
    {
      return ()
    }

    #display a BIG WaRNIng IF to much ID's are seected!
    bExportData <- T
    if( length(indexes) > 25)
    {
      bExportData<-gWidgets2::gconfirm(paste("You are exporting a lot of data. (", length(indexes) ,"mass spectrums )\nThis may take a long time and expend a lot of memory.\nDo you want to store spectra as TXT?"), title = "Warning: Large data export!", icon = "warning")
    }

    mPbar<-.ProgressBarDialog("Exporting spectra to txt...")

    #Create a dir to store all data inside
    store_paths <- file.path(store_paths, paste("Export_", tools::file_path_sans_ext(this$img$name),"_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), sep = "" ))
    dir.create( store_paths )


    #Save mean spectra
    if(class(this$img$mean) == "MassSpectrum")
    {
      #Handling old mean based on MALDIquant
      spc <- matrix(data= c(this$img$mean@mass, this$img$mean@intensity), ncol = 2, byrow = F)
    }
    else
    {
      spc <- matrix(data= c(this$img$mass, this$img$mean), ncol = 2, byrow = F)
    }
    write( x = t(spc), file = file.path( store_paths, "average.txt" ), ncolumns = 2  )


    if(length(indexes) > 0)
    {
      #Save the id file
      id_fname <- file.path( store_paths, "ID.txt")
      write(indexes, file = id_fname, ncolumns = 1)

      if(bExportData)
      {
        #Save each ID in the lit
        dataChunck <- loadImgCunckFromIds(this$img, indexes)

        for( i in 1:length(indexes))
        {
          spc <- matrix(data= c(this$img$mass, dataChunck[i, ]), ncol = 2, byrow = F)
          spc_fname<-  file.path( store_paths, paste("ID_", indexes[i] ,".txt", sep =""))
          write( x = t(spc), file = spc_fname, ncolumns = 2  )
          if( !mPbar$setValue( 100* i/length(indexes) ))
          {
            #Aborted by user
            rm(dataChunck)
            unlink( store_paths, recursive = T )
            gc()
            return()
          }
        }
      }
    }

    mPbar$close()

    rm(dataChunck)
    gc()
    #Show a message of export completed
    smsg<- paste("Export complete!\nYour data is stored at:\n\t", store_paths , sep = "")
    if(bExportData)
    {
      if(length(indexes) > 0)
      {
        smsg <- paste(smsg,"\nAverage spectrum and", length(indexes), "spectrums exported.")
      }
      else
      {
        smsg <- paste(smsg,"\nOnly average spetrum exported.")
      }

    }
    else
    {
      smsg <- paste(smsg, "\nLarge data export disabled. Ony average spectrum and ID list exported.")
    }
    gWidgets2::gmessage(smsg, title = "Export complete!", icon = "info")

  }

  #==================================================================================================
  OnPixelSelection <- function( evt, ...)
  {
    X_left<-round(min(evt$x))
    X_right<-round(max(evt$x))
    Y_bottom<-round(min(evt$y)) #Transform raster coords to image coords (only Y axis is affected)
    Y_top<-round(max(evt$y)) #Transform raster coords to image coords (only Y axis is affected)

    #Apply rotation!
    if(this$Rotation == 0)
    {
      this$ROI <- c(X_left + 1, X_right, this$img$size["y"] - Y_top + 1, this$img$size["y"] - Y_bottom )
    }
    if(this$Rotation == 90)
    {
      this$ROI <- c(Y_bottom + 1, Y_top, X_left + 1, X_right)
    }
    if(this$Rotation == 180)
    {
      this$ROI <- c( this$img$size["x"] - X_right + 1, this$img$size["x"] - X_left, Y_bottom + 1, Y_top)
    }
    if(this$Rotation == 270)
    {
      this$ROI <- c( this$img$size["x"] - Y_top + 1, this$img$size["x"] - Y_bottom, this$img$size["y"] - X_right + 1, this$img$size["y"] - X_left)
    }

    #Set it to ROI spinbuttons
    this$Spin_Xmin$remove_handlers()
    this$Spin_Xmax$remove_handlers()
    this$Spin_Ymin$remove_handlers()
    this$Spin_Ymax$remove_handlers()
    svalue(this$Spin_Xmin) <- this$ROI[1]
    svalue(this$Spin_Xmax) <- this$ROI[2]
    svalue(this$Spin_Ymin) <- this$ROI[3]
    svalue(this$Spin_Ymax) <- this$ROI[4]
    this$Spin_Xmin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Xmax$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymax$add_handler_changed(this$SpinImageRangeChanged)

    this$RedrawMSImage()
    gWidgets2::enabled(this$Btn_RoiZoom) <- T
    gWidgets2::enabled(this$Frame_RoiCtl) <- T
  }

  #==================================================================================================
  SaveImg2Png <- function( ... )
  {
    mass_sel <- c()
    tol_sel <-c ()
    if( svalue(this$Btn_RedEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[1])
      tol_sel<-c(tol_sel, mz_tolerance[1])
    }
    if( svalue(this$Btn_GreenEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[2])
      tol_sel<-c(tol_sel, mz_tolerance[2])
    }
    if( svalue(this$Btn_BlueEnable))
    {
      mass_sel<-c(mass_sel, mz_selected[3])
      tol_sel<-c(tol_sel, mz_tolerance[3])
    }
    if(length(mass_sel) == 0)
    {
      gWidgets2::gmessage("No channel enabled, nothing is exported.", icon = "info")
      return()
    }


    fname<-gWidgets2::gfile("Save current MSI plot to png file", type="save", multi = F, filter =  c("tiff"="tiff", "svg"="svg"), initial.dir = path.expand(getwd()))
    if(length(fname) == 0)
    {
      return ()
    }

    #Auto append the image file extension
    if(grepl(".svg", basename(fname)))
    {
      fname<-paste(fname, ".svg", sep = "")
      svg( filename = fname , width = 1200, height = 500)
    }
    else
    {
      if(!grepl(".tiff", basename(fname)))
      {
        fname<-paste(fname, ".tiff", sep = "")
      }
      tiff( filename = fname , width = 1200, height = 500, compression = "none", res = 160)
    }

    plotMassImageByPeak(this$img,  mass.peak = mass_sel, tolerance = tol_sel,
                        XResLevel = switch(gWidgets2::svalue(this$Combo_Xres), x1 = 1, x2 = 2, x3 = 3, x4 = 4, x5 = 5),
                        rotation = this$Rotation, vlight= gWidgets2::svalue(this$Scale_light),
                        crop_area = ZOOM_win, intensity_limit = this$IntLimits, NormalizationCoefs = this$NormalizationCoefs)
    dev.off()

    gWidgets2::gmessage(paste("Image saved at:", fname), icon = "info")
  }

  #==================================================================================================
  SpinImageRangeChanged <- function( ... )
  {
    #Set ROI from spinbuttons
    this$ROI[1]<- svalue(this$Spin_Xmin)
    this$ROI[2] <- svalue(this$Spin_Xmax)
    this$ROI[3] <- svalue(this$Spin_Ymin)
    this$ROI[4] <- svalue(this$Spin_Ymax)
    this$RedrawMSImage()
  }

  #==================================================================================================
  SliderLightChanged<- function( ... )
  {
    #Set it to ROI spinbuttons
    this$PlotMassImageRGB()
  }

  #==================================================================================================
  BtnRotateCCW <- function ( ... )
  {
    this$Rotation <- this$Rotation + 90
    if(this$Rotation == 360)
    {
      this$Rotation <- 0
    }
    this$RotateImage(this$Rotation)
  }

  #==================================================================================================
  BtnRotateCW <- function ( ... )
  {
    this$Rotation <- this$Rotation - 90
    if(this$Rotation == -90)
    {
      this$Rotation <- 270
    }
    this$RotateImage(this$Rotation)
  }

  #==================================================================================================
  RotateImage <- function( angle )
  {
    rotateLabel <- this$Rotation
    if( rotateLabel == 90 ) {rotateLabel<-270}
    else if(rotateLabel == 270){ rotateLabel<-90}
    Lbl_Rotation$set_value(paste("Rotation:", rotateLabel))

    #Plot rotated image
    this$RedrawMSImage()
  }

  #==================================================================================================
  ComboBox_XRes_Changed <- function( ... )
  {
    this$PlotMassImageRGB()
  }

  #==================================================================================================
  ComboBox_Norm_Changed <- function( ... )
  {
    if(gWidgets2::svalue(this$Combo_Norm) == "RAW")
    {
      this$NormalizationCoefs <- rep(1, nrow(this$img$pos))
    }
    else
    {
      this$NormalizationCoefs <- this$img$normalization[[gWidgets2::svalue(this$Combo_Norm)]]
    }

    #Re-Build the raster with the new Normalization factor
    for( ich in 1:length(this$mz_selected))
    {
      this$ImgBuildFun(ich, this$mz_selected[ich], this$mz_tolerance[ich] )
    }
    this$PlotMassImageRGB()

      #Re-Plot spectra using norm factors
    if(!is.null(this$GetSpectraInfo_ptr))
    {
      plotted_data <- this$GetSpectraInfo_ptr(this$myName)

      if(!is.null(this$ClearSpectraPlot_ptr))
      {
        this$ClearSpectraPlot_ptr(plotted_data$ID, this$myName)
        this$plotSpectra(plotted_data$ID, plotted_data$color)
      }
    }
  }

  #==================================================================================================
  ROI_Deleted <-function (...)
  {
    #Set it to ROI spinbuttons
    this$Spin_Xmin$remove_handlers()
    this$Spin_Xmax$remove_handlers()
    this$Spin_Ymin$remove_handlers()
    this$Spin_Ymax$remove_handlers()
    svalue(this$Spin_Xmin) <- 1
    svalue(this$Spin_Xmax) <- img$size["x"]
    svalue(this$Spin_Ymin) <- 1
    svalue(this$Spin_Ymax) <- img$size["y"]
    this$Spin_Xmin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Xmax$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymax$add_handler_changed(this$SpinImageRangeChanged)

    this$ROI <- NULL
    this$RedrawMSImage()
    gWidgets2::enabled(this$Frame_RoiCtl) <- F
    if(is.null(this$ZOOM_win))
    {
      gWidgets2::enabled(this$Btn_RoiZoom) <- F
    }
  }

  #==================================================================================================
  ROI_Zoom <- function( ... )
  {
    this$ZOOM_win <- switch(svalue(this$Btn_RoiZoom) , this$ROI)
    this$RedrawMSImage()

    if( is.null(this$ZOOM_win) )
    {
      this$Btn_RoiZoom[] <-"Zoom in ROI"
    }
    else
    {
      this$Btn_RoiZoom[] <-"Zoom Out"
    }

    #No roi defined and zoom out
    if(is.null(this$ZOOM_win) && is.null(this$ROI))
    {
      gWidgets2::enabled(this$Btn_RoiZoom) <- F
    }
  }

  #==================================================================================================
  ROI_GetSpectra <- function( ... )
  {
    if(!is.null(this$ROI))
    {
      Left <- this$ROI[1]
      Right <- this$ROI[2]
      Bottom <- this$ROI[3]
      Top <- this$ROI[4]

      #Limits to image size
      Left<-max(1, Left)
      Right<-max(1, Right)
      Bottom<-max(1, Bottom)
      Top<-max(1, Top)

      Left<-min(this$img$size["x"], Left)
      Right<-min(this$img$size["x"], Right)
      Bottom<-min(this$img$size["y"], Bottom)
      Top<-min(this$img$size["y"], Top)

      ID<-c()
      X<-c()
      Y<-c()
      Colour<-c()
      Zpos <-complex(real = this$img$pos[,"x"], imaginary = this$img$pos[,"y"]) #Convert positions to a complex numbers vector to find pointed coords fast and easy
      currID <- data.frame(this$Tbl_spotList$get_items())$ID
      for(xi in Left:Right)
      {
        for(yi in Bottom:Top)
        {
          preID <- which( Zpos == complex(real = xi, imaginary = yi) )

          if(length(preID) > 0)
          {
            if(!(preID %in% currID))
            {
              X<-c(X, xi)
              Y<-c(Y, yi)
              Colour<- c(Colour ,hsv( h = preID/length(Zpos), s = 0.7, v = 1))
              ID<-c(ID, preID)
            }
          }
        }
      }
      this$Tbl_spotList$set_items(rbind(data.frame(this$Tbl_spotList$get_items()), data.frame(ID, X, Y, Colour)))
      gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0),
                                 gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]],
                                 background = 3
      )

      render<-gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0))[[1]]
      render$set( font = "bold")
      gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 0), render)

      gtkTreeViewGetColumn(getToolkitWidget(this$Tbl_spotList), 3)$set(visible = F)
    }
  }

  #==================================================================================================
  IntensityScale_EnableClicked <- function( evt, ...)
  {
    this$PlotMassImageRGB()
  }

  #==================================================================================================
  ROI_IntensityLimit <- function( ... )
  {
    if( is.null(this$ROI))
    {
      gWidgets2::gmessage("To apply intensity limitation you must define a ROI ", icon = "info")
      return()
    }
    else
    {
      this$IntLimit_ROI <- this$ROI
      #Limit Roi to image range
      this$IntLimit_ROI[1] <- max( c(this$IntLimit_ROI[1], 1 ) )
      this$IntLimit_ROI[2] <- min( c(this$IntLimit_ROI[2], this$img$size["x"] ) )
      this$IntLimit_ROI[3] <- max( c(this$IntLimit_ROI[3], 1 ) )
      this$IntLimit_ROI[4] <- min( c(this$IntLimit_ROI[4], this$img$size["y"] ) )
    }

    #Re-Build the raster with the new intensity limit
    for( ich in 1:length(this$mz_selected))
    {
      this$ImgBuildFun(ich, this$mz_selected[ich], this$mz_tolerance[ich] )
    }
    this$PlotMassImageRGB()

    gWidgets2::enabled(this$Btn_RoiIntUnLimit) <- T
    gWidgets2::svalue(this$Btn_RoiIntUnLimit)<- "Remove Intensity Limit"
    this$Btn_RoiIntUnLimit$add_handler_changed(this$ROI_IntensityUnLimit)
  }

  #==================================================================================================
  ROI_IntensityUnLimit <- function( ... )
  {
    #Disable intensity limitation
    this$IntLimit_ROI <- NULL

    #Re-Build the raster with the new intensity limit
    for( ich in 1:length(this$mz_selected))
    {
      this$ImgBuildFun(ich, this$mz_selected[ich], this$mz_tolerance[ich] )
    }
    this$PlotMassImageRGB()

    this$Btn_RoiIntUnLimit$remove_handlers()
    gWidgets2::enabled(this$Btn_RoiIntUnLimit) <- F
    gWidgets2::svalue(this$Btn_RoiIntUnLimit)<- "No Intensity Limit"
  }

  #===================================================================================================
  ShowSpectraList <- function ( ... )
  {
    gWidgets2::add(this$Grp_SidePanel, this$spectraListFrame, expand=TRUE, fill = TRUE)
    gWidgets2::svalue(this$Panel_Img) <- this$SidePanel_position
    RGtk2::gtkWidgetHide(gWidgets2::getToolkitWidget(this$Btn_Show))
  }

  #===================================================================================================
  HideSpectraList <- function ( ... )
  {
    gWidgets2::delete(this$Grp_SidePanel, this$spectraListFrame)
    this$SidePanel_position <- gWidgets2::svalue(this$Panel_Img)
    gWidgets2::svalue(this$Panel_Img) <- 0
    RGtk2::gtkWidgetShow(gWidgets2::getToolkitWidget(this$Btn_Show))
  }

  #Build the GUI
  Top_frm <- gWidgets2::gframe( text =  paste("Image:",img$name), container = parent, spacing = 5 )
  Panel_Img<- gWidgets2::gpanedgroup(horizontal = T, container = Top_frm,  expand=TRUE )
  Grp_SidePanel<- gWidgets2::ggroup(container = Panel_Img, expand=TRUE, fill = TRUE)
  spectraListFrame<-gWidgets2::gframe("Spectra List", container = Grp_SidePanel,  fill = T, spacing = 5, expand = T )
  Grp_Tbl <- gWidgets2::ggroup(horizontal = F, container = spectraListFrame,  expand=TRUE, fill = TRUE)
  Btn_Hide <- gWidgets2::gbutton("<<< Hide", container = Grp_Tbl, handler = this$HideSpectraList )
  ID<-0
  X<-0
  Y<-0
  Tbl_spotList<-gWidgets2::gtable( data.frame(ID,X,Y, Colour = meanSpectrumColor), container = Grp_Tbl, multiple = T, chosen.col = 1)
  size( Tbl_spotList )<- c(120, -1)
  ##Set table style using colors
  RGtk2::gtkTreeViewSetGridLines(getToolkitWidget(Tbl_spotList), 3)
  RGtk2::gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0),
                             gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0))[[1]],
                             background = 3
                              )

  render<-RGtk2::gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0))[[1]]
  render$set( font = "bold")
  gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0), render)

  gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 3)$set(visible = F)

  Grp_BtmTbl <-gWidgets2::ggroup(horizontal = F, container =Grp_Tbl)
  Btn_PlotSelSpotList<-gWidgets2::gbutton("Plot", container= Grp_BtmTbl,  handler = this$SpectraListSelChange)
  Btn_ClearSpotList<-gWidgets2::gbutton("Clear", container= Grp_BtmTbl,  handler = this$BtnClearSpotList)
  Btn_ExportSpotList<-gWidgets2::gbutton("Export", container= Grp_BtmTbl,  handler = this$BtnExportSpotList)

  Grp_TopImg <- gWidgets2::ggroup(horizontal = F, container = Panel_Img, expand = T)
  Grp_Buttons <- gWidgets2::ggroup(horizontal = T, container = Grp_TopImg)
  Btn_Show <- gWidgets2::gbutton(">>>Spectra List", container = Grp_Buttons, handler = this$ShowSpectraList )
  RGtk2::gtkWidgetHide(gWidgets2::getToolkitWidget(this$Btn_Show))
  Lbl_Rotation<- gWidgets2::glabel(text = "Rotation: 0", container = Grp_Buttons)
  Btn_rotate_CCW <- gWidgets2::gbutton("", container = Grp_Buttons, handler = this$BtnRotateCCW)
  RGtk2::gtkImageSetFromFile( getToolkitWidget(Btn_rotate_CCW)$image, filename = file.path(system.file(package = "rMSI", "icons"),"Rotate_CCW.png") )
  Btn_rotate_CW <- gWidgets2::gbutton("", container = Grp_Buttons, handler = this$BtnRotateCW)
  RGtk2::gtkImageSetFromFile( getToolkitWidget(Btn_rotate_CW)$image, filename = file.path(system.file(package = "rMSI", "icons"),"Rotate_CW.png") )
  Lbl_Xres<- gWidgets2::glabel(text = "Interpolation:", container = Grp_Buttons)
  Combo_Xres <- gWidgets2::gcombobox( items = c("x1","x2","x3","x4","x5"), selected = 2, container = Grp_Buttons, handler = this$ComboBox_XRes_Changed)
  Lbl_Normalitzation <- gWidgets2::glabel(text = "Normalization:", container = Grp_Buttons)
  Combo_Norm <- gWidgets2::gcombobox( items = c("RAW", names(img$normalizations)), selected = 1, container = Grp_Buttons, handler = this$ComboBox_Norm_Changed)
  gWidgets2::size(Combo_Norm) <- c(100, -1)
  gWidgets2::glabel("Light:", container = Grp_Buttons)
  Scale_light <- gWidgets2::gslider( from = 0.6, to = 10, by = 0.2, value = 3, horizontal = T, handler = this$SliderLightChanged, container =  Grp_Buttons)
  gWidgets2::addSpring(Grp_Buttons)
  Btn_plot2file<- gWidgets2::gbutton("Save in image file", container = Grp_Buttons, handler = this$SaveImg2Png)

  Grp_ImgTop<-gWidgets2::ggroup( horizontal = T, container =  Grp_TopImg,  fill = T, expand = T)
  Grp_ImgRoi<-gWidgets2::ggroup( horizontal = F, container =  Grp_ImgTop,  fill = T, expand = T)
  imaging_dev <- gWidgets2::ggraphics(spacing = 5 )
  size( imaging_dev )<- c(400, 280)
  gWidgets2::addHandlerSelectionChanged( imaging_dev, handler = this$OnPixelSelection, action = this)
  gWidgets2::add(obj = Grp_ImgRoi, child = imaging_dev,  fill = T, expand = T)

  Grp_ScalesV <- gWidgets2::ggroup( horizontal = F, container =  Grp_ImgTop,  fill = T, expand = F)
  Grp_ScalesH <- gWidgets2::ggroup( horizontal = T, container =  Grp_ScalesV,  fill = T, expand = T)
  Grp_ZoomLimit<-gWidgets2::ggroup(horizontal = T, container = Grp_ScalesV)
  Btn_RoiZoom<-gWidgets2::gcheckbox("Zoom in ROI", checked = F, use.togglebutton = T, container = Grp_ZoomLimit, handler = this$ROI_Zoom)
  Btn_RoiIntUnLimit<-gWidgets2::gbutton("No Intensity Limit", container = Grp_ZoomLimit)

  #Red Color Scale
  Grp_RedScale<-gWidgets2::ggroup( horizontal = F, container = Grp_ScalesH)
  scaleRed_dev <- gWidgets2::ggraphics(spacing = 5 )
  size( scaleRed_dev )<- c(120, 250)
  gWidgets2::add(obj = Grp_RedScale, child = scaleRed_dev,  fill = T, expand = T)
  Grp_RedCtl <- gWidgets2::ggroup( horizontal = T, container = Grp_RedScale)
  Btn_RedEnable<-gWidgets2::gcheckbox("On", container = Grp_RedCtl, use.togglebutton = T, checked = T,  handler = this$IntensityScale_EnableClicked, action = "R")
  Lbl_RedMz <- gWidgets2::glabel("", container = Grp_RedCtl)

  #Green Color scale
  Grp_GreenScale<-gWidgets2::ggroup( horizontal = F, container = Grp_ScalesH)
  scaleGreen_dev <- gWidgets2::ggraphics(spacing = 5 )
  size( scaleGreen_dev )<- c(120, 250)
  gWidgets2::add(obj = Grp_GreenScale, child = scaleGreen_dev,  fill = T, expand = T)
  Grp_GreenCtl <- gWidgets2::ggroup( horizontal = T, container = Grp_GreenScale)
  Btn_GreenEnable<-gWidgets2::gcheckbox("On", container = Grp_GreenCtl, use.togglebutton = T, checked = F, handler = this$IntensityScale_EnableClicked, action = "G")
  Lbl_GreenMz <- gWidgets2::glabel("", container = Grp_GreenCtl)

  #Blue Color scale
  Grp_BlueScale<-gWidgets2::ggroup( horizontal = F, container = Grp_ScalesH)
  scaleBlue_dev <- gWidgets2::ggraphics(spacing = 5 )
  size( scaleBlue_dev )<- c(120, 250)
  gWidgets2::add(obj = Grp_BlueScale, child = scaleBlue_dev,  fill = T, expand = T)
  Grp_BlueCtl <- gWidgets2::ggroup( horizontal = T, container = Grp_BlueScale)
  Btn_BlueEnable<-gWidgets2::gcheckbox("On", container = Grp_BlueCtl, use.togglebutton = T, checked = F, handler = this$IntensityScale_EnableClicked, action = "B")
  Lbl_BlueMz <- gWidgets2::glabel("", container = Grp_BlueCtl)

  #Modify Red, Green, Blue buttons labels to add colors
  .setCheckBoxText(Btn_RedEnable, " ON ", background = "red", foreground = "white", font_size = "large", font_weight = "heavy")
  .setCheckBoxText(Btn_GreenEnable, " ON ", background = "darkgreen", foreground = "grey", font_size = "large", font_weight = "heavy")
  .setCheckBoxText(Btn_BlueEnable, " ON ", background = "darkblue", foreground = "grey", font_size = "large", font_weight = "heavy")

  #ROI CTL
  Frame_RoiCtl<-gWidgets2::gframe("ROI Controls", container = Grp_ImgRoi )
  Grp_RoiCtl<-gWidgets2::ggroup(horizontal = F, container = Frame_RoiCtl)
  Grp_RoiCtl_1stRow<-gWidgets2::ggroup(horizontal = T, container = Grp_RoiCtl)
  Grp_RoiCtl_2ndRow<-gWidgets2::ggroup(horizontal = T, container = Grp_RoiCtl)
  Btn_RoiDelete<-gWidgets2::gbutton("Delete", container = Grp_RoiCtl_1stRow, handler = this$ROI_Deleted)
  Btn_RoiGetSpectra<-gWidgets2::gbutton("Get Spectra", container = Grp_RoiCtl_1stRow, handler = this$ROI_GetSpectra)
  gWidgets2::addSpring(Grp_RoiCtl_1stRow)
  Lbl_XImgRange<- gWidgets2::glabel(text = "X range:", container = Grp_RoiCtl_1stRow)
  Spin_Xmin<- gWidgets2::gspinbutton(from = 1, to =  img$size["x"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_1stRow)
  Spin_Xmax<- gWidgets2::gspinbutton(from = 1, to = img$size["x"], digest = 0, by = 1 , value = img$size["x"], handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_1stRow)
  Btn_RoiIntLimit<-gWidgets2::gbutton("Apply Intensity Limit", container = Grp_RoiCtl_2ndRow, handler = this$ROI_IntensityLimit)
  gWidgets2::addSpring(Grp_RoiCtl_2ndRow)
  Lbl_YImgRange<- gWidgets2::glabel(text = "Y range:", container = Grp_RoiCtl_2ndRow)
  Spin_Ymin<- gWidgets2::gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_2ndRow)
  Spin_Ymax<- gWidgets2::gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = img$size["y"], handler = this$SpinImageRangeChanged, container = Grp_RoiCtl_2ndRow)

  gWidgets2::enabled(Btn_RoiZoom) <- F
  gWidgets2::enabled(Frame_RoiCtl) <- F
  gWidgets2::enabled(Btn_RoiIntUnLimit) <- F

  ## Set the name for the class
  class(this) <- append(class(this),"MSImagePlotWidget")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
