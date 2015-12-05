#MALDI Image reconstruction by selected ion Top Windows
##########################################################

#Open img from Hdd directly, the easy way...
OpenMSI<-function()
{
  require(cairoDevice)
  require(gWidgets2)
  require(gWidgets2RGtk2)
  require(fields)

  fname<-gfile("Select an MSI file to open", type="open", multi = F, filter =  c("tar"="tar"), initial.dir = path.expand("~/"))
  if(length(fname) == 0)
  {
    return ()
  }

  this <- environment()
  bLoading<-T

  #Progress Bar windows close handler
  ProgressBar_Closed<-function( evt, ...)
  {
    if( bLoading )
    {
      gmessage("The MS image loading has been aborted by user. ", icon = "info", title = "MSI loading aborted")
      this$bLoading<-F
    }
  }

  #Create a progress bar
  pwin<-gwindow("", width = 400, height = 50, visible = F, parent=c(500, 500))
  addHandlerDestroy(pwin, handler = ProgressBar_Closed)
  grp_lay<-ggroup(horizontal = F, container = pwin)
  lbl<-glabel("Loading data please wait...", container = grp_lay)
  pbar<-gprogressbar(value = 0, container = grp_lay)
  visible(pwin)<-T
  Sys.sleep(0.1) #This forces a redraw on progressbar

  #Progress Bar update function
  setPbarValue<-function(progress)
  {
    svalue(this$pbar)<-progress
    Sys.sleep(0.1) #This forces a redraw on progressbar
    return(this$bLoading)
  }

  #Preloading 2 speedup
  raw<-LoadMsiData(data_file = fname, restore_path = file.path(dirname(fname), paste("ramdisk",basename(fname), sep = "_")), fun_progress_event = this$setPbarValue)

  if(is.null(raw))
  {
    #Process Aborted By User, return
    unlink( file.path(dirname(fname), paste("ramdisk",basename(fname), sep = "_")), recursive = T)
    return()
  }

  #Close the progressBar when data is loaded
  bLoading<-F
  dispose(pwin)

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

  #Pre copy just to respect paramenters naming
  in_img <-img
  rm(img)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  img <- in_img
  spectra_mass <- img$mass
  spectra_intensity <- img$mean@intensity
  rm(in_img) #Delete the parameter from the environment an keep only the copy
  imaging_dev <- 0
  spectraWidget <- 0
  Spin_Smooth <- 0
  Spin_Xmin <- 0
  Spin_Xmax <- 0
  Spin_Ymin <- 0
  Spin_Ymax <- 0
  Spin_rotation <- 0
  mz_tolerance <- 0
  mz_selected <- 0
  radio_buttons_Xres <- 0
  Tbl_spotList <- 0
  iSel <- 1:nrow(img$pos) #Initially grab the whole image
  image_range <- c( 1, 1, img$size["x"], img$size["y"] )

  #Methods

  #==================================================================================================
  ImgPlotFun <- function( mass, tol )
  {
    this$mz_tolerance <- tol
    this$mz_selected <- mass
    visible(this$imaging_dev)<-TRUE
    return(plotMassImageByPeak(this$img, mass.peak = this$mz_selected, tolerance = this$mz_tolerance, useColors = T, smoothFactor =   svalue(this$Spin_Smooth), XResLevel  = svalue(this$radio_buttons_Xres), selectedPixels = this$iSel, rotation = svalue(this$Spin_rotation)))
  }

  #==================================================================================================
  BtnClearSpotList <- function( mass, tol, ... )
  {
    this$spectraWidget$OnLostFocus() #This is just a test

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
    selected <- svalue(this$Tbl_spotList)
    this$spectraWidget$ClearSpectra()
    df<-data.frame(this$Tbl_spotList$get_items())  #get data frame...
    max_nrow<-nrow(this$img$data[[1]])

    for( i in selected)
    {
      selDf<-df[ df$ID == i,] #Get the correct data row
      if(i == 0)
      {
        this$spectraWidget$AddSpectra(this$spectra_mass, this$spectra_intensity, col = as.character(selDf$Colour))
      }
      else
      {
        icube<-(1+((i-1) %/% max_nrow))
        irow<- (i - (icube -1) * max_nrow)
        this$spectraWidget$AddSpectra(this$img$mass, this$img$data[[icube]][ irow ,] , col =  as.character(selDf$Colour))
      }
    }
  }

  #==================================================================================================
  BtnExportSpotList <- function( ... )
  {
    indexes<-as.vector(data.frame(this$Tbl_spotList$get_items())$ID)
    indexes<-indexes[indexes > 0] #Remove zero!
    store_paths<-gfile("Save current spots to txt files", type="selectdir", multi = F, initial.dir = path.expand("~/"))
    if(length(store_paths) == 0)
    {
      return ()
    }

    #display a BIG WaRNIng IF to much ID's are seected!
    bExportData <- T
    if( length(indexes) > 25)
    {
      bExportData<-gconfirm(paste("You are exporting a lot of data. (", length(indexes) ,"mass spectrums )\nThis may take a long time and expend a lot of memory.\nDo you want to store spectra as TXT?"), title = "Warning: Large data export!", icon = "warning")
    }

    .loadProgressBar("Exporting spectra to txt...")
    Sys.sleep(0.1) #This forces a redraw on progressbar
    .setPbarValue(0)

    #Create a dir to store all data inside
    store_paths <- file.path(store_paths, paste("Export_", tools::file_path_sans_ext(this$img$name),"_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), sep = "" ))
    dir.create( store_paths )


    #Save mean spectra
    spc <- matrix(data= c(this$img$mean@mass, this$img$mean@intensity), ncol = 2, byrow = F)
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
          .setPbarValue( 100* i/length(indexes))
        }
      }
    }

    .closeProgressBar() ####TODO aquesta funcio hauria de ser un metode de la classe
    rm(dataChunk)
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
    gmessage(smsg, title = "Export complete!", icon = "info")

  }

  #==================================================================================================
  OnPixelSelection <- function( evt, ...)
  {
    xResDiv<-switch(svalue(this$radio_buttons_Xres), x1 = 1, x2 = 2, x4 = 4)

    X_left<-round(min(evt$x)/xResDiv)
    X_right<-round(max(evt$x)/xResDiv)
    Y_bottom<-round(min(evt$y)/xResDiv)
    Y_top<-round(max(evt$y)/xResDiv)

    #Apply rotation
    rotate<-svalue(this$Spin_rotation)
    if(rotate >= 0 && rotate < 90)
    {
      Left <- X_left
      Right <- X_right
      Bottom <- Y_bottom
      Top <- Y_top
    }
    if(rotate >= 90 && rotate < 180)
    {
      Left <- Y_bottom
      Right <- Y_top
      Bottom <- this$img$size["y"] - X_right
      Top <- this$img$size["y"] - X_left
    }
    if(rotate >= 180 && rotate < 270)
    {
      Left <- this$img$size["x"] - X_right
      Right <- this$img$size["x"] - X_left
      Bottom <- this$img$size["y"] - Y_top
      Top <- this$img$size["y"] - Y_bottom
    }
    if(rotate >= 270 && rotate < 360)
    {
      Left <- this$img$size["x"] - Y_top
      Right <- this$img$size["x"] - Y_bottom
      Bottom <- X_left
      Top <- X_right
    }

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
        yimg <- 1 + this$img$size["y"] - yi
        preID <- which( Zpos == complex(real = xi, imaginary = yimg) )
        if(!(preID %in% currID))
        {
          X<-c(X, xi)
          Y<-c(Y, yi)
          Colour<- c(Colour ,hsv( h = preID/length(Zpos), s = 0.7, v = 1))
          ID<-c(ID, preID)
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

  #==================================================================================================
  SaveImg2Png <- function( ... )
  {
    fname<-gfile("Save current MSI plot to png file", type="save", multi = F, filter =  c("jpeg"="jpeg"), initial.dir = path.expand(getwd()))
    if(length(fname) == 0)
    {
      return ()
    }

    #Auto append the image file extension
    if(!grepl(".jpeg", basename(fname)))
    {
      fname<-paste(fname, ".jpeg", sep = "")
    }

    visible(this$imaging_dev)<-TRUE
    dev.print(jpeg, filename = fname, quality = "99", width = size(this$imaging_dev)["width"], height = size(this$imaging_dev)["height"])
  }

  #==================================================================================================
  SpinSmoothChanged <- function( ... )
  {
    visible(this$imaging_dev)<-TRUE
    plotMassImageByPeak(this$img, mass.peak = this$mz_selected, tolerance = this$mz_tolerance, useColors = T, smoothFactor =   svalue(this$Spin_Smooth), XResLevel  = svalue(this$radio_buttons_Xres), selectedPixels = this$iSel, rotation = svalue(this$Spin_rotation))
  }

  #==================================================================================================
  SpinImageRangeChanged <- function( ... )
  {
    #Grab new roi to plot
    X_left<-svalue(this$Spin_Xmin)
    X_right<-svalue(this$Spin_Xmax)
    Y_bottom<-svalue(this$Spin_Ymin)
    Y_top<-svalue(this$Spin_Ymax)

    #Apply rotation!
    rotate<-svalue(this$Spin_rotation)
    if(rotate >= 0 && rotate < 90)
    {
      x_min <- X_left
      x_max <- X_right
      y_min <- Y_bottom
      y_max <- Y_top
    }
    if(rotate >= 90 && rotate < 180)
    {
      x_min <- Y_bottom
      x_max <- Y_top
      y_min <- this$img$size["y"] - X_right
      y_max <- this$img$size["y"] - X_left
    }
    if(rotate >= 180 && rotate < 270)
    {
      x_min <- this$img$size["x"] - X_right
      x_max <- this$img$size["x"] - X_left
      y_min <- this$img$size["y"] - Y_top
      y_max <- this$img$size["y"] - Y_bottom
    }
    if(rotate >= 270 && rotate < 360)
    {
      x_min <- this$img$size["x"] - Y_top
      x_max <- this$img$size["x"] - Y_bottom
      y_min <- X_left
      y_max <- X_right
    }

    #Keep a copy of image limits
    this$image_range <- c( x_min, x_max, y_min, y_max )

    #Transform Y coord to fit img space
    y_min_aux<-y_min
    y_min<-1 + this$img$size["y"] -y_max
    y_max<-1 + this$img$size["y"] -y_min_aux

    iSel <- ((this$img$pos[, "x"] >= x_min) & (this$img$pos[, "x"] <= x_max))
    iSel <- iSel & ((this$img$pos[, "y"] >= y_min) & (this$img$pos[, "y"] <= y_max))
    this$iSel <- which(iSel)

    # DEBUG PRINTS
    #   cat("\nImage Range Changed!\n")
    #   cat(paste("x_min = ", x_min, "\n"))
    #   cat(paste("x_max = ", x_max, "\n"))
    #   cat(paste("y_min = ", y_min, "\n"))
    #   cat(paste("y_max = ", y_max, "\n"))
    #   cat("iSel:\n")
    #   print(this$.iSel)

    visible(this$imaging_dev)<-TRUE
    plotMassImageByPeak(this$img, mass.peak = this$mz_selected, tolerance = this$mz_tolerance, useColors = T, smoothFactor =   svalue(this$Spin_Smooth), XResLevel  = svalue(this$radio_buttons_Xres), selectedPixels = this$iSel, rotation = svalue(this$Spin_rotation))
  }

  #==================================================================================================
  SpinRotateChanged <- function( ... )
  {
    #Adjust range controls to fit rotation
    this$Spin_Xmin$remove_handlers()
    this$Spin_Xmax$remove_handlers()
    this$Spin_Ymin$remove_handlers()
    this$Spin_Ymax$remove_handlers()

    rotate<-svalue(this$Spin_rotation)
    if(rotate >= 0 && rotate < 90 || rotate >= 180 && rotate < 270)
    {
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmin), min = 1, max = this$img$size["x"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmax), min = 1, max = this$img$size["x"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymin), min = 1, max = this$img$size["y"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymax), min = 1, max = this$img$size["y"] )
    }
    if(rotate >= 90 && rotate < 180 || rotate >= 270 && rotate < 360)
    {
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmin), min = 1, max = this$img$size["y"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Xmax), min = 1, max = this$img$size["y"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymin), min = 1, max = this$img$size["x"] )
      gtkSpinButtonSetRange(getToolkitWidget(this$Spin_Ymax), min = 1, max = this$img$size["x"] )
    }

    #Apply rotation to image range
    if(rotate >= 0 && rotate < 90)
    {
      svalue(this$Spin_Xmin) <- this$image_range[1]
      svalue(this$Spin_Xmax) <- this$image_range[2]
      svalue(this$Spin_Ymin) <- this$image_range[3]
      svalue(this$Spin_Ymax) <- this$image_range[4]
    }
    if(rotate >= 90 && rotate < 180)
    {
      svalue(this$Spin_Xmin) <- this$img$size["y"] - this$image_range[4] + 1
      svalue(this$Spin_Xmax) <- this$img$size["y"] - this$image_range[3] + 1
      svalue(this$Spin_Ymin) <- this$image_range[1]
      svalue(this$Spin_Ymax) <- this$image_range[2]
    }
    if(rotate >= 180 && rotate < 270)
    {
      svalue(this$Spin_Xmin) <- this$img$size["x"] - this$image_range[2] + 1
      svalue(this$Spin_Xmax) <- this$img$size["x"] - this$image_range[1] + 1
      svalue(this$Spin_Ymin) <- this$img$size["y"] - this$image_range[4] + 1
      svalue(this$Spin_Ymax) <- this$img$size["y"] - this$image_range[3] + 1
    }
    if(rotate >= 270 && rotate < 360)
    {
      svalue(this$Spin_Xmin) <- this$image_range[3]
      svalue(this$Spin_Xmax) <- this$image_range[4]
      svalue(this$Spin_Ymin) <- this$img$size["x"] - this$image_range[2] + 1
      svalue(this$Spin_Ymax) <- this$img$size["x"] - this$image_range[1] + 1
    }

    #Re-connect handlers
    this$Spin_Xmin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Xmax$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymin$add_handler_changed(this$SpinImageRangeChanged)
    this$Spin_Ymax$add_handler_changed(this$SpinImageRangeChanged)

    #Plot rotated image
    visible(this$imaging_dev)<-TRUE
    plotMassImageByPeak(this$img, mass.peak = this$mz_selected, tolerance = this$mz_tolerance, useColors = T, smoothFactor =   svalue(this$Spin_Smooth), XResLevel  = svalue(this$radio_buttons_Xres), selectedPixels = this$iSel, rotation = svalue(this$Spin_rotation))
  }

  #==================================================================================================
  CheckBox_XRes_Changed <- function( ... )
  {
    visible(this$imaging_dev)<-TRUE
    plotMassImageByPeak(this$img, mass.peak = this$mz_selected, tolerance = this$mz_tolerance, useColors = T, smoothFactor =   svalue(this$Spin_Smooth), XResLevel  = svalue(this$radio_buttons_Xres), selectedPixels = this$iSel, rotation = svalue(this$Spin_rotation))
  }

  #GUI builder
  window <- gwindow ( paste("MSI Reconstruction -", img$name) , visible = F , width = 1024, height = 700)

  Grp_Top <- gpanedgroup(horizontal = F, container = window)
  Panel_Img<- gpanedgroup(horizontal = T, container = Grp_Top)

  ID<-0
  X<-0
  Y<-0
  Colour<-"red" #default colour for mean spectra
  Grp_Tbl <- ggroup(horizontal = F, container =Panel_Img)
  Tbl_spotList<-gtable( data.frame(ID,X,Y,Colour), container = Grp_Tbl, multiple = T, chosen.col = 1)

  ##Set table style using colors
  gtkTreeViewSetGridLines(getToolkitWidget(Tbl_spotList), 3)
  gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0),
                             gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0))[[1]],
                             background = 3
  )

  render<-gtkCellLayoutGetCells(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0))[[1]]
  render$set( font = "bold")
  gtkCellLayoutSetAttributes(gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 0), render)

  gtkTreeViewGetColumn(getToolkitWidget(Tbl_spotList), 3)$set(visible = F)

  Grp_BtmTbl <-ggroup(horizontal = T, container =Grp_Tbl)
  Btn_PlotSelSpotList<-gbutton("Plot", container= Grp_BtmTbl,  handler = this$SpectraListSelChange)
  Btn_ClearSpotList<-gbutton("Clear", container= Grp_BtmTbl,  handler = this$BtnClearSpotList)
  Btn_ExportSpotList<-gbutton("Export", container= Grp_BtmTbl,  handler = this$BtnExportSpotList)

  Grp_TopImg <- ggroup(horizontal = F, container = Panel_Img)
  Grp_Buttons <- ggroup(horizontal = T, container = Grp_TopImg)

  Lbl_Smooth<- glabel(text = "Smooth:", container = Grp_Buttons)
  Spin_Smooth<- gspinbutton(from = 0.1, to = 10, digest = 2, by = 0.2 , value = 0.3, handler = this$SpinSmoothChanged, container = Grp_Buttons)
  Lbl_Xres<- glabel(text = "Resolution:", container = Grp_Buttons)
  radio_buttons_Xres <- gradio( c("x1","x2", "x4"), selected = 2, horizontal = T, handler = this$CheckBox_XRes_Changed, action = this, container = Grp_Buttons)
  addSpring(Grp_Buttons)

  Lbl_Rotation<- glabel(text = "Rotation:", container = Grp_Buttons)
  Spin_rotation<- gspinbutton(from = 0, to = 270, digest = 0, by = 90 , value = 0, handler = this$SpinRotateChanged, container = Grp_Buttons)

  Lbl_XImgRange<- glabel(text = "X range:", container = Grp_Buttons)
  Spin_Xmin<- gspinbutton(from = 1, to =  img$size["x"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Spin_Xmax<- gspinbutton(from = 1, to = img$size["x"], digest = 0, by = 1 , value = img$size["x"], handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Lbl_YImgRange<- glabel(text = "Y range:", container = Grp_Buttons)
  Spin_Ymin<- gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = 1, handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Spin_Ymax<- gspinbutton(from = 1, to =  img$size["y"], digest = 0, by = 1 , value = img$size["y"], handler = this$SpinImageRangeChanged, container = Grp_Buttons)
  Btn_plot2file<- gbutton("Save to jpeg", container = Grp_Buttons, handler = this$SaveImg2Png)

  imagingFrame<-gframe("Ion Image", container = Grp_TopImg, fill = T, expand = T, spacing = 5)
  imaging_dev <- ggraphics(spacing = 5)
  addHandlerSelectionChanged( imaging_dev, handler = this$OnPixelSelection, action = this) ##TODO change this with the pointer to enviorment
  add(obj = imagingFrame, child = imaging_dev,  fill = T, expand = T)


  spectraFrame<-gframe("Average Spectra", container = Grp_Top,  fill = T, expand = T, spacing = 5)
  spectraWidget<-.SpectraPlotWidget(parent_widget = spectraFrame, clicFuntion = this$ImgPlotFun, showOpenFileButton = F)
  visible(window)<-TRUE

  spectraWidget$AddSpectra( spectra_mass, spectra_intensity, col = Colour)
  spectraWidget$ZoomResetClicked()


  window$widget$present() ##Tot i forzzar aki un raise-up del top widget en RStudio no deixa!

  ## Set the name for the class
  class(this) <- append(class(this),"MsiWindows")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  #return(this)
}
