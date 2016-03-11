#A GUI to select up to two MS images to load in rMSI data format

LoadTwoMsImages <- function( )
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  img1_obj<-NULL
  img2_obj<-NULL
  load_completed<-F

  #Signal handlers
  fileChooseClicked <- function(h, ...)
  {
    fname<-gWidgets2::gfile("Select an MSI file to open", type="open", multi = F, filter =  c("tar"="tar"), initial.dir = path.expand("~/"))
    if(length(fname) == 0)
    {
      return ()
    }

    gWidgets2::svalue(h$action) <- fname
  }

  cancelClicked <- function(h, ...)
  {
    this$pbar$ProgressBar_Closed()
    gWidgets2::dispose(h$obj)
  }

  Win_Closed <- function( h, ...)
  {
    if(!this$load_completed)
    {
      this$pbar$ProgressBar_Closed()
    }
  }

  okClicked <- function(h, ...)
  {
    #Hide and show some widgets
    RGtk2::gtkWidgetHide(gWidgets2::getToolkitWidget(this$btn_proced))
    gWidgets2::enabled(this$lbl_tl) <- F
    gWidgets2::enabled(this$frm_img1) <- F
    gWidgets2::enabled(this$frm_img2) <- F
    RGtk2::gtkWidgetShow(gWidgets2::getToolkitWidget(this$frm_pbar))

    #Check if images paths are correct and load data
    fpath1<- gWidgets2::svalue( this$txt_img1 )
    fpath2<- gWidgets2::svalue( this$txt_img2 )

    if( file.exists(fpath1) )
    {
      this$pbar$setTitle("Loading image 1")
      this$img1_obj<-LoadMsiData(data_file = fpath1, restore_path = file.path(dirname(fpath1), paste("ramdisk",basename(fpath1), sep = "_")), fun_progress = this$pbar$setValue)
      this$pbar$setValue(100)
    }
    else
    {
      if(fpath1 != "")
      {
        gWidgets2::gmessage("Image 1 file is not valid or does not exists", "Error loading image 1", icon = "error")
      }
    }


    if( file.exists(fpath2) )
    {
      this$pbar$setTitle("Loading image 2")
      this$img2_obj<-LoadMsiData(data_file = fpath2, restore_path = file.path(dirname(fpath2), paste("ramdisk",basename(fpath2), sep = "_")), fun_progress = this$pbar$setValue)
      this$pbar$setValue(100)
    }
    else
    {
      if(fpath2 != "")
    {
        gWidgets2::gmessage("Image 2 file is not valid or does not exists", "Error loading image 2", icon = "error")
      }
    }

    if(is.null(this$img1_obj))
    {
      #Process Aborted By User, return
      unlink( file.path(dirname(fpath1), paste("ramdisk",basename(fpath1), sep = "_")), recursive = T)
    }

    if(is.null(this$img2_obj))
    {
      #Process Aborted By User, return
      unlink( file.path(dirname(fpath2), paste("ramdisk",basename(fpath2), sep = "_")), recursive = T)
    }

    #Restore widgets show state
    if( gWidgets2::isExtant(this$win_tp ) )
    {
      RGtk2::gtkWidgetShow(gWidgets2::getToolkitWidget(this$btn_proced))
      gWidgets2::enabled(this$lbl_tl) <- T
      gWidgets2::enabled(this$frm_img1) <- T
      gWidgets2::enabled(this$frm_img2) <- T
      RGtk2::gtkWidgetHide(gWidgets2::getToolkitWidget(this$frm_pbar))
    }

    if( !file.exists(fpath1) && !file.exists(fpath2) )
    {
      gWidgets2::gmessage("At least one file must be selected.", "Error loading images", icon = "error")
    }
    else
    {
      this$load_completed <- T
      gWidgets2::dispose(h$obj)
    }
  }

  #Build GUI
  win_tp<-gWidgets2::gwindow("Load MSI data", visible = F, width = 600, height = 200, parent = c(0.5*RGtk2::gdkScreenWidth() - 300, 0.5*RGtk2::gdkScreenHeight() - 100))
  box_tp<-gWidgets2::ggroup(horizontal = F,container = win_tp, spacing = 10)
  box_title <- gWidgets2::ggroup(horizontal = T,container = box_tp, spacing = 10, expand = T, fill = T)
  lbl_tl<-gWidgets2::glabel("<span  weight=\"bold\">  Select up to two MS images to load in .tar format</span>", markup = T,  container = box_title)
  gWidgets2::addSpring(box_title)

  frm_img1<-gWidgets2::gframe(spacing = 5, container = box_tp)
  box_img1<-gWidgets2::ggroup(horizontal = T, container = frm_img1, spacing = 10, expand = T, fill = T)
  lbl_img1<-gWidgets2::glabel("  Image 1:", container = box_img1)
  txt_img1<-gWidgets2::gedit(width = 60, container = box_img1, expand = T, fill = T)
  btn_img1<-gWidgets2::gbutton("Select file", container = box_img1, handler = fileChooseClicked, action = txt_img1)

  frm_img2<-gWidgets2::gframe(spacing = 5, container = box_tp)
  box_img2<-gWidgets2::ggroup(horizontal = T, container = frm_img2, spacing = 10, expand = T, fill = T)
  lbl_img2<-gWidgets2::glabel("  Image 2:", container = box_img2)
  txt_img2<-gWidgets2::gedit(width = 60, container = box_img2, expand = T, fill = T)
  btn_img2<-gWidgets2::gbutton("Select file", container = box_img2, handler = fileChooseClicked, action = txt_img2)

  frm_pbar <- gWidgets2::gframe(spacing = 5, container = box_tp)
  pbar <- .ProgressBarDialog( "Loading data...",parent = frm_pbar)

  frm_btns <- gWidgets2::gframe( spacing = 5, container = box_tp)
  box_btns <- gWidgets2::ggroup(horizontal = T, spacing = 10, container = frm_btns, expand = T, fill = T)
  gWidgets2::addSpring(box_btns)
  btn_cancel <- gWidgets2::gbutton("Cancel", container = box_btns, handler = cancelClicked)
  btn_proced <- gWidgets2::gbutton("OK", container = box_btns, handler = okClicked)

  gWidgets2::addHandlerDestroy(win_tp, handler = Win_Closed)
  RGtk2::gtkWidgetHide(gWidgets2::getToolkitWidget(this$frm_pbar))

  gWidgets2::visible(win_tp)<-T
  ## Set the name for the class
  class(this) <- append(class(this),"LoadTwoMsImages")
  gc()

  #Do not return until this windos is disposed...
  while(gWidgets2::isExtant(this$win_tp ))
  {
    Sys.sleep(0.1)
  }

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
