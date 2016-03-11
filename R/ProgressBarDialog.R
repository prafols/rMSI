###A GUI to display a progress bar dialog

.ProgressBarDialog <- function( Title="", parent = NULL )
{
  options(guiToolkit="RGtk2") # Força que toolquit sigu GTK pq fas crides directes a events GTK!!!
  oldWarning<-options()$warn
  options(warn = -1)

  ## Get the environment for this
  ## instance of the function.
  this <- environment()

  ##Class data members
  bKeepLoading<-T
  sTitle <- Title
  rm(Title)

  #Progress Bar update method
  setValue<-function(progress)
  {
    gWidgets2::svalue(this$pbar)<-progress
    this$lbl$set_value(paste( sTitle, "-", round(progress, digits = 2),"%") )
    Sys.sleep(0.1) #This forces a redraw on progressbar
    return(this$bKeepLoading)
  }

  #Set title
  setTitle<-function(Title)
  {
    this$sTitle <- Title
    this$lbl$set_value(paste( sTitle, "-", round(gWidgets2::svalue(this$pbar), digits = 2),"%") )
    Sys.sleep(0.1) #This forces a redraw on progressbar
  }

  #Close method
  close<-function()
  {
    if( this$bKeepLoading )
    {
      this$bKeepLoading<-F
      gWidgets2::dispose(this$pwin)
    }
  }

  #Close event handler
  ProgressBar_Closed<-function( evt, ...)
  {
    if( this$bKeepLoading )
    {
      gWidgets2::gmessage("The loading process has been aborted by the user. ", icon = "info", title = "User aborted")
      this$bKeepLoading<-F
    }
  }

  #Build GUI
  if(is.null(parent))
  {
    pwin<-gWidgets2::gwindow("Working...", width = 400, height = 50, visible = F, parent=c(500, 500))
    gWidgets2::addHandlerDestroy(pwin, handler = ProgressBar_Closed)
  }
  else
  {
    pwin <- parent
  }
  grp_lay<-gWidgets2::ggroup(horizontal = F, container = pwin, expand = T, fill = T)
  lbl<-gWidgets2::glabel(paste( sTitle, "- 0.00 %") , container = grp_lay)
  pbar<-gWidgets2::gprogressbar(value = 0, container = grp_lay, expand = T, fill = T)
  if(is.null(parent))
  {
    gWidgets2::visible(pwin)<-T
  }
  Sys.sleep(0.1) #This forces a redraw on progressbar

  ## Set the name for the class
  class(this) <- append(class(this),"ProgressBarDialog")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
