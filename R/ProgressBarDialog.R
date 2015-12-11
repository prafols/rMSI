###A GUI to display a progress bar dialog

.ProgressBarDialog <- function( Title )
{
  require(gWidgets2)
  require(gWidgets2RGtk2)

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
    svalue(this$pbar)<-progress
    this$lbl$set_value(paste( sTitle, "-", round(progress, digits = 2),"%") )
    Sys.sleep(0.1) #This forces a redraw on progressbar
    return(this$bKeepLoading)
  }

  #Close method
  close<-function()
  {
    if( this$bKeepLoading )
    {
      this$bKeepLoading<-F
      dispose(this$pwin)
    }
  }

  #Close event handler
  ProgressBar_Closed<-function( evt, ...)
  {
    if( this$bKeepLoading )
    {
      gmessage("The prosess has been aborted by the user. ", icon = "info", title = "User aborted")
      this$bKeepLoading<-F
    }
  }

  #Build GUI
  pwin<-gwindow("Working...", width = 400, height = 50, visible = F, parent=c(500, 500))
  addHandlerDestroy(pwin, handler = ProgressBar_Closed)
  grp_lay<-ggroup(horizontal = F, container = pwin)
  lbl<-glabel(paste( sTitle, "- 0.00 %") , container = grp_lay)
  pbar<-gprogressbar(value = 0, container = grp_lay)
  visible(pwin)<-T
  Sys.sleep(0.1) #This forces a redraw on progressbar

  ## Set the name for the class
  class(this) <- append(class(this),"ProgressBarDialog")
  gc()

  #Restore warnings level
  options(warn = oldWarning)
  rm(oldWarning)

  return(this)
}
