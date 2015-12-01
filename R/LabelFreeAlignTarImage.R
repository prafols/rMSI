#Easy Label-Free aligment tool
#Just source this file and select the *.tar image to align

LabelFreeAlign_Wizard<-function()
{
  startWD <- getwd()

  #Ask for the data path
  cat("Select a tar image to align...\n\n")
  path_in_img <- gfile(text = "Select a tar image to align", type = "open", filter = c("tar" = "tar"), multi = T)
  if(length(path_in_img) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }
  setwd(dirname(path_in_img))

  path_out_img <- paste(tools::file_path_sans_ext(path_in_img),"_Aligned",".tar",sep = "" )


  for( i in 1:length(path_in_img))
  {
    cat( paste("Starting aligning image:", path_in_img[i], "  ", i ,"of", length(path_in_img), "\n"))

    #Load the input image
    raw<-LoadImageWithProgressBar(path_in_img[i])

    #Algin the image
    FullDataSetAligment(raw, subDataMemPercent = 1)

    #Saving Data to compressed format
    SaveMsiData(data_file = path_out_img[i], imgData = raw, meanSpcData = raw$mean)

    #Remove the ramdisk
    RemoveImageRamdisk(raw)

    cat( paste("Image:", path_in_img[i], "Aligned to:", path_out_img[i], "\n\n\n"))
  }

  setwd(startWD)
}
