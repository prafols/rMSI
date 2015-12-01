#Easy import Bruker image files from a XMASS data directory and a xml file specifying raster spots in Bruker format
#Just source this file in R and select the XML file, the data dir and where to store the resulting tar.gz file

importBrukerXMASSImg_Wizard <- function()
{
  startWD <- getwd()

  #Ask for the data path
  cat("Select RAW data directory...\n\n")
  path_data <- gfile(text = "Select RAW data directory", type = "selectdir", multi = F)
  if(length(path_data) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }
  setwd(path_data)

  #Where data will be stored, the tar.gz filename will be the same as xml file is
  cat("Select output directory...\n\n")
  path_output <- gfile(text = "Select output directory", type = "selectdir", multi = F)
  if(length(path_output) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }


  #Ask for the XML file
  cat("Select XML file...\n\n")
  path_xml <- gfile(text = "Select XML file", type = "open", filter = c("xml" = "xml"), multi = T)
  if(length(path_xml) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }
  path_output_file<- file.path(path_output, paste(tools::file_path_sans_ext(basename(path_xml)),".tar",sep = "" ))

  #Propmt user to proceed
  cat("Data will be imported with the following settings:\n")
  cat(paste("RAW Data directory:\n\t",path_data, "\n", sep =""))
  cat("Output files:\n")
  for(i in 1:length(path_xml))
  {
    cat(paste("\t[", i,"]\t",path_output_file[i], "\n", sep = ""))
  }

  cat("XML files:\n")
  for(i in 1:length(path_xml))
  {
    cat(paste("\t[", i,"]\t",path_xml[i], "\n", sep = ""))
  }
  cat("\n\n")
  resp <- ""
  while(resp != "y" && resp != "n" && resp != "Y" && resp != "N"){
    resp <- readline(prompt = "Is everything correct? [y, n]:")
    if(resp != "y" && resp != "n" && resp != "Y" && resp != "N")  {
      cat("Invalid response, valid responses are: y, n, Y and N. Try again.\n")
    }
  }

  if( resp == "y" || resp =="Y") {
    for(i in 1:length(path_xml))
    {
      cat(paste("\n\nStarting imporation of:", basename(path_xml[i]), "(file", i, "of", length(path_xml), ")\n"))
      veryStart<-proc.time()
      MSI2Rdata(raw_data_full_path = path_data, spot_selection_xml_file_full_path = path_xml[i], output_data_filename= path_output_file[i])
      cat(paste("Importation of", basename(path_xml[i]), "complete  with the following time statistics:\n"))
      print(proc.time() - veryStart)
    }
  } else {
    cat("Data importation aborted by user, please try again!\n")
  }

  setwd(startWD)
}
