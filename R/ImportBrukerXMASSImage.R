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


#Importation from Bruker XMASS data to custom data format based on ff pacakge

#This is the top level function for RAW datacuve import
#RAW Data strcutre matrix (stored in HDD), nrow = npixels, ncol = length(mz_axis)
# | mz1 | mz2 | ... | mzN |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
# |     |     | ... |     |
# +-----+-----+-----+-----+
#
#POS Data matrix (stored in RAM)
# | X  | Y  |
# +----+----+
# | x1 | y1 |
# +----+----+
# | x2 | y2 |
# +----+----+
# | .. | .. |
# +----+----+
# | xM | yM |
# +----+----+
#
#MZ_AXIS -> vector (stored in RAM)

#' Import a MSI dataset from Bruker XMASS folder and a XML file.
#'
#' This function converts MSI data contained in Bruker XMASS directory into  rMSI data structure.
#' Only spectra contained in the XML file are imported. The ramdisk will be created at the specified
#' directory as a collection of ff files with the maximum specified size (50 MB is the default setting).
#' The rMSI object pointing to ramdisk and containing RAM stored variables is returned.
#'
#' @param raw_data_full_path Where the Bruker XMASS data is located.
#' @param resolution_um The image pixel size un micrometers (sorry it can not be read directly from XMASS).
#' @param xml_file_full_path Full path to the XML file or alternatively an already parsed xml list resulting from ParseBrukerXML() function.
#' @param txt_spectrum_path Full path to a bruker spectrum txt file to extract mass axis (optional).
#' @param selected_img the image number to select in the XML file containing various classes.
#' @param ... Extra parameters to .readBrukerXmassImg ( for example max_ff_file_size_MB as max size in MB of each ff file of the ramdisk).
#'
#' @export
importBrukerXmassImg<-function(raw_data_full_path, resolution_um, xml_file_full_path, txt_spectrum_path = NULL, selected_img = 1, ...)
{
  cat("Importing data to R session...\n")
  pt<-proc.time()
  ff_folder<-file.path(raw_data_full_path, "ffdata")
  raw<-.readBrukerXmassImg(raw_data_folder = raw_data_full_path, xml_file = xml_file_full_path, ff_data_folder = ff_folder, sample_spectrum_path = txt_spectrum_path, selected_xml_class = selected_img, ...)

  pt<-proc.time() - pt
  cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))

  pt<-proc.time()
  raw$mean <-AverageSpectrum(raw)
  pt<-proc.time() - pt
  cat(paste("Average spectrun time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  gc()

  #Append resolution fields to img Object
  raw$pixel_size_um <- resolution_um

  class(raw) <- "rMSIObj"
  return(raw)
}

#Internal function too read bruker XMASS data and create the ramdisk
.readBrukerXmassImg<-function(raw_data_folder, xml_file, ff_data_folder, sample_spectrum_path = NULL, selected_xml_class = 1)
{
  #Internal methods
  ##Stores a directory structure of raw_data
  dataListing<-function(raw_data_folder, sample_spectra_name)
  {
    topDataDir<-dir(raw_data_folder, full.names = F)
    if( sample_spectra_name %in% topDataDir )
    {
      #No sub-data dirs, spectra is directly in raw_data_folder
      rawDirs<-dir(raw_data_folder, full.names = T)
      names(rawDirs)<-topDataDir
    }
    else
    {
      #There exists sub-data dirs
      rawDirs<-unlist(lapply(topDataDir ,function(x) dir(file.path(raw_data_folder, x), full.names = T)))
      names(rawDirs)<-unlist(lapply(topDataDir ,function(x) dir(file.path(raw_data_folder, x), full.names = F)))
    }
    return(rawDirs)
  }

  #Read Spectrum acqu file
  #Returns a vector of all possible m/z values
  generateMzAxisFromAcquFile<-function(acqu_file){
    acqu_con<-file(acqu_file, open="r")
    file_lines<-readLines(acqu_con)
    close(acqu_con)

    #Grab Mz axis generation data
    td<-as.numeric(strsplit(grep("##\\$TD=", file_lines, value=T),"=")[[1]][2])
    dw<-as.numeric(strsplit(grep("##\\$DW=", file_lines, value=T),"=")[[1]][2])
    delay<-as.numeric(strsplit(grep("##\\$DELAY=", file_lines, value=T),"=")[[1]][2])
    ml1<-as.numeric(strsplit(grep("##\\$ML1=", file_lines, value=T),"=")[[1]][2])
    ml2<-as.numeric(strsplit(grep("##\\$ML2=", file_lines, value=T),"=")[[1]][2])
    ml3<-as.numeric(strsplit(grep("##\\$ML3=", file_lines, value=T),"=")[[1]][2])

    tof<-delay + (0:(td-1))*dw
    #tof<-delay + (1:(td))*dw
    A<-ml3
    B<-sqrt(1e+12/ml1)
    C<-ml2 - tof
    mz<-((-B + sqrt((B * B) - (4 * A * C)))/(2 * A))^2

    return(mz)
  }

  #Read FID intensities (Bruker Formated)
  #Returns a vector of intensities
  readFidFile<-function(fid_file, points)
  {
    fid_con<-file(fid_file, open="rb")
    specVect<-as.integer(readBin(fid_con, integer(),n=points ,size = 4, signed=T, endian="little" ))
    close(fid_con)
    return(specVect)
  }


  #1- Read XML input data to get a vector of spectrums
  if( is.character(xml_file))
  {
    spectraListAllRois <- ParseBrukerXML(xml_file)  
    spectraList <- spectraListAllRois[[selected_xml_class]]
  }
  else
  {
    spectraList <- xml_file[[selected_xml_class]]
  }
  
  pb<-txtProgressBar(min = 0, max = length(spectraList$spots), style = 3 )
  #Store directory structure to avoid wasting time in dir() fuctions
  rawDirs<-dataListing(raw_data_folder, spectraList$spots[1])

  #2- Read Spectrum acqu file and generate m/z vector or use sampel txt spectrum if it is available
  mz_axis_acqu<-generateMzAxisFromAcquFile(file.path(rawDirs[spectraList$spots[1]], "1", "1SRef","acqu"))
  if(!is.null(sample_spectrum_path))
  {
    mz_axis_sample<-read.table(sample_spectrum_path)[,1]
    if(length(mz_axis_acqu) != length(mz_axis_sample))
    {
      stop("Error: The mass axis length of acqu file is different of the sample txt spectrum.\n")
    }
    mz_axis <- mz_axis_sample
  }
  else
  {
    #No sample spectrum provided, so use the acqu file
    mz_axis <- mz_axis_acqu
  }

  #3- Read Spectrum FID file and fill data structure if not zero intensity
  dataPos <- matrix(NA, ncol = 2, nrow = length(spectraList$spots))
  colnames(dataPos)<-c("x","y")
  
  UUID <- format(Sys.time(), "%Y%m%d%H%M%S")
  ff_data_folder <- file.path(ff_data_folder, paste0(UUID, "_", spectraList$name))
  dir.create(ff_data_folder, showWarnings = F, recursive = T)
  dataCube<-.CreateEmptyRamdisk(length(mz_axis), length(spectraList$spots), ff_data_folder )
  max_nrow <- nrow(dataCube[[1]])

  for(i in 1:length(spectraList$spots))
  {
    setTxtProgressBar(pb, i)
    i_cube<-(1+((i-1) %/% max_nrow))
    dataCube[[i_cube]][(i - (i_cube -1) * max_nrow), ]<-readFidFile(file.path(rawDirs[spectraList$spots[i]], "1", "1SRef","fid"), length(mz_axis))

    #Extract X Y Coords
    dataPos[i,"x"]<- spectraList$x[i]
    dataPos[i,"y"]<- spectraList$y[i]
  }
  close(pb)

  #4- Remap motor coords to image coords
  BrukerPos <- dataPos
  dataPos <- remap2ImageCoords(dataPos)
  x_size <- max(dataPos[,"x"])
  y_size <- max(dataPos[,"y"])

  #5- Return dataCube, mz_axis, xsize, ysize as a list of elements
  return(list(name = spectraList$name, mass = mz_axis, size = c(x = x_size, y = y_size), data = dataCube, pos=dataPos, posMotors = BrukerPos, uuid = UUID))
}

#' Import a MSI dataset from Bruker XMASS folder and a XML file using a Wizard.
#'
#' This function converts MSI data contained in Bruker XMASS directory into  rMSI data structure.
#' Only spectra contained in the XML file are imported. This method will prompt the user to select
#' the data direcotry, the XML files and output data path. A rMSI compressed images are stored as .tar
#' in specified output directory.
#' @export
importBrukerXMASSImg_Wizard <- function()
{
  startWD <- getwd()

  #Ask for the data path
  cat("Select RAW data directory...\n\n")
  path_data <- gWidgets2::gfile(text = "Select RAW data directory", type = "selectdir", multi = F)
  if(length(path_data) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }
  setwd(path_data)

  #Where data will be stored, the tar.gz filename will be the same as xml file is
  cat("Select output directory...\n\n")
  path_output <- gWidgets2::gfile(text = "Select output directory", type = "selectdir", multi = F)
  if(length(path_output) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }


  #Ask for the XML file
  cat("Select XML file...\n\n")
  path_xml <- gWidgets2::gfile(text = "Select XML file", type = "open", filter = c("xml" = "xml"), multi = T)
  if(length(path_xml) == 0){
    setwd(startWD)
    stop("Data importation aborted by user, please source this file again to restart it!\n", call. = F)
  }
  path_output_file<- file.path(path_output, paste(tools::file_path_sans_ext(basename(path_xml)),".tar",sep = "" ))

  #Ask for resolution
  resolutions <- rep(NA, length(path_xml))
  for( i in 1:length(path_xml))
  {
    while(T)
    {
      resp <- readline(prompt = paste("Pixel resolution in um of ", basename(path_xml[i]), ": ", sep = ""))
      resolutions[i] <- as.numeric(resp)
      if( !is.na(resolutions[i]) )
      {
        break
      }
      cat("Please provide a valid number for the resolution. Try again.")
    }

  }

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
    cat(paste("\t[", i,"]\t",path_xml[i], " resolution: ", resolutions[i], " um\n", sep = ""))
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
      rawImg <- importBrukerXmassImg(raw_data_full_path = path_data, resolution_um =  resolutions[i] , xml_file_full_path = path_xml[i])
      SaveMsiData(rawImg, path_output_file[i])
      DeleteRamdisk(rawImg)
      rm(rawImg)
      gc()
      cat(paste("Importation of", basename(path_xml[i]), "complete  with the following time statistics:\n"))
      print(proc.time() - veryStart)
    }
  } else {
    cat("Data importation aborted by user, please try again!\n")
  }

  gc()

  setwd(startWD)
}


