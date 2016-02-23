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
#' @param xml_file_full_path Full path to the XML file.
#' @param output_data_filename Where the compressed image .tar will be stored
#' @param ... Extra parameters to .readBrukerXmassImg ( for example max_ff_file_size_MB as max size in MB of each ff file of the ramdisk).
#'
#' @export
importBrukerXmassImg<-function(raw_data_full_path, resolution_um, xml_file_full_path, output_data_filename, ...)
{
  cat("Importing data to R session...\n")
  pt<-proc.time()
  ff_folder<-file.path(raw_data_full_path, "ffdata")
  dir.create(ff_folder)
  raw<-.readBrukerXmassImg(raw_data_folder = raw_data_full_path, xml_file = xml_file_full_path, ff_data_folder = ff_folder, ...)

  pt<-proc.time() - pt
  cat(paste("Importing time:",round(pt["elapsed"], digits = 1),"seconds\n"))

  cat("Calculating average spectrum...\n")
  pt<-proc.time()
  raw$mean <-AverageSpectrum(raw)
  pt<-proc.time() - pt
  cat(paste("Average spectrun time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  gc()

  cat("Packaging R data objects...\n")
  pt<-proc.time()

  #Append esolution fields to img Object
  raw$pixel_size_um <- resolution_um

  #Store the img to hdd
  SaveMsiData(raw, output_data_filename)
  pt<-proc.time() - pt
  cat(paste("Data saving time:", round(pt["elapsed"], digits = 1),"seconds\n"))

  cat("Done. Data is saved in:\n",  output_data_filename, "\n")
  gc()

  #Remove ff file created on HDD
  lapply(raw$data, ff::delete)
  rm(raw)
  unlink(ff_folder, recursive = T) #Remove folder containing ffdata
}


#Internal function too read bruker XMASS data and create the ramdisk
.readBrukerXmassImg<-function(raw_data_folder, xml_file, ff_data_folder, max_ff_file_size_MB = 50)
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

  #Read XML input data
  #Returns a vector of selected spectrums in a given ROI
  readXmlSpectraList<-function(xml_file){
    xml_data<-XML::xmlToList(xml_file)
    #spots<-xml_data$Class[[1]]["Spot"]
    spots<-c()
    for (i in 1:length(xml_data$Class))
    {
      if(!is.na(xml_data$Class[[i]]["Spot"]))
      {
        spots<-c(spots,xml_data$Class[[i]]["Spot"])
      }
    }
    names(spots)<-NULL
    return(spots)
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
  spectraList<-readXmlSpectraList(xml_file)
  pb<-txtProgressBar(min = 0, max = length(spectraList), style = 3 )
  #Store directory structure to avoid wasting time in dir() fuctions
  rawDirs<-dataListing(raw_data_folder, spectraList[1])

  #2- Read Spectrum acqu file and generate m/z vector
  mz_axis<-generateMzAxisFromAcquFile(file.path(rawDirs[spectraList[1]], "1", "1SRef","acqu"))

  #3- Read Spectrum FID file and fill data structure if not zero intensity
  dataPos <- matrix(NA, ncol = 2, nrow = length(spectraList))
  colnames(dataPos)<-c("x","y")
  dataCube<-.CreateEmptyRamdisk(length(mz_axis), length(spectraList), ff_data_folder, max_ff_file_size_MB )
  max_nrow <- nrow(dataCube[[1]])

  for(i in 1:length(spectraList))
  {
    setTxtProgressBar(pb, i)
    i_cube<-(1+((i-1) %/% max_nrow))
    dataCube[[i_cube]][(i - (i_cube -1) * max_nrow), ]<-readFidFile(file.path(rawDirs[spectraList[i]], "1", "1SRef","fid"), length(mz_axis))

    #Extract X Y Coords
    dataPos[i,"x"]<-as.integer(strsplit(strsplit(spectraList[i],"X")[[1]][2], "Y")[[1]][1])
    dataPos[i,"y"]<-as.integer(strsplit(strsplit(spectraList[i],"Y")[[1]][2], "X")[[1]][1])
  }
  close(pb)

  #4- Calc offsets and subtract it
  x_offset<-min(dataPos[,"x"])
  y_offset<-min(dataPos[,"y"])
  for(i in 1:length(spectraList))
  {
    dataPos[i, "x"] <- dataPos[i, "x"] - x_offset + 1
    dataPos[i, "y"] <- dataPos[i, "y"] - y_offset + 1
  }

  #5- Compute Motor coords range
  x_size<-max(dataPos[,"x"])
  y_size<-max(dataPos[,"y"])

  #6- Map MALDI motor coords to image cords (1-pixels steps)
  #It is important to map MALDI motors coords to image coords.
  #Otherwise, null extra pixels may be added leading to bad reconstruction
  px_map <- matrix( 0, nrow = x_size, ncol = y_size)
  for(i in 1:nrow(dataPos))
  {
    xi <- dataPos[i, "x"]
    yi <- dataPos[i, "y"]
    px_map[xi, yi]<- i
  }

  colNull <- which( base::colSums(px_map) == 0)
  rowNull <- which( base::rowSums(px_map) == 0)
  remap<-FALSE
  if( length(colNull) > 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) > 0 && length(rowNull) == 0 )
  {
    px_map_ <- px_map[ , -colNull ]
    remap<-TRUE
  }
  if( length(colNull) == 0 && length(rowNull) > 0 )
  {
    px_map_ <- px_map[ -rowNull , ]
    remap<-TRUE
  }

  if(remap)
  {
    for(ix in 1:nrow(px_map_))
    {
      for(iy in 1:ncol(px_map_))
      {
        if(px_map_[ix, iy] > 0)
        {
          dataPos[px_map_[ix, iy], "x"] <- ix
          dataPos[px_map_[ix, iy], "y"] <- iy
        }
      }
    }

    x_size <- max(dataPos[,"x"])
    y_size <- max(dataPos[,"y"])
  }

  #7- Return dataCube, mz_axis, xsize, ysize as a list of elements
  return(list(mass = mz_axis, size = c(x = x_size, y = y_size), data = dataCube, pos=dataPos))
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
      importBrukerXmassImg(raw_data_full_path = path_data, resolution_um =  resolutions[i] , xml_file_full_path = path_xml[i], output_data_filename= path_output_file[i])
      cat(paste("Importation of", basename(path_xml[i]), "complete  with the following time statistics:\n"))
      print(proc.time() - veryStart)
    }
  } else {
    cat("Data importation aborted by user, please try again!\n")
  }

  setwd(startWD)
}


