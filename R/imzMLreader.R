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


#' import_imzML.
#'
#' @param imzML_File full path to .imzML file.
#' @param ibd_File path to the binary file (default the same as imzML file but with .ibd extension)
#' @param ramdisk_path where the ramdisk will be created.
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param fun_text This is a callback function to update the label widget of loading data. See details for more information.
#' @param close_signal function to be called if loading process is abored.
#' @param verifyChecksum if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#'
#'  Imports an imzML image to an rMSI data object.
#'  It is recomanded to use rMSI::LoadMsiData directly instead of this function.
#'
#' @return an rMSI data object.
#' @export
#'
import_imzML <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ), ramdisk_path,  fun_progress = NULL, fun_text = NULL, close_signal = NULL, verifyChecksum = T)
{
  setPbarValue<-function(progress)
  {
    setTxtProgressBar(pb, progress)
    return(T)
  }

  if(is.null(fun_progress))
  {
    pb<-txtProgressBar(min = 0, max = 100, style = 3 )
    fun_progress <- setPbarValue
    cat("\n")
  }
  else
  {
    pb<-NULL
  }

  #1- Parse XML data
  if(is.null(fun_text))
  {
    cat("Parsing XML data in imzML file...\n")
  }
  else
  {
    fun_text("Parsing XML data in imzML file...")
  }

  xmlRes <- CimzMLParse(path.expand(imzML_File))
  if( !is.null(xmlRes$Error))
  {
    .controlled_loadAbort(paste0(xmlRes$Error, "\n"), close_signal)
  }
  
  #TODO test Checksum!
  
  #2- Create a connection to read binary file
  bincon <- file(description = ibd_File, open = "rb")

  #3- Test the UUID in binary file (the first 16 bytes are always UUID (in XML file are in hex codes))
  binUUID <- paste(sprintf("%.2X", readBin(bincon, integer(), 16, size = 1, signed = F)), collapse = "")
  if(binUUID != xmlRes$UUID)
  {
    close(bincon)
    .controlled_loadAbort("ERROR: UUID in imzML file does not match UUID in ibd file\n", close_signal)
  }

  sizeInBytesFromDataType <- function(str_dataType)
  {
    if(str_dataType == "int" || str_dataType == "float")
    {
      return(4)
    }
    if(str_dataType == "long" || str_dataType == "double")
    {
      return(8)
    }
  }

  #4- Obtain de m/z axis (this is only valid for imzML continuous mode)
  pt<-proc.time()
  if(xmlRes$mz_dataType == "int" || xmlRes$mz_dataType == "long")
  {
    readDataTypeMz <- integer()
  }
  if(xmlRes$mz_dataType == "float" || xmlRes$mz_dataType == "double")
  {
    readDataTypeMz <- numeric()
  }
  bytes2ReadMz <- sizeInBytesFromDataType(xmlRes$mz_dataType)
  
  if(xmlRes$continuous_mode)
  {
    mzAxis <- readBin(bincon, readDataTypeMz, xmlRes$run_data[1, "mzLength"], size = bytes2ReadMz, signed = T)
  }
  else
  {
    #Processed mode, so a common mass axis must be calculated and stored in mzAxis var
    cat("Calculating the new mass axis...\n")
    if(!is.null(fun_text))
    {
      fun_text("Process mode, re-calculating mass axis...")
    }
    ppStep<-100/nrow(xmlRes$run_data)
    pp<-0
    #Update progress bar
    if( !fun_progress(pp) )
    {
      return(NULL) #progress bar function must return true if the loading process is to be continued.
    }
    
    #Fill the initial spectrum as the larger one
    seli <- which.max(xmlRes$run_data$mzLength)
    seek(bincon, rw = "read", where = xmlRes$run_data[seli, "mzOffset"] )
    mzAxis <- readBin(bincon, readDataTypeMz, xmlRes$run_data[seli, "mzLength"], size = bytes2ReadMz, signed = T)
    mzAxis <- unique(mzAxis) #Avoid duplicates
    for( i in 1:nrow(xmlRes$run_data))
    {
      #Read mass axis for the current spectrum 
      seek(bincon, rw = "read", where = xmlRes$run_data[i, "mzOffset"] )
      mzdd <- readBin(bincon, readDataTypeMz, xmlRes$run_data[i, "mzLength"], size = bytes2ReadMz, signed = T)
      mzdd <- unique(mzdd) #Avoid duplicates
      
      #Combine the two mass axis using Cpp method
      mzAxis <- MergeMassAxis(mzAxis, mzdd)
      
      #Update progress bar
      pp_ant<-pp
      pp<-pp+ppStep
      if(!is.null(fun_progress) && (round(pp) > round(pp_ant)) )
      {
        #Update progress bar
        if( !fun_progress(pp) )
        {
          return(NULL) #progress bar function must return true if the loading process is to be continued.
        }
      }
    }
  }
  pt<-proc.time() - pt
  cat(paste("\nMass axis calculation time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  cat(paste("The re-sampled mass axis contains", length(mzAxis), "data points\n"))

  #5- Map imzML possible data types to ff packages available data types
  if(xmlRes$int_dataType == "int")
  {
    ffDataType <- "integer" #32 bit signed integer with NA.
    readDataTypeInt <- integer()
  }
  if(xmlRes$int_dataType == "long")
  {
    ffDataType <- "single" #32 bit signed integer is not available in ff so I map it to 32 bits float to allow enough range.
    readDataTypeInt <- integer()
  }
  if(xmlRes$int_dataType == "float")
  {
    ffDataType <- "single"
    readDataTypeInt <- numeric()
  }
  if(xmlRes$int_dataType == "double")
  {
    ffDataType <- "double"
    readDataTypeInt <- numeric()
  }
  bytes2ReadInt <- sizeInBytesFromDataType(xmlRes$int_dataType)

  #6- Create the ramdisk
  if(!is.null(fun_text))
  {
    fun_text("Loading data in the ramdisk...")
  }
  dir.create(ramdisk_path, recursive = T)
  datacube <- CreateEmptyImage(num_of_pixels = nrow(xmlRes$run_data), mass_axis = mzAxis, pixel_resolution = xmlRes$pixel_size_um,
                               img_name = basename(imzML_File),
                               ramdisk_folder = ramdisk_path,
                               data_type = ffDataType,
                               uuid = binUUID
                               )

  #7- Read all spectra
  pt <- proc.time()
  cat("\nReading spectra from binary file...\n")
  ppStep<-100/nrow(xmlRes$run_data)
  pp<-0
  #Update progress bar
  if( !fun_progress(pp) )
  {
    return(NULL) #progress bar function must return true if the loading process is to be continued.
  }

  #Ramdisk data accessors for full ff matrix acces (should be faster than accessing singles spectrum)
  currentDataCube <- 1
  currentDataIRow <- 1
  bufferMatrix <- matrix(nrow = nrow(datacube$data[[1]]), ncol = length(datacube$mass))
  for(i in 1:nrow(xmlRes$run_data))
  {
    #Read intensity of current spectrum
    seek(bincon, rw = "read", where = xmlRes$run_data[i, "intOffset"] )
    dd <- readBin(bincon, readDataTypeInt, xmlRes$run_data[i, "intLength"], size = bytes2ReadInt, signed = T)
    
    if(!xmlRes$continuous_mode)
    {
      #Read mass axis for the current spectrum 
      seek(bincon, rw = "read", where = xmlRes$run_data[i, "mzOffset"] )
      mzdd <- readBin(bincon, readDataTypeMz, xmlRes$run_data[i, "mzLength"], size = bytes2ReadMz, signed = T)
      
      #Delete duplicates and possible zero-drops errors (fixing Bruker's bugs in imzML)
      idup <- which(duplicated(mzdd))
      if(length(idup) > 0)
      {
        for( i in 1:length(idup))
        {
          dd[idup[i] - 1] <- max( dd[idup[i]], dd[idup[i] - 1])
        }
        mzdd <- mzdd[-idup]
        dd <- dd[-idup]
      }
      
      #Apply re-sampling (only if needed...)
      if( ! identical(mzdd, mzAxis))
      {
        dd <- (approx( x = mzdd, y = dd, xout = mzAxis, ties = "ordered", yleft = 0, yright = 0))$y
        dd[which(is.na(dd))] <- 0 #Remove any possible NA
      } 
    }
    
    #Store the intensities to the ramdisk by blocks using a buffer matrix for better performance
    bufferMatrix[currentDataIRow, ] <- dd
    currentDataIRow <- currentDataIRow  + 1
    if( currentDataIRow > nrow(datacube$data[[currentDataCube]]))
    {
      datacube$data[[currentDataCube]][,] <- bufferMatrix
      currentDataCube <-  currentDataCube + 1
      currentDataIRow <- 1
      if(currentDataCube <= length(datacube$data))
      {
        bufferMatrix <- matrix(nrow = nrow(datacube$data[[currentDataCube]]), ncol = length(datacube$mass))
      }
    }
 
    datacube$pos[i, "x"] <- xmlRes$run_data[i, "x"]
    datacube$pos[i, "y"] <- xmlRes$run_data[i, "y"]

    #Update progress bar
    pp_ant<-pp
    pp<-pp+ppStep
    if(!is.null(fun_progress) && (round(pp) > round(pp_ant)) )
    {
      #Update progress bar
      if( !fun_progress(pp) )
      {
        return(NULL) #progress bar function must return true if the loading process is to be continued.
      }
    }
  }
  pt<-proc.time() - pt
  cat(paste("\nBinary file reading time:",round(pt["elapsed"], digits = 1),"seconds\n"))

  #8- Close the bin file connection
  close(bincon)
  cat("\n")

  #9- Compute image size and arrange coords to avoid holes
  datacube$posMotors <- datacube$pos
  datacube$pos <- remap2ImageCoords(datacube$pos)
  datacube$size["x"] <- max(datacube$pos[,"x"])
  datacube$size["y"] <- max(datacube$pos[,"y"])

  #10- Compute average spectrum
  datacube$mean <- AverageSpectrum(datacube)

  #11- Save the data cube for further fast access
  save(datacube, file = file.path(ramdisk_path, "datacube.RImg"))

  #11- And it's done, just return de rMSI object
  if(!is.null(pb))
  {
    close(pb)
  }
  gc()
  if(!is.null(fun_text))
  {
    fun_text("Done")
  }
  
  class(datacube) <- "rMSIObj"
  return(datacube)
}
