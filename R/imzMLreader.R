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
#' @param verifyChecksum if the binary file checksum must be verified, it is disabled by default for convenice with really big files.
#' @param subImg_rename alternative image name, a new ramdisk will be created with the given name.
#' @param subImg_Coords a Complex vector with the motors coordinates to be included in the ramdisk.
#' @param convertProcessed2Continuous if true (the default) an imzML file in processed mode will be converted to a continuous mode.
#'
#'  Imports an imzML image to an rMSI data object.
#'  It is recomanded to use rMSI::LoadMsiData directly instead of this function.
#'
#' @return an rMSI data object.
#' @export
#'
import_imzML <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ), ramdisk_path,  fun_progress = NULL, fun_text = NULL, close_signal = NULL, verifyChecksum = F, subImg_rename = NULL, subImg_Coords = NULL, convertProcessed2Continuous = T)
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
  if(verifyChecksum)
  {
    if( xmlRes$SHA != "" )
    {
      cat("\nChecking binary data checksum using SHA-1 key... ")
      res <- toupper(digest::digest( ibd_File, algo = "sha1", file = T))
      if( res == xmlRes$SHA )
      {
        cat("OK\n")
      }
      else
      {
        cat(paste("NOK\nChecksums don't match\nXML key:", xmlRes$SHA, "\nBinary file key:", res,"\n"))
        #.controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
        #Disableing the abort, just showing a warning here as it seams that there is bug in brukers checksum imzml file...
        cat("WARNING: MS data my be corrupt!\n")
      }
    }
    if( xmlRes$MD5 != "")
    {
      cat("Checking binary data checksum using MD5 key... ")
      res <- toupper(digest::digest( ibd_File, algo = "md5", file = T))
      if( res == xmlRes$MD5 )
      {
        cat("OK\n")
      }
      else
      {
        cat(paste("NOK\nChecksums don't match\nXML key:", xmlRes$MD5, "\nBinary file key:", res,"\n"))
        #.controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
        #Disableing the abort, just showing a warning here as it seams that there is bug in brukers checksum imzml file...
        cat("WARNING: MS data my be corrupt!\n")
      }
    }  
  }
  else
  {
    cat("WARNING: Checksum validation is disabled, data may be corrupt!\n")
  }
  
  #Select specific pixels in the dataset
  if( !is.null(subImg_Coords))
  {
    keepIds <-  which( complex( real = xmlRes$run_data$x, imaginary = xmlRes$run_data$y) %in% subImg_Coords)

    if(length(keepIds) == 0)
    {
      .controlled_loadAbort("ERROR: no subImg_Coords found in current imzML data.\n", close_signal)
    }
    
    xmlRes$run_data <- xmlRes$run_data[keepIds,]
  }
  
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
  
  bCreateRamdisk <- F #Start assuming data in processed and peak list
  if(xmlRes$continuous_mode)
  {
    mzAxis <- readBin(bincon, readDataTypeMz, xmlRes$run_data[1, "mzLength"], size = bytes2ReadMz, signed = T)
    bCreateRamdisk <- T
  }
  else
  {
    if(convertProcessed2Continuous)
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
      mzMergeErrorCount <- 0 #Count merge mass axis errors
      for( i in 1:nrow(xmlRes$run_data))
      {
        #Read mass axis for the current spectrum 
        seek(bincon, rw = "read", where = xmlRes$run_data[i, "mzOffset"] )
        mzdd <- readBin(bincon, readDataTypeMz, xmlRes$run_data[i, "mzLength"], size = bytes2ReadMz, signed = T)
        mzdd <- unique(mzdd) #Avoid duplicates
        
        #Combine the two mass axis using Cpp method
        resMZMerge <- MergeMassAxis(mzAxis, mzdd)
        if(resMZMerge$error)
        {
          mzMergeErrorCount <- mzMergeErrorCount + 1
        }
        else
        {
          mzAxis <- resMZMerge$mass
        }
        
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
      
      #Check the resulting mass axis looking at errors:
      if ( mzMergeErrorCount >= nrow(xmlRes$run_data))
      {
        stop("Error: The mass axis of two vectors to merge is not compatible because they do not share a common range.");    
      }
        
      bCreateRamdisk <- T
    }
  }
  pt<-proc.time() - pt
  if(bCreateRamdisk)
  {
    cat(paste("\nMass axis calculation time:",round(pt["elapsed"], digits = 1),"seconds\n"))
    cat(paste("The re-sampled mass axis contains", length(mzAxis), "data points\n"))
  }
  
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

  #6- Create the ramdisk only if data is going to be coerced to continuous mode
  if( is.null(subImg_rename))
  {
    subImg_rename <-  basename(imzML_File)
  }
  else
  {
    ramdisk_path <- paste0(ramdisk_path, "_",subImg_rename )
  }
  if(bCreateRamdisk)
  {
    if(!is.null(fun_text))
    {
      fun_text("Loading data in the ramdisk...")
    }
    dir.create(ramdisk_path, recursive = T, showWarnings = F)
    datacube <- CreateEmptyImage(num_of_pixels = nrow(xmlRes$run_data), mass_axis = mzAxis, pixel_resolution = xmlRes$pixel_size_um,
                                 img_name = subImg_rename,
                                 ramdisk_folder = ramdisk_path,
                                 data_type = ffDataType,
                                 uuid = binUUID
                                 )
    #Ramdisk data accessors for full ff matrix acces (should be faster than accessing singles spectrum)
    currentDataCube <- 1
    currentDataIRow <- 1
    bufferMatrix <- matrix(nrow = nrow(datacube$data[[1]]), ncol = length(datacube$mass))
  }
  else
  {
    #No ramdisk is created then datacube will be a dummy list to store positions and peak lists
    datacube <- list(name = subImg_rename, 
                     size = c(NA,NA),
                     uuid = binUUID,
                     pos = matrix(NA, nrow = nrow(xmlRes$run_data), ncol = 2),
                     pixel_size_um = xmlRes$pixel_size_um,
                     peaks = list())
    colnames(datacube$pos) <- c("x", "y")
    names(datacube$size) <- c("x", "y")
  }
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
      
      if(bCreateRamdisk)
      {
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
          if(length(mzdd) == length(dd))
          {
            dd <- (approx( x = mzdd, y = dd, xout = mzAxis, ties = "ordered", yleft = 0, yright = 0))$y
            dd[which(is.na(dd))] <- 0 #Remove any possible NA
          }
          else
          {
            cat(paste0("WARNING: spectra at X = ", xmlRes$run_data$x[i], "  Y = ", xmlRes$run_data$y[i], " corrupt!\n" ))
            dd <- rep(0, length(mzAxis))  
          }
        } 
      }
    }
    
    if(bCreateRamdisk)
    {
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
    }
    else
    {
      datacube$peaks[[i]]<-list(mass = mzdd, intensity = dd)
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

  if(bCreateRamdisk)
  {
    #10- Compute average spectrum
    datacube$mean <- AverageSpectrum(datacube)
  
    #11- Save the data cube for further fast access
    save(datacube, file = file.path(ramdisk_path, "datacube.RImg"))
    class(datacube) <- "rMSIObj"
  }
  else
  {
    class(datacube) <- "peakList"
  }
  
  #12- And it's done, just return de rMSI object
  if(!is.null(pb))
  {
    close(pb)
  }
  gc()
  if(!is.null(fun_text))
  {
    fun_text("Done")
  }
  
  return(datacube)
}
