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
    cat("Parsing XML data in imzML file...")
  }
  else
  {
    fun_text("Parsing XML data in imzML file...")
  }
  
  xmlRes <- imzMLparse( imzML_File, ibd_File, fun_progress =  fun_progress, close_signal = close_signal, verifyChecksum = verifyChecksum)
  
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

#' imzMLparse.
#' Parses a imzML file and obtain all relevant data for rMSI data object.
#'
#' @param fpath_xml full path to imzML file (XML data).
#' @param fpath_bin full path to ibd file (binary data).
#' @param fun_progress callback to a pbar function.
#' @param close_signal function to call if error.
#' @param verifyChecksum if the binary file checksum must be verified, it can be disabled for convenice with really big files.
#'
#' @return a named list containing all relevant information for rMSI data object.
#' UUID: A string containing the UUID to be verified also in binary file.
#' continuous_mode: boolean indicating if binary data is in continuous mode (TRUE) or processed mode (FALSE).
#' compression_mz: boolean indicating if binary mz data is compressed (TRUE) or not (FALSE).
#' compression_int: boolean indicating if binary intensity data is compressed (TRUE) or not (FALSE).
#' mz_dataType: a string respresenting the data type used to encode mz axis (int, long, float or double).
#' int_dataType: a string respresenting the data type used to encode intensity data (int, long, float or double).
#' pixel_size_um: per pixel resolution in micrometers.
#' run_data: a data frame containing information of pixel location in image ans offsets in binary data file.
#'
imzMLparse <- function( fpath_xml, fpath_bin, fun_progress = NULL, close_signal = NULL, verifyChecksum = T)
{
  xmld <- XML::xmlTreeParse(fpath_xml)
  xmltop <- XML::xmlRoot(xmld)
  xmlch <- XML::xmlChildren(xmltop)

  #Obtain UUID, SHA, and mode
  xmlFDesc <- XML::xmlChildren(xmlch$fileDescription)
  UUID <- NULL
  SHA <- NULL
  MD5 <- NULL
  continuous_mode <- NULL
  for( i in 1:(XML::xmlSize(xmlFDesc$fileContent)))
  {
    id <- XML::xmlGetAttr(xmlFDesc$fileContent[[i]], "accession")
    value <- XML::xmlGetAttr(xmlFDesc$fileContent[[i]], "value")

    if(id == "IMS:1000080")
    {
      UUID <- toupper(paste(unlist(strsplit(unlist(strsplit(unlist(strsplit(value, "{", fixed = T)), "}", fixed = T)), "-", fixed = T)), collapse = ""))
    }
    if( id == "IMS:1000091")
    {
      SHA <- toupper(value)
    }
    if( id == "IMS:1000090")
    {
      MD5 <- toupper(value)
    }
    if( id == "IMS:1000030")
    {
      continuous_mode <- T
    }
    if( id == "IMS:1000031")
    {
      continuous_mode <- F
    }
  }

  #Check for NULL data
  err_string <- character(0)
  if(is.null(UUID))
  {
    err_string <- paste(err_string, "No UUID field found\n", sep = "")
  }
  if(is.null(continuous_mode))
  {
    err_string <- paste(err_string, "No continuos/processed field found\n", sep = "")
  }
  if(length( err_string) > 0)
  {
    cat(paste("imzML reading has been aborted due the following errors:\n", err_string, sep = ""))
    .controlled_loadAbort("imzML XML parese ERROR\n", close_signal)
  }

  #Check the checksum
  if(verifyChecksum)
  {
    if( !is.null(SHA))
    {
      cat("\nChecking binary data checksum using SHA-1 key... ")
      res <- toupper(digest::digest( fpath_bin, algo = "sha1", file = T))
      if( res == SHA )
      {
        cat("OK\n")
      }
      else
      {
        cat(paste("NOK\nChecksums don't match\nXML key:", SHA, "\nBinary file key:", res,"\n"))
        #.controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
        #Disableing the abort, just showing a warning here as it seams that there is bug in brukers checksum imzml file...
        cat("WARNING: MS data my be corrupt!\n")
      }
    }
    if( !is.null(MD5))
    {
      cat("Checking binary data checksum using MD5 key... ")
      res <- toupper(digest::digest( fpath_bin, algo = "md5", file = T))
      if( res == MD5 )
      {
        cat("OK\n")
      }
      else
      {
        cat(paste("NOK\nChecksums don't match\nXML key:", MD5, "\nBinary file key:", res,"\n"))
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
  
  #Obtain data size per bin
  parseXmlArrayDataType <- function(xmlArrNode)
  {
    no_compression <- F
    units <- NULL
    dataType <- NULL
    isMassArray <- F
    isIntensityArray <- F
    xmlSubDataType <- XML::xmlChildren(xmlArrNode)
    for( j in 1:(XML::xmlSize(xmlSubDataType)))
    {
      id <- XML::xmlGetAttr(xmlSubDataType[[j]], "accession")

      if( id == "MS:1000576") #test compression
      {
        no_compression <- T
      }

      if( id == "MS:1000514") #m/z array units
      {
        isMassArray <- T
        units <-  XML::xmlGetAttr(xmlSubDataType[[j]], "unitAccession")
      }

      if( id == "MS:1000515") #intensity array units
      {
        isIntensityArray <- T
        units <-  XML::xmlGetAttr(xmlSubDataType[[j]], "unitAccession")
      }

      #8 bits integer data array
      #TODO This is not encoded in mzML CV obo file but is specified as valid in imzML spec!

      #16 bits integer data array
      #TODO This is not encoded in mzML CV obo file but is specified as valid in imzML spec!

      #32 bits integer data array
      if(id == "IMS:1000141")
      {
        dataType <- "int"
      }

      #64 bits integer data array
      if(id == "IMS:1000142")
      {
        dataType <- "long"
      }

      #32 bits float data array
      if(id == "MS:1000521")
      {
        dataType <- "float"
      }

      #64 bits double data array
      if(id == "MS:1000523")
      {
        dataType <- "double"
      }

    }
    return(list( no_compression = no_compression, units = units, dataType = dataType, 
                 massArray = isMassArray, intensityArray = isIntensityArray ))
  }

  mzArrayDesc <- NULL
  intArrayDesc <- NULL
  xmlDataType <- XML::xmlChildren( xmlch$referenceableParamGroupList)
  for ( i in 1:(XML::xmlSize(xmlDataType)))
  {
    auxArrayDesc <- parseXmlArrayDataType(xmlDataType[[i]])
    auxArrayDesc$ref <- XML::xmlGetAttr(xmlDataType[[i]], "id")
    
    if(auxArrayDesc$massArray)
    {
      mzArrayDesc <- auxArrayDesc
    }
    if(auxArrayDesc$intensityArray)
    {
      intArrayDesc <- auxArrayDesc
    }
    rm(auxArrayDesc)
  }

  #Check for emty data
  err_string <- character(0)
  if( is.null(mzArrayDesc$units ))
  {
    err_string <- paste(err_string, "No m/z array units declared\n", sep = "")
  }
  if( is.null(intArrayDesc$units ))
  {
    cat("\nWarning: No intensity array units declared. Assuming number of counts as default unit.\n") #Just display a warining because some implementation do not provide data units for intensity array
    intArrayDesc$units <- "MS:1000131"
  }
  if( is.null(mzArrayDesc$dataType))
  {
    err_string <- paste(err_string, "No m/z array data type declared\n", sep = "")
  }
  if( is.null(intArrayDesc$dataType))
  {
    err_string <- paste(err_string, "No intensity array data type declared\n", sep = "")
  }
  if(length( err_string) > 0)
  {
    cat(paste("imzML reading has been aborted due the following errors:\n", err_string, sep = ""))
    .controlled_loadAbort("imzML XML parese ERROR\n", close_signal)
  }

  #Check data definitions fits the rMSI requirments
  if( !mzArrayDesc$no_compression || !intArrayDesc$no_compression)
  {
    .controlled_loadAbort("Data is compressed, compressed imzML images are still not supported for rMSI\n", close_signal)
  }

  if (mzArrayDesc$units != "MS:1000040")
  {
    .controlled_loadAbort("m/z Array is not in m/z format, rMSI only supports m/z units\n", close_signal)
  }

  if( intArrayDesc$units != "MS:1000131")
  {
    .controlled_loadAbort("intensity Arrays are not in number of detector counts, rMSI does not support it\n", close_signal)
  }

  #Obtain pixel resolution (I'm reading directly the accession becasue I've found some mismatch beetween obo and imzML example)
  pixel_size_um <- NULL
  xmlScanSettings <- XML::xmlChildren(xmlch$scanSettingsList)
  for( i in 1:(XML::xmlSize(xmlScanSettings$scanSettings)))
  {
    accession <- XML::xmlGetAttr(xmlScanSettings$scanSettings[[i]], "accession")
    if( accession == "IMS:1000046")
    {
      pixel_size_um <- sqrt(as.numeric(XML::xmlGetAttr(xmlScanSettings$scanSettings[[i]], "value"))) #Getting the sqrt value since the newest imzML standard povides the area of the pixel in microns
    }
  }

  #Obtain the RUN data
  imgData <- xmlParseSpectra(XML::xmlChildren(xmlch$run), mzArrayDesc$ref, intArrayDesc$ref,  fun_progress = fun_progress)

  return( list( UUID = UUID, continuous_mode = continuous_mode,
                compression_mz = !mzArrayDesc$no_compression, compression_int = !intArrayDesc$no_compression,
                mz_dataType = mzArrayDesc$dataType, int_dataType = intArrayDesc$dataType,
                pixel_size_um = pixel_size_um, run_data = imgData))
}

#' xmlParseSpectra.
#'
#' Extracts all important information to read the imzML binary file from the XML run node.
#' Internal use only.
#'
#' @param xmlRunNode an XML node containing the imzML run tree.
#' @param fun_progress a callback to use another progress bar instead of the internaly used text pbar.
#' @param mzArrayRef a string with the value declared to search for mz arrays.
#' @param inArrayRef a string with the value declared to search for intensity arrays.
#'
#' @return a named data frame containing pixels positions in image (x, y) and spectra poisitions in binary file (offsets).
#'
xmlParseSpectra <- function (xmlRunNode, mzArrayRef, inArrayRef, fun_progress = NULL)
{
  updatePbar <- function( value )
  {
    if(!is.null(fun_progress))
    {
      if(!fun_progress(value))
      {
        return(NULL) #progress bar function must return true if the loading process is to be continued.
      }
    }
    else
    {
      setTxtProgressBar(pb, i)
    }
  }

  if(is.null(fun_progress))
  {
    pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  }

  pt<-proc.time()
  imgData <- data.frame(x = numeric(0), y = numeric(0), mzLength = integer(0), mzOffset = integer(0), intLength = integer(0), intOffset = integer(0))
  cat("Reading imzML spectra data...\n")
  ppStep<-100/(XML::xmlSize(xmlRunNode$spectrumList))
  pp<-0
  for( i in 1:(XML::xmlSize(xmlRunNode$spectrumList)))
  {
    imgData[i, ] <- unlist(xmlParseSpectrum(xmlRunNode$spectrumList[[i]], mzArrayRef, inArrayRef))
    pp_ant<-pp
    pp<-pp+ppStep
    if(round(pp) > round(pp_ant))
    {
      updatePbar( pp )
    }
  }

  if(is.null(fun_progress))
  {
    close(pb)
  }
  pt<-proc.time() - pt
  cat(paste("\nimzML XML parse time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  return( imgData )
}

#' xmlParseSpectrum.
#'
#' Parses an XML node in imzML format containing a single spectrum. All meaning data is returned as named list.
#' Internal use only.
#'
#' @param xmlSingleSpectrum a xml node containing data for a single spectrum.
#' @param mzSRef a string with the naming reference used to search mz array.
#' @param inSRef a string with the naming reference used to search intensity array.
#'
#' @return a named list containing all parameters of the spectrum.
#'
xmlParseSpectrum <- function( xmlSingleSpectrum, mzSRef, inSRef)
{
  #Read child data in spectrum node
  ttscan <- XML::xmlChildren(xmlSingleSpectrum)

  #Read pixel position
  ttscanPos <- XML::xmlChildren(ttscan$scanList)

  if( is.null((XML::xmlChildren(ttscanPos$scan))$referenceableParamGroupRef ))
  {
    ttscanPosXY <- ttscanPos$scan
  }
  else
  {
    ttscanPosXY <- XML::removeChildren(ttscanPos$scan, "referenceableParamGroupRef") #Remove unused node to facilitat reading
  }

  spectrumInfo <- list (x = NULL, y = NULL, mzLength = NULL, mzOffset = NULL, intLength = NULL, intOffset = NULL)
  for( i in 1:(XML::xmlSize(ttscanPosXY) ))
  {
    if(XML::xmlGetAttr(ttscanPosXY[[i]], "accession") == "IMS:1000050") #position x
    {
      spectrumInfo$x <- as.numeric(XML::xmlGetAttr(ttscanPosXY[[i]], "value"))
    }
    if(XML::xmlGetAttr(ttscanPosXY[[i]], "accession") ==  "IMS:1000051") #position y
    {
      spectrumInfo$y <- as.numeric(XML::xmlGetAttr(ttscanPosXY[[i]], "value"))
    }
    if( !is.null(spectrumInfo$x) && !is.null(spectrumInfo$y))
    {
      break
    }
  }

  #Read binary data info.
  ttscanBin <- XML::xmlChildren(ttscan$binaryDataArrayList)
  for( i in 1:(XML::xmlSize(ttscanBin)))
  {
    subNode <- XML::xmlChildren(ttscanBin[[i]])
    subNodeRef <- XML::xmlGetAttr(subNode$referenceableParamGroupRef, "ref")
    subNode <- XML::xmlChildren(XML::removeChildren(ttscanBin[[i]], "referenceableParamGroupRef", "binary")) #Remove unused node to facilitat reading

    for( j in 1:(XML::xmlSize(subNode)) )
    {
      access <- XML::xmlGetAttr(subNode[[j]], "accession")
      val <- as.numeric(XML::xmlGetAttr(subNode[[j]], "value"))

      if (subNodeRef == mzSRef)
      {
        if( access == "IMS:1000103") #external array length
        {
          spectrumInfo$mzLength <- val
        }
        if( access == "IMS:1000102") #external offset
        {
          spectrumInfo$mzOffset <- val
        }
      }

      if (subNodeRef == inSRef)
      {
        if( access == "IMS:1000103")  #external array length
        {
          spectrumInfo$intLength <- val
        }
        if( access == "IMS:1000102") #external offset
        {
          spectrumInfo$intOffset <- val
        }
      }

      if(  !is.null(spectrumInfo$mzLength) && !is.null(spectrumInfo$mzOffset) && !is.null(spectrumInfo$intLength) && !is.null(spectrumInfo$intOffset))
      {
        #I have all need data so stop now.
        break
      }
    }
    if(  !is.null(spectrumInfo$mzLength) && !is.null(spectrumInfo$mzOffset) && !is.null(spectrumInfo$intLength) && !is.null(spectrumInfo$intOffset))
    {
      #I have all need data so stop now, break the top loop.
      break
    }
  }
  return(spectrumInfo)
}

