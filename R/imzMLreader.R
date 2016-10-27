#' import_imzML.
#'
#' @param imzML_File full path to .imzML file.
#' @param ibd_File path to the binary file (default the same as imzML file but with .ibd extension)
#' @param ramdisk_path where the ramdisk will be created.
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param fun_text This is a callback function to update the label widget of loading data. See details for more information.
#' @param close_signal function to be called if loading process is abored.
#'
#'  Imports an imzML image to an rMSI data object.
#'  It is recomanded to use rMSI::LoadMsiData directly instead of this function.
#'
#' @return an rMSI data object.
#' @export
#'
import_imzML <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ), ramdisk_path,  fun_progress = NULL, fun_text = NULL, close_signal = NULL)
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
  }
  else
  {
    pb<-NULL
  }

  #1- Parse XML data
  if(!is.null(fun_text))
  {
    fun_text("Parsing XML data in imzML file...")
  }
  xmlRes <- imzMLparse( imzML_File, ibd_File, fun_progress =  fun_progress, close_signal = close_signal)

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
  if(xmlRes$mz_dataType == "int" || xmlRes$mz_dataType == "long")
  {
    readDataType <- integer()
  }
  if(xmlRes$mz_dataType == "float" || xmlRes$mz_dataType == "double")
  {
    readDataType <- numeric()
  }
  mzAxis <- readBin(bincon, readDataType, xmlRes$run_data[1, "mzLength"], size = sizeInBytesFromDataType(xmlRes$mz_dataType), signed = T)

  #5- Map imzML possible data types to ff packages available data types
  if(xmlRes$int_dataType == "int")
  {
    ffDataType <- "integer" #32 bit signed integer with NA.
    readDataType <- integer()
  }
  if(xmlRes$int_dataType == "long")
  {
    ffDataType <- "single" #32 bit signed integer is not available in ff so I map it to 32 bits float to allow enough range.
    readDataType <- integer()
  }
  if(xmlRes$int_dataType == "float")
  {
    ffDataType <- "single"
    readDataType <- numeric()
  }
  if(xmlRes$int_dataType == "double")
  {
    ffDataType <- "double"
    readDataType <- numeric()
  }
  bytes2Read <- sizeInBytesFromDataType(xmlRes$int_dataType)

  #6- Create the ramdisk
  if(!is.null(fun_text))
  {
    fun_text("Loading data in the ramdisk...")
  }
  dir.create(ramdisk_path, recursive = T)
  datacube <- CreateEmptyImage(num_of_pixels = nrow(xmlRes$run_data), mass_axis = mzAxis, pixel_resolution = xmlRes$pixel_size_um,
                               img_name = basename(imzML_File),
                               ramdisk_folder = ramdisk_path,
                               data_type = ffDataType)

  #7- Read all spectra
  binReadOffset <- 16 + ((xmlRes$run_data[1, "mzLength"]) * sizeInBytesFromDataType(xmlRes$mz_dataType))
  cat("\nReading spectra from binary file...\n")
  ppStep<-100/nrow(xmlRes$run_data)
  pp<-0
  for( i in 1:nrow(xmlRes$run_data))
  {
    iCol <- which(xmlRes$run_data[, "intOffset"] == binReadOffset)
    dd <- readBin(bincon, readDataType, xmlRes$run_data[iCol, "intLength"], size = bytes2Read, signed = T)
    binReadOffset <- binReadOffset + (bytes2Read * xmlRes$run_data[iCol, "intLength"])
    saveImgCunckAtIds(datacube, Ids =  i, dm = matrix(dd, nrow = 1, ncol = length(dd)))
    datacube$pos[i, "x"] <- xmlRes$run_data[iCol, "x"]
    datacube$pos[i, "y"] <- xmlRes$run_data[iCol, "y"]

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

  #8- Close the bin file connection
  close(bincon)

  #9- Compute image size and arrange coords to avoid holes
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
  return(datacube)
}

#' imzMLparse.
#' Parses a imzML file and obtain all relevant data for rMSI data object.
#'
#' @param fpath_xml full path to imzML file (XML data).
#' @param fpath_bin full path to ibd file (binary data).
#' @param fun_progress callback to a pbar function.
#' @param close_signal function to call if error.
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
imzMLparse <- function( fpath_xml, fpath_bin, fun_progress = NULL, close_signal = NULL)
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
    name <- XML::xmlGetAttr(xmlFDesc$fileContent[[i]], "name")
    id <- XML::xmlGetAttr(xmlFDesc$fileContent[[i]], "accession")
    value <- XML::xmlGetAttr(xmlFDesc$fileContent[[i]], "value")

    if(name == "universally unique identifier" && id == "IMS:1000080")
    {
      UUID <- toupper(paste(unlist(strsplit(unlist(strsplit(unlist(strsplit(value, "{", fixed = T)), "}", fixed = T)), "-", fixed = T)), collapse = ""))
    }
    if( name == "ibd SHA-1" && id == "IMS:1000091")
    {
      SHA <- toupper(value)
    }
    if( name == "ibd MD5" && id == "IMS:1000090")
    {
      MD5 <- toupper(value)
    }
    if( name == "continuous" && id == "IMS:1000030")
    {
      continuous_mode <- T
    }
    if( name == "processed" && id == "IMS:1000031")
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

  #Check data in continuos mode, processed is currently not supported for rMSI
  if( !continuous_mode )
  {
    .controlled_loadAbort("Binary data is in imzML processed mode. Sorry this mode is not supported yet\n", close_signal)
  }

  #Check the checksum
  bChecked <- F
  if( !is.null(SHA))
  {
    cat("Checking binary data checksum using SHA-1 key... ")
    res <- toupper(digest::digest( fpath_bin, algo = "sha1", file = T))
    if( res == SHA )
    {
      cat("OK\n")
      bChecked <- T
    }
    else
    {
      cat(paste("NOK\nChecksums don't match\nXML key:", SHA, "\nBinary file key:", res,"\n"))
      .controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
    }
  }
  if( !is.null(MD5))
  {
    cat("Checking binary data checksum using MD5 key... ")
    res <- toupper(digest::digest( fpath_bin, algo = "md5", file = T))
    if( res == MD5 )
    {
      cat("OK\n")
      bChecked <- T
    }
    else
    {
      cat(paste("NOK\nChecksums don't match\nXML key:", MD5, "\nBinary file key:", res,"\n"))
      .controlled_loadAbort("ERROR: possible data corruption\n", close_signal)
    }
  }
  if(!bChecked)
  {
    cat("WARNING: No checksum available to test data integrity\n")
  }

  #Obtain data size per bin
  parseXmlArrayDataType <- function(xmlArrNode)
  {
    no_compression <- F
    units <- NULL
    dataType <- NULL
    xmlSubDataType <- XML::xmlChildren(xmlArrNode)
    for( j in 1:(XML::xmlSize(xmlSubDataType)))
    {
      name <- XML::xmlGetAttr(xmlSubDataType[[j]], "name")
      id <- XML::xmlGetAttr(xmlSubDataType[[j]], "accession")
      value <- XML::xmlGetAttr(xmlSubDataType[[j]], "value")

      if(name == "no compression" && id == "MS:1000576") #test compression
      {
        no_compression <- T
      }

      if(name == "m/z array" && id == "MS:1000514") #m/z array units
      {
        units <-  XML::xmlGetAttr(xmlSubDataType[[j]], "unitName")
      }

      if(name == "intensity array" && id == "MS:1000515") #intensity array units
      {
        units <-  XML::xmlGetAttr(xmlSubDataType[[j]], "unitName")
      }

      #8 bits integer data array
      #TODO This is not encoded in mzML CV obo file but is specified as valid in imzML spec!

      #16 bits integer data array
      #TODO This is not encoded in mzML CV obo file but is specified as valid in imzML spec!

      #32 bits integer data array
      if(name == "32-bit integer" && id == "MS:1000519")
      {
        dataType <- "int"
      }

      #64 bits integer data array
      if(name == "64-bit integer" && id == "MS:1000522")
      {
        dataType <- "long"
      }

      #32 bits float data array
      if(name == "32-bit float" && id == "MS:1000521")
      {
        dataType <- "float"
      }

      #64 bits double data array
      if(name == "64-bit float" && id == "MS:1000523")
      {
        dataType <- "double"
      }

    }
    return(list( no_compression = no_compression, units = units, dataType = dataType ))
  }

  mzArrayDesc <- NULL
  intArrayDesc <- NULL
  xmlDataType <- XML::xmlChildren( xmlch$referenceableParamGroupList)
  for ( i in 1:(XML::xmlSize(xmlDataType)))
  {
    if ( XML::xmlGetAttr(xmlDataType[[i]], "id") == "mzArray")
    {
      mzArrayDesc <- parseXmlArrayDataType(xmlDataType[[i]])
    }

    if ( XML::xmlGetAttr(xmlDataType[[i]], "id") == "intensityArray")
    {
      intArrayDesc <- parseXmlArrayDataType(xmlDataType[[i]])
    }
  }

  #Check for emty data
  err_string <- character(0)
  if( is.null(mzArrayDesc$units ))
  {
    err_string <- paste(err_string, "No m/z array units declared\n", sep = "")
  }
  if( is.null(intArrayDesc$units ))
  {
    cat("\nWarning: No intensity array units declared. Assuming number of counts as default unit.\n") #Just display a wrining because some implementation do not provide data units for intensity array
    intArrayDesc$units <- "number of counts"
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

  if (mzArrayDesc$units != "m/z")
  {
    .controlled_loadAbort("m/z Array is not in m/z format, rMSI only supports m/z units\n", close_signal)
  }

  if( intArrayDesc$units != "number of counts")
  {
    .controlled_loadAbort("intensity Arrays are not in number of counts units, rMSI does not support it\n", close_signal)
  }

  #Obtain pixel resolution (I'm reading directly the accession becasue I've found some mismatch beetween obo and imzML example)
  pixel_size_um <- NULL
  xmlScanSettings <- XML::xmlChildren(xmlch$scanSettingsList)
  for( i in 1:(XML::xmlSize(xmlScanSettings$scanSettings)))
  {
    accession <- XML::xmlGetAttr(xmlScanSettings$scanSettings[[i]], "accession")
    if( accession == "IMS:1000046")
    {
      pixel_size_um <- as.numeric(XML::xmlGetAttr(xmlScanSettings$scanSettings[[i]], "value"))
    }
  }

  #Obtain the RUN data
  imgData <- xmlParseSpectra(XML::xmlChildren(xmlch$run), fun_progress = fun_progress)

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
#'
#' @return a named data frame containing pixels positions in image (x, y) and spectra poisitions in binary file (offsets).
#'
xmlParseSpectra <- function (xmlRunNode, fun_progress = NULL)
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
    imgData[i, ] <- unlist(xmlParseSpectrum(xmlRunNode$spectrumList[[i]]))
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
  cat(paste("imzML XML parse time:",round(pt["elapsed"], digits = 1),"seconds\n"))
  return( imgData )
}

#' xmlParseSpectrum.
#'
#' Parses an XML node in imzML format containing a single spectrum. All meaning data is returned as named list.
#' Internal use only.
#'
#' @param xmlSingleSpectrum a xml node containing data for a single spectrum.
#'
#' @return a named list containing all parameters of the spectrum.
#'
xmlParseSpectrum <- function( xmlSingleSpectrum )
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
    if(XML::xmlGetAttr(ttscanPosXY[[i]], "name") ==  "position x")
    {
      spectrumInfo$x <- as.numeric(XML::xmlGetAttr(ttscanPosXY[[i]], "value"))
    }
    if(XML::xmlGetAttr(ttscanPosXY[[i]], "name") ==  "position y")
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
      name <- XML::xmlGetAttr(subNode[[j]], "name")
      val <- as.numeric(XML::xmlGetAttr(subNode[[j]], "value"))

      if (subNodeRef == "mzArray")
      {
        if( name == "external array length")
        {
          spectrumInfo$mzLength <- val
        }
        if( name == "external offset")
        {
          spectrumInfo$mzOffset <- val
        }
      }

      if (subNodeRef == "intensityArray")
      {
        if( name == "external array length")
        {
          spectrumInfo$intLength <- val
        }
        if( name == "external offset")
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

