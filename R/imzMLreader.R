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

#' fixImzMLDuplicates.
#' Delete duplicates and possible zero-drop errors (fixing Bruker's bugs in imzML).
#' 
#' @param mass spectrum mass axis.
#' @param intensity spectrum intensity.
#'
#' @return a list with non-duplicated mass and intensity vectors.
#'
fixImzMLDuplicates <- function(mass, intensity)
{
  idup <- which(duplicated(mass))
  if(length(idup) > 0)
  {
    for( i in 1:length(idup))
    {
      intensity[idup[i] - 1] <- max( intensity[idup[i]], intensity[idup[i] - 1])
    }
    mass <- mass[-idup]
    intensity <- intensity[-idup]
  }
  return(list(mass=mass, intensity=intensity))
}

#' import_imzML.
#'
#' @param imzML_File full path to .imzML file 
#' @param ibd_File path to the binary file (default the same as imzML file but with .ibd extension)
#' @param fun_progress This is a callback function to update the progress of loading data. See details for more information.
#' @param fun_text This is a callback function to update the label widget of loading data. See details for more information.
#' @param close_signal function to be called if loading process is abored.
#' @param verifyChecksum if the binary file checksum must be verified, it is disabled by default for convenice with really big files.
#' @param createImgStream true if the ion image stream must be created.
#' @param subImg_rename alternative image name, new rMSI files will be created with the given name.
#' @param subImg_Coords a Complex vector with the motors coordinates to be included in the rMSI data.
#' @param convertProcessed2Continuous if true (the default) an imzML file in processed mode will be converted to a continuous mode.
#'
#'  Imports an imzML image to an rMSI data object.
#'  It is recomanded to use rMSI::LoadMsiData directly instead of this function.
#'
#' @return an rMSI data object.
#' @export
#' 
import_imzML <- function(imzML_File, ibd_File =  paste(sub("\\.[^.]*$", "", imzML_File), ".ibd", sep = "" ),
                         fun_progress = NULL, fun_text = NULL, close_signal = NULL, 
                         verifyChecksum = F, createImgStream = T,
                         subImg_rename = NULL, subImg_Coords = NULL, convertProcessed2Continuous = T)
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
    
    #Provide a new name for the sub image if it has not been provided
    if(is.null(subImg_rename))
    {
      subImg_rename <- paste0(basename(imzML_File), "_subset")
    }
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

  #4- Obtain de m/z axis (for continuous mode just read the mass axis, for processed mode calculate a common mass axis)
  pt<-proc.time()
  dataPointEncoding_Mz  <- dataPointBinaryEncoding(xmlRes$mz_dataType)
  dataPointEncoding_Int <- dataPointBinaryEncoding(xmlRes$int_dataType)
  bNoNeed2Resample <- T #Start assuming there is no need to resample the data
  
  if(xmlRes$continuous_mode)
  {
    mzAxis <- readBin(bincon, dataPointEncoding_Mz$dataType, xmlRes$run_data[1, "mzLength"], size = dataPointEncoding_Mz$bytes, signed = T)
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
      
      icurr <- 1
      bLoad <- TRUE
      MergedSpc <- list()
      
      #Read only the first mass axis to compare if others are identical, this is the case for Bruker FTICR
      seek(bincon, rw = "read", where = xmlRes$run_data[1, "mzOffset"] )
      firstMassAxis <- unique(readBin(bincon, dataPointEncoding_Mz$dataType, xmlRes$run_data[1, "mzLength"], size = dataPointEncoding_Mz$bytes, signed = T))
      
      while(TRUE)
      {
        if(bLoad)
        {
          #Read mass axis for the current spectrum 
          seek(bincon, rw = "read", where = xmlRes$run_data[icurr, "mzOffset"] )
          mzdd <- readBin(bincon, dataPointEncoding_Mz$dataType, xmlRes$run_data[icurr, "mzLength"], size = dataPointEncoding_Mz$bytes, signed = T)
          
          #Read intensity of current spectrum
          seek(bincon, rw = "read", where = xmlRes$run_data[icurr, "intOffset"] )
          dd <- readBin(bincon, dataPointEncoding_Int$dataType, xmlRes$run_data[icurr, "intLength"], size = dataPointEncoding_Int$bytes, signed = T)
          
          #Fix duplicates and zero drops
          CurrSpectrumFixed <- fixImzMLDuplicates(mzdd, dd)
          
          #Get Bin size at peaks
          LoadMass <- CurrSpectrumFixed$mass
          bMassMerge <- T
          if(identical(firstMassAxis, LoadMass))
          {
            bMassMerge <- F
            LoadBins <- NULL
          }
          else
          {
            LoadBins <- CalcMassAxisBinSize( LoadMass, CurrSpectrumFixed$intensity)
          }
          
          bNoNeed2Resample <- bNoNeed2Resample & (!bMassMerge)
          
          #Shift register
          if(length(MergedSpc) > 0)
          {
            for(i in length(MergedSpc):1)
            {
              MergedSpc[[i+1]] <- MergedSpc[[i]]
            }
          }
          
          MergedSpc[[1]] <- list( level = 0, mass = LoadMass, bins =  LoadBins, merge = bMassMerge )
          icurr <- icurr + 1
          bLoad <- FALSE
          
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
        
        if(length(MergedSpc) > 1)
        {
          if(MergedSpc[[1]]$level == MergedSpc[[2]]$level || icurr > nrow(xmlRes$run_data))
          { 
            #Merge!
            if( MergedSpc[[2]]$merge)
            {
              mam <- rMSI::MergeMassAxis(MergedSpc[[1]]$mass, MergedSpc[[1]]$bins, MergedSpc[[2]]$mass, MergedSpc[[2]]$bins )
            }
            else
            {
              #Both mass axes are identical so there is no need to merge them, this is the case for Bruker FTICR data
              mam <- list(mass = MergedSpc[[1]]$mass, bins = MergedSpc[[1]]$bins )
            }
            MergedSpc[[1]] <- list( level = MergedSpc[[1]]$level + 1, mass = mam$mass, bins =  mam$bins, merge = MergedSpc[[2]]$merge )
            
            #Shift register
            if(length(MergedSpc) > 2)
            {
              for(i in 2:(length(MergedSpc)-1))
              {
                MergedSpc[[i]] <- MergedSpc[[i+1]]
              }
            }
            MergedSpc[[length(MergedSpc)]] <- NULL #Delete the last element
            
            #End Condition
            if(icurr > nrow(xmlRes$run_data) && length(MergedSpc) == 1)
            {
              break
            }
          }
          else
          {
            bLoad <- TRUE
          }
        }
        else
        {
          bLoad <- TRUE
        }
      }
      
      mzAxis <-  MergedSpc[[1]]$mass
      rm(MergedSpc)
      gc()
      pt<-proc.time() - pt
      cat(paste("\nMass axis calculation time:",round(pt["elapsed"], digits = 1),"seconds\n"))
      cat(paste("The re-sampled mass axis contains", length(mzAxis), "data points\n"))
    }
  }
  close(bincon)
 
  #5- Create the rMSIXBin
  if( is.null(subImg_rename))
  {
    subImg_rename <-  basename(imzML_File)
  }
  if(!is.null(fun_text))
  {
    fun_text("Loading data...")
  }
  img <- CreateEmptyImage(num_of_pixels = nrow(xmlRes$run_data), mass_axis = mzAxis, pixel_resolution = xmlRes$pixel_size_um,
                               img_name = subImg_rename,
                               rMSIXBin_path = path.expand(dirname(imzML_File)),
                               uuid = binUUID
                                )

  #Fill missing data
  img$data$rMSIXBin$file <- sub("\\.[^.]*$", "", basename(imzML_File))
  #img$data$rMSIXBin$uuid has been set by CreateEmptyImage
  #TODO fill all rMSIXBin data!
  
  img$data$imzML$file <- sub("\\.[^.]*$", "", basename(imzML_File))
  #img$data$imzML$uuid has been set by CreateEmptyImage
  if( xmlRes$SHA != "" )
  {
    img$data$imzML$SHA <- xmlRes$SHA
  }
  if( xmlRes$MD5 != "" )
  {
    img$data$imzML$MD5 <- xmlRes$MD5
  }
  img$data$imzML$continuous_mode <- bNoNeed2Resample #Use this bool instead of xmlRes data becaus Bruker FTICR use processed mode but the mass axis is actually replicated
  img$data$imzML$mz_dataType <- xmlRes$mz_dataType
  img$data$imzML$int_dataType <-xmlRes$int_dataType
  img$data$imzML$run <- xmlRes$run_data
  
  #Compute image size and arrange coords to avoid holes
  img$posMotors[, "x"] <- xmlRes$run_data[, "x"]
  img$posMotors[, "y"] <- xmlRes$run_data[, "y"]
  img$pos <- remap2ImageCoords(img$posMotors)
  img$size["x"] <- max(img$pos[,"x"])
  img$size["y"] <- max(img$pos[,"y"])
  
  
  #6- Read all spectra and creat the ImgStrem
  if(createImgStream && (xmlRes$continuous_mode || (!xmlRes$continuous_mode && convertProcessed2Continuous) ))
  {
    pt <- proc.time()
    cat("\nReading spectra from binary file...\n")
    
    #TODO call the function to create the imgStream here, such function still does not exsit
    #The rMSIXBin files (.XrMSI and .BrMSI) must be present at this point! so here, both file will be filled with the imgstream
    
    pt<-proc.time() - pt
    cat(paste("\nBinary file reading time:",round(pt["elapsed"], digits = 1),"seconds\n\n"))
  
    #Compute average spectrum
    #TODO aixo va abans d crear el imgStream ja que caal guardarho en el rMSIXBin, utilitzar imzML directament per calcularho, 
    #TODO     crec k pots passar del Multithreading pq es molt rapid el calcul i lectura disc es sequencial, 
    #TODO     la propia classe imzMLBin pot tenir un metode per calcul espectre mig (aixi ho pot cridar el rMSIXBin directament).
    #img$mean <- AverageSpectrum(img) #TODO si bull nomes peak list aixo no s'ha de cridar
  }
  
  #7- Just reading the peak lists
  if(!xmlRes$continuous_mode && !convertProcessed2Continuous)
  {
    pt <- proc.time()
    cat("\nReading peak lists from binary file...\n")
    img <- read_imzML_peaklists(img)
    pt<-proc.time() - pt
    cat(paste("\nBinary file reading time:",round(pt["elapsed"], digits = 1),"seconds\n\n"))
  }

  #8- And it's done, just return de rMSI object
  if(!is.null(pb))
  {
    close(pb)
  }
  gc()
  if(!is.null(fun_text))
  {
    fun_text("Done")
  }
  return(img)
}

#' read_imzML_peaklists.
#' 
#' Reads the peak list in an imzML file containing data in processed mode.
#' Use this function carefully since all peak lists are loaded in memory.
#'
#' @param img a rMSIObj with valid imzML data in processed mode.
#'
#' @return a rMSIObj with the peak list in the img$data$peaklist field.
#'
read_imzML_peaklists <- function(img)
{
  if(img$data$imzML$continuous_mode)
  {
    stop("Error: read_imzML_peaklists() is only available with imzML in processed mode.\n")
  }
    
  dataPointEncoding_Int <- dataPointBinaryEncoding(img$data$imzML$int_dataType)
  dataPointEncoding_Mz  <- dataPointBinaryEncoding(img$data$imzML$mz_dataType)
  bincon <- file(description = file.path(img$data$path, paste0(img$data$imzML$file, ".ibd")), open = "rb")
  
  ppStep<-100/nrow(img$data$imzML$run)
  pp<-0
  pb<-txtProgressBar(min = 0, max = 100, style = 3 )
  setTxtProgressBar(pb, pp)
  
  for(i in 1:nrow(img$data$imzML$run))
  {
    #Read intensity of current spectrum
    seek(bincon, rw = "read", where = img$data$imzML$run[i, "intOffset"] )
    dd <- readBin(bincon, dataPointEncoding_Int$dataType, img$data$imzML$run[i, "intLength"], size = dataPointEncoding_Int$bytes, signed = T)
    
    #Read mass axis for the current spectrum 
    seek(bincon, rw = "read", where = img$data$imzML$run[i, "mzOffset"] )
    mzdd <- readBin(bincon, dataPointEncoding_Mz$dataType, img$data$imzML$run[i, "mzLength"], size = dataPointEncoding_Mz$bytes, signed = T) 
    
    img$data$peaklist[[i]]<-list(mass = mzdd, intensity = dd)

    #Update progress bar
    pp<-pp+ppStep
    setTxtProgressBar(pb, pp)
  }
  close(pb)
  close(bincon)
  return(img)
}

#' dataPointBinaryEncoding
#' 
#' Get the number of bytes and data type that must be readed in binary streams for each encoded data point.
#'
#' @param str_dataType a string with the C-style data type: int, long, float and double.
#'
#' @return a list containing an integer with the number of bytes to read for each data point and the R data type.
#'
dataPointBinaryEncoding <- function(str_dataType)
{
  result <- list( bytes = NULL,  dataType = NULL)
  if(str_dataType == "int" || str_dataType == "float")
  {
    result$bytes <- 4
  }
  if(str_dataType == "long" || str_dataType == "double")
  {
    result$bytes <- 8
  }
  
  if(str_dataType == "int" || str_dataType == "long")
  {
    result$dataType <- integer()
  }
  if(str_dataType == "float" || str_dataType == "double")
  {
    result$dataType <- numeric()
  }

  return(result)
}

##TODO revise and reuse the following code for data interpolation!!! ()
reuseThisFunctionCodeProperly<-function()
{
bCreaterMSIXBin <- T

  
  if(bCreaterMSIXBin)
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