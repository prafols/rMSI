#########################################################################
#     rMSI - R package for MSI data handling and visualization
#     Copyright (C) 2021 Pere Rafols Soler
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

#' ProcessImages.
#' 
#' Process a single image or multiple images with the complete processing workflow.
#' 
#'
#' @param proc_params a ProcParams object containing the processing parameters
#' @param data_description a DataInfo object containing the imzML files paths and parsed ROIs.
#' @param verifyImzMLChecksums a boolean indicating whether imzML checksums must be validated or not (default is not since it takes a long time).
#' @param numOfThreads the number number of threads used to process the data.
#' @param memoryPerThreadMB maximum allowed memory by each thread. The total number of trehad will be two times numOfThreads, so the total memory usage will be: 2*numOfThreads*memoryPerThreadMB.
#' 
#' @return
#' @export
#'
#' @examples
ProcessImages <- function(proc_params,
                          data_description,
                          verifyImzMLChecksums = F,
                          numOfThreads = min(parallel::detectCores()/2, 6),
                          memoryPerThreadMB = 100 )
{
  if(class(proc_params) != "ProcParams")
  {
    stop("ERROR: proc_params argument must be an object of class \"ProcParams\". Use the rMSI::ProcessingParameters() function to create a valid proc_params\n")
  }
  
  if(class(data_description) != "DataInfo")
  {
    stop("ERROR: data_description argument must be an object of class \"DataInfo\". Use the rMSI::ImzMLDataDescription() function to create a valid data_description.\n")
  }
  
  pt <- proc.time()
  CalibrationWindowElapsedTime <- 0 #Keep track of the elapsed time during the calibration GUI
  
  #Start by parsing ROI information
  data_description$parseROIs() 
  
  #Check if data output path is set and create it
  if(length(proc_params$outputpath) == 0)
  {
    stop("ERROR: Output data path is empty. Use the setOutputPath() method to set it for the processing parameters object.\n")
  }
  dir.create(proc_params$outputpath, showWarnings = F, recursive = T)
  
  # Load the data (only imzML parsing)
  img_lst <- list()
  for( i in 1:data_description$getNumberOfImages())
  {
    cat(paste0("Loading image ", i , " of ", data_description$getNumberOfImages(), "\n" ))
    img_desc <- data_description$getImgPathPos(i)
    if(is.null(img_desc$ROIpos))
    {
      #Handle images without ROIs
      img_coords <- NULL 
    }
    else
    {
      #Handle images with ROIs
      img_coords <- complex(real = img_desc$ROIpos$x, imaginary = img_desc$ROIpos$y)  
    }
    img_lst[[i]] <- import_imzML(imzML_File = img_desc$imzML, 
                                 verifyChecksum = verifyImzMLChecksums, 
                                 subImg_rename =  img_desc$name,  
                                 subImg_Coords = img_coords,
                                 convertProcessed2Continuous = TRUE) #TODO consider giving the option to set this to false to process peaklists
  }
  
  # At this point, img_lst contains a list with all the images to process. If various images must be extracted form the same imzML file, then there will be an item for
  # each image(aka ROI) in the img_lst all of them will point to the same imzML but with different saptial coordinates.
  
  if(proc_params$getMergedProcessing())
  {
    #Merge processing, process all images at once
    result <- RunPreProcessing(proc_params,
                        img_lst,
                        numOfThreads,
                        memoryPerThreadMB)
    
    #Get the time elapsed during calibration GUI
    CalibrationWindowElapsedTime <- result$CalibrationElapsedTime 
  }
  else
  {
    #Image by image processing, non-merging
    result <- list()
    for(i in 1:length(img_lst))
    {
      cat(paste0("\nProcessing image ", i, " of ", length(img_lst), " \n"))
      result[[i]] <- RunPreProcessing(proc_params,
                          list(img_lst[[i]]),
                          numOfThreads,
                          memoryPerThreadMB)
      
      #Get the time elapsed during calibration GUI
      CalibrationWindowElapsedTime <- CalibrationWindowElapsedTime + result[[i]]$CalibrationElapsedTime
    }
  }
  
  #TODO Intensity Normalizations are only calculated automatically when the alginemnt is enbaled (to calculate the reference spectrum). 
  #TODO So, at the end reuse them according to the desired output: if rMSIXBin must be exported reuse or calculate, if not just forget about normalizations.
  #TODO think about returning alignment lags here. Currently, im returning it but this may be confusing for the end user. Also the Calibration time... I dont need any of these!
  
  #Display the used processing time
  elap <- proc.time() - pt - CalibrationWindowElapsedTime
  display_processing_time(elap, "Total data processing time")
  
  return(result)
}

#' RunPreProcessing
#' 
#' Process a single image or multiple images with the complete processing workflow.
#' 
#'
#' @param proc_params a ProcParams object containing the processing parameters
#' @param img_lst a rMSI objects lis to process.
#' @param numOfThreads the number number of threads used to process the data.
#' @param memoryPerThreadMB maximum allowed memory by each thread. The total number of trehad will be two times numOfThreads, so the total memory usage will be: 2*numOfThreads*memoryPerThreadMB.
#'
#' @return 
RunPreProcessing <- function(proc_params,
                                img_lst,
                                verifyImzMLChecksums = F,
                                numOfThreads = min(parallel::detectCores()/2, 6),
                                memoryPerThreadMB = 200 )
{
  # Check if the mass axis is the same
  common_mass <- img_lst[[1]]$mass
  identicalMassAxis <- TRUE
  if( length(img_lst) > 1)
  {
    for( i in 2:length(img_lst))
    {
      identicalMassAxis <- identicalMassAxis & identical(common_mass, img_lst[[i]]$mass)
    }
  }
  
  # Calculate the new common mass axis 
  if(!identicalMassAxis)
  {
    for( i in 2:length(img_lst))
    {
      massMergeRes <- MergeMassAxisAutoBinSize(common_mass, img_lst[[i]]$mass)
      if(massMergeRes$error)
      {
        stop("ERROR: The mass axis of the images to merge is not compatible because they do not share a common range.\n")
      }
      common_mass <- massMergeRes$mass
    }
  }
    
  if(proc_params$preprocessing$alignment$enable)
  {  
    #Calculate normalizations, TIC normalization is needd for internal reference calculation so, when alginemtn is used normalizations will be precalculated
    Normalizations <- CNormalizations(img_lst, numOfThreads, memoryPerThreadMB, common_mass)
    #Replace each image normalizations
    for( i in 1:length(img_lst))
    {
      img_lst[[i]]$normalizations <- Normalizations[[i]]
    }
    
    #Get the 25% and 75% quantiles of TIC norms
    allTICs <- unlist(lapply(Normalizations, function(x){ x$TIC }))
    TICquantiles <- quantile(allTICs)
    ticMin <- TICquantiles[2] # 25%
    ticMax <- TICquantiles[4] # 75%
    rm(TICquantiles)
    rm(allTICs)
    rm(Normalizations)
    
    #Calculate the internal reference for alignment
    AverageSpectrum <- COverallAverageSpectrum(img_lst, numOfThreads, memoryPerThreadMB, common_mass, ticMin, ticMax) 
    
    refSpc <- InternalReferenceSpectrumMultipleDatasets(img_lst, AverageSpectrum, common_mass)
    cat(paste0("Pixel with ID ", refSpc$ID, " from image indexed as ", refSpc$imgIndex, " (", img_lst[[ refSpc$imgIndex]]$name, ") selected as internal reference.\n"))
    refSpc <- refSpc$spectrum
    
    #TODO refSpc must be baseline corrected the same as the rest of the data
    
    if(proc_params$preprocessing$smoothing$enable) #Apply smoothing to the reference spectrum if needed
    {
      refSpc <- Smoothing_SavitzkyGolay(refSpc, proc_params$preprocessing$smoothing$kernelSize)  
    }
  }
  else
  {
    #I need to supply a reference spectrum even if alignment is not enabled, so just feed it with zeros
    refSpc <- rep(0.0, length(img_lst[[1]]$mass)) 
  }

  #Calc new UUID's for the processed imzML files
  if( proc_params$preprocessing$smoothing$enable ||
      proc_params$preprocessing$alignment$enable ||
      proc_params$preprocessing$massCalibration  ) #TODO add basline condition here
  {
    uuids_new <- c()
    out_imzML_fnames <- c()
    for( i in 1:length(img_lst))
    {
      uuids_new <- c(uuids_new, uuid_timebased())
      out_imzML_fnames <- c(out_imzML_fnames, paste0(img_lst[[i]]$name, "-proc"))
    }
  }
    
  #Run the preprocessing
  if( proc_params$preprocessing$smoothing$enable ||
      proc_params$preprocessing$alignment$enable ) #TODO add basline condition here
  {
    
    result <-  CRunPreProcessing( img_lst, numOfThreads, memoryPerThreadMB, 
                                  proc_params$preprocessing, refSpc, 
                                  uuids_new, proc_params$outputpath, out_imzML_fnames, 
                                  common_mass)
  }

  #Calculate the calibration model
  if(proc_params$preprocessing$massCalibration)
  {
    pt <- proc.time() #do not take into account the user-time during the calibration GUI!
    if(proc_params$preprocessing$massCalibration)
    {
      calModel <- CalibrationWindow(common_mass, refSpc) 
      #TODO the ref spectrum will be set to zero if alignment is disabled! check if zero and supply an average instead
      #TODO what about adding a menu in the calibration window to allow selecting from multiple calibration sources (the average of each image, the skyline etc..)
      
      if(is.null( calModel$model))
      {
        #Calibration aborted by user
        stop("Calibration aborted") #TODO improve the abort sequence: maybe removing already created ibd files... asking for confirmation... think about it
      }
    }
    calibrationElapsedTime <- proc.time() - pt
  }
  else
  {
    calibrationElapsedTime <- 0 
  }
  
  #Set the preprocessed data using resulting offsets and original data info
  if( proc_params$preprocessing$smoothing$enable ||
      proc_params$preprocessing$alignment$enable ||
      proc_params$preprocessing$massCalibration  )  #TODO add basline condition here
  {
    img_lst_proc <- list()
    for( i in 1:length(img_lst))
    {
      img_lst_proc[[i]] <- img_lst[[i]] #copy the original data
      img_lst_proc[[i]]$mass <- common_mass
      img_lst_proc[[i]]$name <- paste0(img_lst_proc[[i]]$name, "-proc")
      img_lst_proc[[i]]$data$path <- proc_params$outputpath
      img_lst_proc[[i]]$data$rMSIXBin$uuid <- uuid_timebased()
      img_lst_proc[[i]]$data$rMSIXBin$file <- img_lst_proc[[i]]$name 
      img_lst_proc[[i]]$data$imzML$uuid <- uuids_new[i]
      img_lst_proc[[i]]$data$imzML$file <- out_imzML_fnames[i]
      img_lst_proc[[i]]$data$imzML$SHA <- NULL #Remove posible SHA checksum since it must be recalculated
      img_lst_proc[[i]]$data$imzML$MD5 <- NULL #Remove posible MD5 checksum since it must be recalculated
      
      #Part only available if preprocessing enabled
      if(proc_params$preprocessing$smoothing$enable ||
         proc_params$preprocessing$alignment$enable)  #TODO add basline condition here
      {
        img_lst_proc[[i]]$mean <- result$AverageSpectra[[i]]
        img_lst_proc[[i]]$base <- result$BaseSpectra[[i]] 
        img_lst_proc[[i]]$data$imzML$run[, -c(1,2)] <- result$Offsets[[i]]
      }
      
      #Apply mass recalibration here to all images
      if(proc_params$preprocessing$massCalibration)
      {
        img_lst_proc[[i]] <- rMSI::applyMassCalibrationImage(img_lst_proc[[i]], calModel$model)
      }
      
      #The XML part is stored after the mass re-calibration since the mass axis overwrittening process will change the results of the checksums
      cat(paste0("Calculating MD5 checksum for image ", img_lst_proc[[i]]$name, "...\n"))
      img_lst_proc[[i]]$data$imzML$MD5 <- toupper(digest::digest( file.path( img_lst_proc[[i]]$data$path, paste0(img_lst_proc[[i]]$data$imzML$file, ".ibd")),
                                                                  algo = "md5",
                                                                  file = T))
      
      #Store the xml part
      cat(paste0("Writing the .imzML file for image ", img_lst_proc[[i]]$name, "...\n"))
      if(!CimzMLStore( path.expand(file.path( img_lst_proc[[i]]$data$path, paste0(img_lst_proc[[i]]$data$imzML$file, ".imzML"))), 
                       list( UUID = img_lst_proc[[i]]$data$imzML$uuid,
                             continuous_mode = img_lst_proc[[i]]$data$imzML$continuous_mode,
                             compression_mz = F,
                             compression_int = F,
                             MD5 = img_lst_proc[[i]]$data$imzML$MD5,
                             SHA = "",
                             mz_dataType = img_lst_proc[[i]]$data$imzML$mz_dataType,
                             int_dataType = img_lst_proc[[i]]$data$imzML$int_dataType,
                             pixel_size_um = img_lst_proc[[i]]$pixel_size_um,
                             run_data = img_lst_proc[[i]]$data$imzML$run
                       )))
      {
        stop(paste0("ERROR: imzML exported for image ", img_lst_proc[[i]]$name, " failed. Aborting...\n" ))
      }
    }
    
    cat("\nPre-processing completed\n")
  }
  else
  {
    cat("\nPre-processing Bypassed\n")
    img_lst_proc <- img_lst
    LagLow <- NULL
    LagHigh <- NULL
  }
  
  #Peak-picking 
  if(proc_params$preprocessing$peakpicking$enable)
  {
    #Calc new UUID's for the peaklists
    uuids_peakLists <- c()
    out_imzMLpeakLists_fnames <- c()
    for( i in 1:length(img_lst))
    {
      uuids_peakLists <- c(uuids_peakLists, uuid_timebased())
      out_imzMLpeakLists_fnames <- c(out_imzMLpeakLists_fnames, paste0(img_lst[[i]]$name, "-peaks"))
    }
    
    peakListsOffsets <- CRunPeakPicking(img_lst_proc, numOfThreads, memoryPerThreadMB, 
                                        proc_params$preprocessing, 
                                        uuids_peakLists, proc_params$outputpath, out_imzMLpeakLists_fnames, 
                                        common_mass)
    
    #Store peak lists imzML files and keep references to them in peaklists_lst
    peaklists_lst <- list()
    for( i in 1:length(img_lst_proc))
    {
      #Prepare the current peak list offset info
      currentRunData <- img_lst_proc[[i]]$data$imzML$run
      currentRunData[, -c(1,2)] <- peakListsOffsets[[i]]
      
      #The XML part is stored after the mass re-calibration since the mass axis overwrittening process will change the results of the checksums
      cat(paste0("Calculating MD5 checksum for the peaks lists of image ", img_lst_proc[[i]]$name, "...\n"))
      peaklistMD5 <- toupper(digest::digest( file.path( proc_params$outputpath, paste0(out_imzMLpeakLists_fnames[i], ".ibd")),
                                                                  algo = "md5",
                                                                  file = T))
      
      #Store peak list imzML info in a new list
      peaklists_lst[[i]] <- list( UUID = uuids_peakLists[i],
                                  continuous_mode = F,
                                  compression_mz = F,
                                  compression_int = F,
                                  MD5 = peaklistMD5,
                                  SHA = "",
                                  mz_dataType = "double",
                                  int_dataType = "double",
                                  pixel_size_um = img_lst_proc[[i]]$pixel_size_um,
                                  run_data = currentRunData,
                                  path = proc_params$outputpath,
                                  file = out_imzMLpeakLists_fnames[i],
                                  rMSIpeakList = T #Peak list created with rMSI format
                                  )
      
      #Store the xml part of the peak list specifing it is a peak list in rMSI format
      cat(paste0("Writing the .imzML file the peaks lists of image ", img_lst_proc[[i]]$name, "...\n"))
      if(!CimzMLStore( path.expand(file.path( peaklists_lst[[i]]$path , paste0(peaklists_lst[[i]]$file, ".imzML"))), 
                       peaklists_lst[[i]],
                        "rMSIpeakList"))
      {
        stop(paste0("ERROR: imzML exported for image ", img_lst_proc[[i]]$name, " failed. Aborting...\n" ))
      }
      
      #Add peak list description to each rMSI object
      img_lst_proc[[i]]$data$peaklist <- peaklists_lst[[i]]
    }
  }
  
  #Run the peakbining
  if(proc_params$preprocessing$peakbinning$enable)
  {
    peakMatrix <- CRunPeakBinning(peaklists_lst,  proc_params$preprocessing, numOfThreads)
    
    #Execute the fillpeaks after running the binning routine
    if(proc_params$preprocessing$peakbinning$fillpeaks) #TODO fill peaks must only be executed when spectral data is available.... force the fillpeaks param to FALSE when no spectral data is abailable (processing peaklists)
    {
      CRunFillPeaks(img_lst_proc, numOfThreads, memoryPerThreadMB, proc_params$preprocessing, common_mass, peakMatrix)
    }              
                   
    #Append normalizations to the peak matrix
    peakMatrix$normalizations <- img_lst_proc[[1]]$normalizations
    if(length(img_lst_proc) > 1)
    {
      for( i in 2:length(img_lst_proc))
      {
        peakMatrix$normalizations <- rbind(peakMatrix$normalizations, img_lst_proc[[i]]$normalizations)
      }
    }
    
    #Add a copy of img$pos to the peakMatrix
    mergedNames <- unlist(lapply(img_lst_proc, function(x){ return(x$name) }))
    mergedNumPixels <- unlist(lapply(img_lst_proc, function(x){ return(nrow(x$pos)) }))
    mergedPos <- matrix(ncol = 2, nrow = sum(mergedNumPixels))
    mergedMotors <- matrix(ncol = 2, nrow = sum(mergedNumPixels))
    mergedUUIDs <- unlist(lapply(img_lst_proc, function(x){ return(x$data$imzML$uuid) }))
    colnames(mergedPos) <- c("x", "y")
    colnames(mergedMotors) <- c("x", "y")
    istart <- 1
    for( i in 1:length(img_lst_proc))
    {
      istop <- istart + nrow(img_lst_proc[[i]]$pos) - 1
      mergedPos[ istart:istop , "x"] <- img_lst_proc[[i]]$pos[, "x"]
      mergedPos[ istart:istop , "y"] <- img_lst_proc[[i]]$pos[, "y"]
      
      mergedMotors[ istart:istop , "x"] <- img_lst_proc[[i]]$posMotors[, "x"]
      mergedMotors[ istart:istop , "y"] <- img_lst_proc[[i]]$posMotors[, "y"]

      istart <- istop + 1 
    }
    peakMatrix <- FormatPeakMatrix(peakMatrix, mergedPos,  mergedNumPixels, mergedNames, mergedUUIDs, mergedMotors) 
    
  }
  
  return( list( processed_data = img_lst_proc, LagLow = result$LagLow, LagHigh = result$LagHigh, CalibrationElapsedTime = calibrationElapsedTime, PeakMatrix = peakMatrix ))
}

#' FormatPeakMatrix.
#' Formats a C style peak matrix generated by MTPeakPicking::BinPeaks() to a rMSIprocPeakMatrix.
#'
#' @param cPeakMatrix a peak matrix with the same format as retured by MTPeakPicking::BinPeaks().
#' @param posMat a rMSI image pos matrix.
#' @param numPixels a vector including the number of pixels of each sample.
#' @param names a vector of strings with the name of each sample.
#' @param uuid a vector of img UUID to be also stored in peak matrices
#' @param posMotors a rMSI image original motros coordinates matrix.
#'
#' @return the formated matrix.
#'
FormatPeakMatrix <- function (cPeakMatrix, posMat, numPixels, names, uuid, posMotors)
{
  cPeakMatrix$pos <- posMat
  cPeakMatrix$numPixels <- numPixels
  cPeakMatrix$names <- names
  cPeakMatrix$uuid <- uuid
  cPeakMatrix$posMotors <- posMotors
  class(cPeakMatrix) <- "rMSIprocPeakMatrix"
  return(cPeakMatrix)
}

