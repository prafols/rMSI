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
  pt <- Sys.time()
  CalibrationWindowElapsedTime <- 0 #Keep track of the elapsed time during the calibration GUI
  
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
  
  #TODO peak-picking and binning
  
  elap <- Sys.time() - pt
  elap <- elap - CalibrationWindowElapsedTime #Substract the calbration GUI elapsed time
  cat("Total used processing time:\n")
  print(elap)
  
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

  #Run the preprocessing
  if( proc_params$preprocessing$smoothing$enable ||
      proc_params$preprocessing$alignment$enable ||
      proc_params$preprocessing$massCalibration  ) #TODO add basline condition here
  {
    #Calc new UUID's'
    uuids_new <- c()
    out_imzML_fnames <- c()
    for( i in 1:length(img_lst))
    {
      uuids_new <- c(uuids_new, uuid_timebased())
      out_imzML_fnames <- c(out_imzML_fnames, paste0(img_lst[[i]]$name, "-proc"))
    }
    
    result <-  CRunPreProcessing( img_lst, numOfThreads, memoryPerThreadMB, 
                                  proc_params$preprocessing, refSpc, 
                                  uuids_new, proc_params$outputpath, out_imzML_fnames, 
                                  common_mass)
    
    #Calculate the calibration model
    pt <- Sys.time() #do not take into account the user-time during the calibration GUI!
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
    calibrationElapsedTime <- Sys.time() - pt
    
    
    #Set the preprocessed data using resulting offsets and original data info
    img_lst_proc <- list()
    for( i in 1:length(img_lst))
    {
      img_lst_proc[[i]] <- img_lst[[i]] #copy the original data
      img_lst_proc[[i]]$mass <- common_mass
      img_lst_proc[[i]]$name <- paste0(img_lst_proc[[i]]$name, "-proc")
      img_lst_proc[[i]]$mean <- result$AverageSpectra[[i]]
      img_lst_proc[[i]]$base <- result$BaseSpectra[[i]]
      img_lst_proc[[i]]$data$path <- proc_params$outputpath
      img_lst_proc[[i]]$data$rMSIXBin$uuid <- uuid_timebased()
      img_lst_proc[[i]]$data$rMSIXBin$file <- img_lst_proc[[i]]$name 
      img_lst_proc[[i]]$data$imzML$uuid <- uuids_new[i]
      img_lst_proc[[i]]$data$imzML$file <- out_imzML_fnames[i]
      img_lst_proc[[i]]$data$imzML$SHA <- NULL #Remove posible SHA checksum since it must be recalculated
      img_lst_proc[[i]]$data$imzML$MD5 <- NULL #Remove posible MD5 checksum since it must be recalculated
      img_lst_proc[[i]]$data$imzML$run[, -c(1,2)] <- result$Offsets[[i]] 
      
      #Apply mass recalibration here to all images
      if(proc_params$preprocessing$massCalibration)
      {
        img_lst_proc[[i]] <- rMSI::applyMassCalibrationImage(img_lst_proc[[i]], calModel$model)
      }
      
      #The XML part is stored after the mass re-calibration since the mass axis overwrittening process will change the results of the checksums
      cat(paste0("Calculating MD5 checksum for image ", img_lst_proc[[i]]$name, "...\n"))
      img_lst_proc[[i]]$data$imzML$MD5 <- toupper(digest::digest( file.path( img_lst_proc[[1]]$data$path, paste0(img_lst_proc[[i]]$data$imzML$file, ".ibd")),
                                                                  algo = "md5",
                                                                  file = T))
      
      #Store the xml part
      cat(paste0("Writing the .imzML file for image ", img_lst_proc[[i]]$name, "...\n"))
      if(!CimzMLStore( path.expand(file.path( img_lst_proc[[1]]$data$path, paste0(img_lst_proc[[i]]$data$imzML$file, ".imzML"))), 
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
    return( list( processed_data = img_lst_proc, LagLow = result$LagLow, LagHigh = result$LagHigh, CalibrationElapsedTime = calibrationElapsedTime ))
  }
  else
  {
    cat("\nPre-processing Bypassed\n")
    return( list( raw_data = img_lst))
  }
}
