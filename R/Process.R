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
                          memoryPerThreadMB = 200 )
{
  pt <- Sys.time()
  
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
    #TODO this is only needed when alignment is enabled
    common_mass <- img_lst[[1]]$mass
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
      
      #Replace each image mass axis with the common
      for( i in 1:length(img_lst))
      {
        img_lst[[i]]$mass <- common_mass
      }
    }
    
    #Calculate the internal reference for alignment
    AverageSpectrum <- COverallAverageSpectrum(img_lst, numOfThreads, memoryPerThreadMB) #TODO  it is crashing here!
    
  }

  #Calculate reference spectrum for label free alignment
  if(proc_params$preprocessing$alignment$enable)
  {
    #TODO error here! if mergeprocessing is disables then the reference cannot be calcualted!
    refSpc <- InternalReferenceSpectrumMultipleDatasets(img_lst, AverageSpectrum)
    cat(paste0("Pixel with ID ", refSpc$ID, " from image indexed as ", refSpc$imgIndex, " (", img_lst[[ refSpc$imgIndex]]$name, ") selected as internal reference.\n"))
    refSpc <- refSpc$spectrum
  }
  else
  {
    refSpc <- rep(0.0, length(img_lst[[1]]$mass)) #I need to supply a reference spectrum even if alignment is not enables, so just feed it with zeros
  }
  
  #TODO run the preprocssing!
  procData <- RunPreProcessing( img_lst, numOfThreads, memoryPerThreadMB, 
                                proc_params$preprocessing, refSpc)
  
  #TODO return the processed data instead of the average!
  return(AverageSpectrum)
  
  
  #return(img_lst)
}

