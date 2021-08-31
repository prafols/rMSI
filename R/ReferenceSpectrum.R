#########################################################################
#     rMSIproc - R package for MSI data processing
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

#' InternalReferenceSpectrum.
#' 
#' Calculates the dataset reference spectrum to use in label free alignment.
#' The reference spectrum is the spectrum in the dataset with the best correlation to average spectrum and with high TIC.
#'
#' @param img MSI image to calculate internal reference spectrum.
#' @param reference the spectrum to use as reference for correlations.
#' @param commonMassAxis the common mass axis for spectra interpolation.
#'
#' @return a list with the intensity vector corresponding to the reference spectrum, the score and the pixel ID selected as reference.
#'
InternalReferenceSpectrum <- function(img, reference, commonMassAxis = img$mass)
{
  pb <- txtProgressBar(min = 0,  max = nrow(img$pos), style = 3)
  maxScore <- 0
  maxId <- 0
  selID <- NA
  for( id in 1:nrow(img$pos))
  {
    spc <- rMSI::loadImgChunkFromIds(img, id, commonMassAxis)[1,]
    if(var(spc) > 0)
    {
      score <- cor(reference, spc )
      if( score > maxScore )
      {
        maxScore <- score
        maxId <- id
        selID <- id
      }
    }
    setTxtProgressBar(pb, id)
  }
  
  close(pb)
  return( list(spectrum = rMSI::loadImgChunkFromIds(img, maxId, commonMassAxis)[1,], score = maxScore, ID = selID ))
}

#' InternalReferenceSpectrumMultipleDatasets.
#' 
#' Calculates the dataset reference spectrum to use in label free alignment.
#' The reference spectrum is the spectrum in the dataset with the best correlation to average spectrum and with high TIC.
#'
#' @param img_list a list of various rMSI objects.
#' @param reference the spectrum to use as reference for correlations.
#' @param commonMasAxis the common mass axis for spectra interpolation.
#'
#' @return a list with the intensity vector corresponding to the reference spectrum, the score, the image index which contains the refernce spectrum and the pixel ID selected as reference.
#' 
InternalReferenceSpectrumMultipleDatasets <- function(img_list, reference, commonMasAxis)
{
  #Calculate correlations
  bestRefs <- list()
  for( i in 1:length(img_list) )
  {
    cat(paste0("Calculating internal reference spectrum ", i, "/",length(img_list),"...\n"))
    bestRefs[[i]] <- InternalReferenceSpectrum(img_list[[i]], reference, commonMasAxis)
  }
  
  #Select the best correlation in all datasets
  bestID <- which.max(unlist(lapply(bestRefs, function(x){ return(x$score) })))

  return(  list(spectrum = bestRefs[[bestID]]$spectrum, score = bestRefs[[bestID]]$score, imgIndex = bestID, ID =  bestRefs[[bestID]]$ID ) )
}
