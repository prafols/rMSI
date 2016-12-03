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

#' NormalizeTIC: Calculates the TIC normalizatin of each pixel as the sum of all intensities.
#'
#' @param img the rMSI image object.
#' @param remove_empty_pixels boolean detailing if pixels detected to not contain data must be removed from normalization (smaller than mean-*sd).
#'
#' @return  a rMSI image containing the normalizations$TIC field or TICne if remove_empty_pixels is true.
#' @export
#'
NormalizeTIC <- function(img, remove_empty_pixels = FALSE)
{
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)
  TICs <- c()
  for( i in 1:length(img$data))
  {
    TICs <- c(TICs, rowSums(img$data[[i]][,]))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  if(remove_empty_pixels)
  {
    minallowedTIC <- mean(TICs) - sd(TICs)
    TICs[which(TICs < minallowedTIC)] <- Inf #Remove pixels divinding them by infinite
    img <- AppendNormalizationCoefs(img, "TICne", TICs)
  }
  else
  {
    img <- AppendNormalizationCoefs(img, "TIC", TICs)
  }
  return(img)
}

#' NormalizeMAX: Calculates the MAX normalizatin of each pixel as the maximum of all intensities.
#'
#' @param img the rMSI image object.
#' @param remove_empty_pixels boolean detailing if pixels detected to not contain data must be removed from normalization (smaller than mean-*sd).
#'
#' @return  a rMSI image containing the normalizations$MAX field or MAXne if remove_empty_pixels is true..
#' @export
#'
NormalizeMAX <- function(img, remove_empty_pixels = FALSE)
{
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)
  MAXs <- c()
  for( i in 1:length(img$data))
  {
    MAXs <- c(MAXs, apply(img$data[[i]][,], 1, max))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  if(remove_empty_pixels)
  {
    minallowedMAX <- mean(MAXs) - sd(MAXs)
    MAXs[which(MAXs < minallowedMAX)] <- Inf #Remove pixels divinding them by infinite
    img <- AppendNormalizationCoefs(img, "MAXne", MAXs)
  }
  else
  {
    img <- AppendNormalizationCoefs(img, "MAX", MAXs)
  }
  return(img)
}

#' NormalizeTargetPeaks: Calculates the normalization of each pixel as the sum of areas of specified peaks.
#' The spectra that not contains target mass peaks will be removed by normalizing them by Inf.
#'
#' @param img the rMSI image object.
#' @param mz_vector a vector containing m/z values of peaks to be used for normalization.
#' @param mz_window_size specifies the with of peak gaussian. Set to 0.1 Da by default which is OK for reflector TOF MALDI instrument
#' @param norm_name a name to be used for the normalization. Default is set to "SelPK".
#'
#' @return  a rMSI image containing the normalizations$<name> field
#' @export
#'
NormalizeTargetPeaks <- function(img, mz_vector, mz_window_size = 0.1, norm_name = "SelPK")
{
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)

  #Obtain columns numbers to be used for normalization
  iCols <- c()
  for( i in 1:length(mz_vector))
  {
    iBot <- which.min(abs(img$mass - (mz_vector[i] - 0.5*mz_window_size)))
    iTop <- which.min(abs(img$mass - (mz_vector[i] + 0.5*mz_window_size)))
    iCols <- c(iCols, iBot:iTop)
  }

  #Calc normalizations
  Norms <- c()
  for( i in 1:length(img$data))
  {
    Norms <- c(Norms, rowSums(img$data[[i]][,iCols]))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #Remove spectra that no contains target peaks as it is not normalizable.
  minallowedNorm <- mean(Norms) - sd(Norms)
  Norms[which(Norms < minallowedNorm)] <- Inf #Remove pixels divinding them by infinite
  img <- AppendNormalizationCoefs(img, norm_name, Norms)

  return(img)
}

#' NormalizeByAcqDegradation: Normalizes an rMSI image to compensate ionization source degradation during acquisition.
#'
#' @param img  the rMSI object to normalize.
#' @param winSize the window size use for smoothing (0 to 1).
#'
#' @return a rMSI image containing the normalizations$AcqTic field.
#' @export
#'
NormalizeByAcqDegradation <- function( img, winSize = 0.1 )
{
  #Calc all tics
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)
  TICs <- c()
  for( i in 1:length(img$data))
  {
    TICs <- c(TICs, rowSums(img$data[[i]][,]))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #Sort TICs according acquisition
  idxArray <- matrix( c(1:nrow(img$pos), img$pos[,"x"], img$pos[,"y"]), nrow = nrow(img$pos), ncol = 3, byrow = F )
  colnames(idxArray) <- c("id", "x", "y")
  idxArray <- idxArray[order(idxArray[,"y"], idxArray[,"x"]), ]
  TICs <- matrix( c( idxArray[,"id"], TICs[idxArray[,"id"]]), ncol = 2)
  colnames(TICs) <- c("id", "TIC")

  smTICs <- TICs
  winLen <- floor(0.5*nrow(TICs) * winSize)
  for( i in 1:nrow(TICs))
  {
    iStart <- i - winLen
    iStop <- i + winLen
    if( iStart < 1)
    {
      iStop <- iStop + (1 - iStart)
      iStart <- 1
    }
    if( iStop > nrow(TICs))
    {
      iStart <- iStart + (nrow(TICs) - iStop)
      iStop <- nrow(TICs)
    }

    discardLows <- mean(TICs[iStart:iStop, "TIC"]) - sd(TICs[iStart:iStop, "TIC"])
    dWind <- TICs[iStart:iStop, "TIC"]
    dWind <- dWind[which(dWind > discardLows)]
    smTICs[i, "TIC"] <- mean(dWind)
    #TODO pensar que fer si es NAN
  }

  #Order smTICs according Id's and append it to image normalization
  smTICs <- smTICs[order(smTICs[,"id"]), ]
  img <- AppendNormalizationCoefs(img, "AcqTic", smTICs[,"TIC"])
  return(img)
}

#' NormalizeByAcqDegradationOnTargetPeaks: Calculates the normalization of each pixel as the sum of areas of specified peaks.
#' Each normalization coeficient is smoothed by its neighbours to compensate for ionization degradation during acqusition.
#'
#' @param img the rMSI image object.
#' @param mz_vector a vector containing m/z values of peaks to be used for normalization.
#' @param mz_window_size specifies the with of peak gaussian. Set to 0.1 Da by default which is OK for reflector TOF MALDI instrument
#' @param norm_name a name to be used for the normalization. Default is set to "SelPK".
#' @param winSize the window size use for smoothing (0 to 1).
#'
#' @return  a rMSI image containing the normalizations$<name> field
#' @export
#'
NormalizeByAcqDegradationOnTargetPeaks <- function(img, mz_vector, mz_window_size = 0.1, norm_name = "AcqSelPK", winSize = 0.1)
{
  pb<-txtProgressBar(min = 0, max = length(img$data), style = 3 )
  setTxtProgressBar(pb, 0)

  #Obtain columns numbers to be used for normalization
  iCols <- c()
  for( i in 1:length(mz_vector))
  {
    iBot <- which.min(abs(img$mass - (mz_vector[i] - 0.5*mz_window_size)))
    iTop <- which.min(abs(img$mass - (mz_vector[i] + 0.5*mz_window_size)))
    iCols <- c(iCols, iBot:iTop)
  }

  #Calc normalizations
  Norms <- c()
  for( i in 1:length(img$data))
  {
    Norms <- c(Norms, rowSums(img$data[[i]][,iCols]))
    setTxtProgressBar(pb, i)
  }
  close(pb)

  #Sort Norms according acquisition
  idxArray <- matrix( c(1:nrow(img$pos), img$pos[,"x"], img$pos[,"y"]), nrow = nrow(img$pos), ncol = 3, byrow = F )
  colnames(idxArray) <- c("id", "x", "y")
  idxArray <- idxArray[order(idxArray[,"y"], idxArray[,"x"]), ]
  Norms <- matrix( c( idxArray[,"id"], Norms[idxArray[,"id"]]), ncol = 2)
  colnames(Norms) <- c("id", "Norm")

  smNorm <- Norms
  winLen <- floor(0.5*nrow(Norms) * winSize)
  for( i in 1:nrow(Norms))
  {
    iStart <- i - winLen
    iStop <- i + winLen
    if( iStart < 1)
    {
      iStop <- iStop + (1 - iStart)
      iStart <- 1
    }
    if( iStop > nrow(Norms))
    {
      iStart <- iStart + (nrow(Norms) - iStop)
      iStop <- nrow(Norms)
    }

    discardLows <- mean(Norms[iStart:iStop, "Norm"]) - sd(Norms[iStart:iStop, "Norm"])
    dWind <- Norms[iStart:iStop, "Norm"]
    dWind <- dWind[which(dWind > discardLows)]
    smNorm[i, "Norm"] <- mean(dWind)
    #TODO pensar que fer si es NAN
  }

  #Order smNorms according Id's and append it to image normalization
  smNorm <- smNorm[order(smNorm[,"id"]), ]
  img <- AppendNormalizationCoefs(img, norm_name, smNorm[,"Norm"])

  return(img)
}

#' AppendNormalizationCoefs: Appends a new normalization coefinients to rMSI object.
#'
#' @param img the rMSI object to append normalization coeficients.
#' @param normName the given name of the normalization.
#' @param normCoefs the normalizations coeficients.
#'
#' @return  a rMSI image containing the new normalization field.
#' @export
#'
AppendNormalizationCoefs <- function(img, normName, normCoefs)
{
  if(is.null(img$normalizations))
  {
    img$normalizations <- list()
  }

  img$normalizations[[normName]] <- normCoefs

  return(img)
}
