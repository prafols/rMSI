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

#' DetectPeaks'
#'
#' @param mass a vector containing the mass axis.
#' @param intensity a vector containinf the spectrum intensities.
#' @param SNR the minimum signal to noise ratio of retained peaks
#' @param WinSize the used windows size for peak detection
#' @param OverSampling the used oversampling value for interpolating the peak shape and improve mass and area calculation.
#'
#' @return a list containing mass, intensity, SNR, area and the binSize arround peak fields of detected peaks.
#' @export
#'
DetectPeaks <- function(mass, intensity, SNR = 5, WinSize = 20, OverSampling = 10)
{
  pm <- DetectPeaks_C(mass, intensity,SNR, WinSize, OverSampling)
  peaks <- list()
  peaks$mass <- pm["mass", ]
  peaks$intensity <- pm["intensity", ]
  peaks$SNR <- pm["SNR", ]
  peaks$area <- pm["area", ]
  peaks$binSize <- pm["binSize", ]
  return(peaks)
}
