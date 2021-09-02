/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2014 Pere Rafols Soler
 * 
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 * 
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 * 
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 **************************************************************************/

#ifndef PEAKBINNING_H
#define PEAKBINNING_H
#include <Rcpp.h>
#include "peakpicking.h"

class PeakBinning
{
public:
  //TODO replace PeakPicking::Peaks **peakList, int pixelCount with a data structure refencing the imzML peaklist
  PeakBinning(PeakPicking::Peaks **peakList, int pixelCount, double binTolerance, bool toleranceInppm, double binFilter);
  ~PeakBinning();
  
  //Perfomr peak binning, this is mono-thread implemented.
  Rcpp::List BinPeaks(); 

private:
  bool binSizeInppm;
  double binFilter;
  double tolerance;
  bool tolerance_in_ppm; //If true the binning tolerance is specified in  ppm, if false then the number of datapoints per peak is used instead
  int numOfPixels;
  
  //TODO add datastructure to handle peak list in imzML an Tic normalizations... maybe I can add the peak list info inside each rMSIobject...
  
  typedef struct
  {
    double intensity;
    double SNR;
    double area;
  }TBin;
};
#endif
