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
#include "imzMLBin.h"

class PeakBinning
{
public:
  
  //Constructor arguments:
  // preProcessingParams: An R reference class with the pre-processing parameters.
  // numberOfThreads: Total number of threads to use during processing
  PeakBinning(Rcpp::Reference preProcessingParams, int numberOfThreads);
  ~PeakBinning();
  
  //Appends an image to be processed.
  // - imzMLDescriptor: An R list object describing the peakList file. Same format as outputed by the imzML parser
  void appedImageData(Rcpp::List imzMLDescriptor);
  
  //Perfomr peak binning, this is mono-thread implemented. //TODO work on the multithreaded implementation using C++11 tasks and paralelizin mass searches
  Rcpp::List BinPeaks(); 

private:
  double binFilter;
  double tolerance;
  bool tolerance_in_ppm; //If true the binning tolerance is specified in  ppm, if false then the number of datapoints per peak is used instead
  int totalNumOfPixels;
  
  std::vector<ImzMLBinRead*> imzMLReaders;  //Pointers to multiple imzMLReadrs initialized with openIbd = false to avoid exiding the maximum open files.
  
  typedef struct
  {
    double intensity;
    double SNR;
    double area;
  }TBin;
};
#endif
