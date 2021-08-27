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

#ifndef MT_ALIGN_H
  #define MT_ALIGN_H
#include <Rcpp.h>

#include "smoothing.h"
#include "labelfreealign.h"
#include "threadingmsiproc.h"
#include "noiseestimation.h"

class MTPreProcessing : public ThreadingMsiProc 
{
  public:
    
    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    // preProcessingParams: An R reference class with the pre-processing parameters.
    // mass: a numeric vector with the common mass axis
    // reference: a reference spectrum for the alignment
    // uuid: a string vector with new UUID for the output imzML files
    // outputImzMLPath: an existing target path to store imzML files with the processed data
    // outputImzMLfnames: a string vector with the file names for the output imzML files
    // forceDataResampling: if data in continuous mode and forceDataResampling is set to true interpolation will be used.
    // bitDepthReductionNoiseWindows: The noise estimation windows used by the bitdepth reduction
    MTPreProcessing(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB,
                    Rcpp::Reference preProcessingParams, Rcpp::NumericVector reference,
                    Rcpp::StringVector uuid, Rcpp::String outputImzMLPath, Rcpp::StringVector outputImzMLfnames, 
                    bool forceDataResampling,
                    int bitDepthReductionNoiseWindows = 16);
    ~MTPreProcessing();
 
    //Exectue a full imatge processing using threaded methods and returns the used shifts in the first iteration
    Rcpp::List Run(); 
    
    //Single spectrum bit depth reduction
    void BitDepthReduction(double *data, int dataLength, int noiseModelThreadSlot);

  private:
    bool bEnableSmoothing; //Set to true if smoothing must be performed
    bool bEnableAlignment; //Set to true if alignment must be performed

    Smoothing **smoothObj;
    LabelFreeAlign **alngObj;
    LabelFreeAlign::TLags *mLags; //A place to store alignment lags
    
    //Bit depth reduction data
    NoiseEstimation **noiseModel;
    int NoiseWinSize;
    unsigned long long maskLUT_double[52];
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
    
};
#endif
