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

#include <Rcpp.h>
#include <cmath>
#include "mtpreprocessing.h"
using namespace Rcpp;

MTPreProcessing::MTPreProcessing(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB,
                                 Rcpp::Reference preProcessingParams, Rcpp::NumericVector reference,
                                 Rcpp::StringVector uuid,  Rcpp::String outputImzMLPath, Rcpp::StringVector outputImzMLfnames, 
                                 Rcpp::NumericVector commonMassAxis,
                                 int bitDepthReductionNoiseWindows) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis, true, uuid, outputImzMLPath, outputImzMLfnames), 
  NoiseWinSize(bitDepthReductionNoiseWindows)
{
  //TODO add baseline params here!
  
  //Get the smoothing parameters
  Rcpp::Reference smoothingParams = preProcessingParams.field("smoothing");
  bEnableSmoothing = smoothingParams.field("enable");
  int smoothinKernelSize = smoothingParams.field("kernelSize");
  
  //Get the alignment params
  Rcpp::Reference alignmentParams = preProcessingParams.field("alignment");
  bEnableAlignment = alignmentParams.field("enable");
  bool bilinear = alignmentParams.field("bilinear");
  int alignIterations = alignmentParams.field("iterations");
  double maxShiftppm = alignmentParams.field("maxShiftppm");
  double lagRefLow = alignmentParams.field("refLow");
  double lagRefMid = alignmentParams.field("refMid");
  double lagRefHigh = alignmentParams.field("refHigh");
  int fftOverSampling = alignmentParams.field("overSampling");
  double winSizeRelative = alignmentParams.field("winSizeRelative");
  
  smoothObj = new Smoothing*[numOfThreadsDouble];
  alngObj = new LabelFreeAlign*[numOfThreadsDouble];
  noiseModel = new NoiseEstimation*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    smoothObj[i] = new Smoothing(smoothinKernelSize);
    
    alngObj[i] = new LabelFreeAlign(massAxis.begin(), reference.begin(), massAxis.length(), bilinear, 
                                    alignIterations, lagRefLow, lagRefMid, lagRefHigh,
                                    maxShiftppm,  fftOverSampling, winSizeRelative);
    
    noiseModel[i] = new NoiseEstimation(massAxis.length()); //Used by the bitdepth reduction
  }

  mLags  =  new LabelFreeAlign::TLags[numPixels]; 
  
  //Fill the bit depth reduction LUT with all possible resolutions from 1 to 52 bits
  maskLUT_double[0] = 0xFFF8000000000000;
  for( int i = 1; i < 52; i++)
  {
    maskLUT_double[i] = (maskLUT_double[i-1] >> 1) | 0x8000000000000000; 
  }
}

MTPreProcessing::~MTPreProcessing()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete smoothObj[i];
    delete alngObj[i];
    delete noiseModel[i];
  }
  delete[] smoothObj;
  delete[] alngObj;
  delete[] noiseModel;
  delete[] mLags;
}

List MTPreProcessing::Run()
{
  Rcpp::Rcout<<"Spectral pre-processing...\n";
  
  //Run preprocessing in mutli-threading
  runMSIProcessingCpp();

  //Retrun first iteration lags
  NumericVector LagsLow(numPixels);
  NumericVector LagsHigh(numPixels);
  for( int i = 0; i < numPixels; i++)
  {
    LagsLow[i] = mLags[i].lagLow;
    LagsHigh[i] = mLags[i].lagHigh;
  }
  
  //Get the imzMLwriters offset and average/base spectrums
  Rcpp::List offsetLst;
  Rcpp::List averageSpectraLst;
  Rcpp::List baseSpectraLst;
  for( unsigned int i = 0; i < ioObj->get_images_count(); i++)
  {
    offsetLst.push_back(ioObj->get_OffsetsLengths(i));
    averageSpectraLst.push_back(ioObj->get_AverageSpectrum(i));
    baseSpectraLst.push_back(ioObj->get_BaseSpectrum(i));
  }
  
  return List::create( Named("LagLow") = LagsLow, 
                       Named("LagHigh") = LagsHigh, 
                       Named("Offsets") = offsetLst, 
                       Named("AverageSpectra") = averageSpectraLst, 
                       Named("BaseSpectra") = baseSpectraLst
                        );
}

void MTPreProcessing::ProcessingFunction(int threadSlot)
{
  //Process each spectrum in the current loaded cube
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
   
   //TODO add the baseline reduction
   
    if(bEnableSmoothing)
    {
      smoothObj[threadSlot]->smoothSavitzkyGolay( massAxis.begin(),
                                                  cubes[threadSlot]->dataInterpolated[j],
                                                  massAxis.length(), 
                                                  cubes[threadSlot]->dataOriginal[j].imzMLmass.data(),
                                                  cubes[threadSlot]->dataOriginal[j].imzMLintensity.data(),
                                                  cubes[threadSlot]->dataOriginal[j].imzMLmass.size()
                                                  ); 
      
    }
    
    if(bEnableAlignment)
    {
      mLags[cubes[threadSlot]->dataOriginal[j].pixelID] = alngObj[threadSlot]->AlignSpectrum( cubes[threadSlot]->dataInterpolated[j], 
                                                                                              cubes[threadSlot]->dataOriginal[j].imzMLmass.data(),
                                                                                              cubes[threadSlot]->dataOriginal[j].imzMLintensity.data(),
                                                                                              cubes[threadSlot]->dataOriginal[j].imzMLmass.size()
                                                                                              );
    }
    
   
    if(bEnableSmoothing || bEnableAlignment)
    {
      if(cubes[threadSlot]->dataOriginal[j].imzMLmass.size() == 0)
      {
        //Continuous mode
        BitDepthReduction(cubes[threadSlot]->dataInterpolated[j], cubes[threadSlot]->ncols, threadSlot);
      }
      else
      {
        //Processed mode
        if(cubes[threadSlot]->dataOriginal[j].imzMLintensity.size() <= cubes[threadSlot]->ncols)
        {
          BitDepthReduction(cubes[threadSlot]->dataOriginal[j].imzMLintensity.data(), cubes[threadSlot]->dataOriginal[j].imzMLintensity.size(), threadSlot);
        }
      }
    }
   
  }
}

#define MIN_MANTISSA_BITS 4
#define NOISE_THRESHOLD_LOWER 0.5
#define NOISE_THRESHOLD_UPPER 10.0

void MTPreProcessing::BitDepthReduction(double *data, int dataLength, int noiseModelThreadSlot)
{
  int resolution_bits;
  double *noise_floor = new double[dataLength];
  memcpy(noise_floor, data, sizeof(double)*dataLength);
  noiseModel[noiseModelThreadSlot]->NoiseEstimationFFTExpWin(noise_floor, dataLength, NoiseWinSize);
  
  double m, n;
  unsigned long long *ptr;
  for( int i = 0; i < dataLength; i++)
  {
    //Calculate required bit-depth according to the distance to the noise floor
    m = (52 - MIN_MANTISSA_BITS)/( noise_floor[i] * ( NOISE_THRESHOLD_UPPER - NOISE_THRESHOLD_LOWER ) );
    n = MIN_MANTISSA_BITS - m*NOISE_THRESHOLD_LOWER*noise_floor[i];
    resolution_bits = (int)round( m*data[i] + n );
    
    //Limit resolution bits to a double mantissa valid range
    resolution_bits = resolution_bits > 52 ? 52 : resolution_bits;
    resolution_bits = resolution_bits <  1 ?  1 : resolution_bits;
    
    //Apply the mask to double's mantissa
    ptr = (unsigned long long*) (data + i);
    *ptr &= maskLUT_double[resolution_bits-1];  
    
  }
  
  delete[] noise_floor;
}

// [[Rcpp::export]]
List CRunPreProcessing( Rcpp::List rMSIObj_list,int numOfThreads, double memoryPerThreadMB, 
                     Rcpp::Reference preProcessingParams, Rcpp::NumericVector reference, 
                     Rcpp::StringVector uuid, Rcpp::String outputDataPath, Rcpp::StringVector imzMLoutFnames,
                     Rcpp::NumericVector commonMassAxis)
{
   MTPreProcessing myPreProcessing(rMSIObj_list, numOfThreads, memoryPerThreadMB,
                                  preProcessingParams, reference, 
                                  uuid, outputDataPath, imzMLoutFnames, 
                                  commonMassAxis);
   return myPreProcessing.Run();
}
