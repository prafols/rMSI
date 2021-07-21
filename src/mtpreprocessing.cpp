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
                                 Rcpp::StringVector uuid,  Rcpp::String outputImzMLPath, Rcpp::StringVector outputImzMLfnames) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, true, uuid, outputImzMLPath, outputImzMLfnames)
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
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    smoothObj[i] = new Smoothing(smoothinKernelSize);
    
    alngObj[i] = new LabelFreeAlign(massAxis.begin(), reference.begin(), massAxis.length(), bilinear, 
                                    alignIterations, lagRefLow, lagRefMid, lagRefHigh,
                                    maxShiftppm,  fftOverSampling, winSizeRelative);
  }

  mLags  =  new LabelFreeAlign::TLags[numPixels]; 
}

MTPreProcessing::~MTPreProcessing()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete smoothObj[i];
    delete alngObj[i];
  }
  delete[] smoothObj;
  delete[] alngObj;
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
    
   //TODO add the average accumulation! don't do it here its a mess!

   
    //TODO add the bitdepth reduction!
  }
}

// [[Rcpp::export]]
List CRunPreProcessing( Rcpp::List rMSIObj_list,int numOfThreads, double memoryPerThreadMB, 
                     Rcpp::Reference preProcessingParams, Rcpp::NumericVector reference, 
                     Rcpp::StringVector uuid, Rcpp::String outputDataPath, Rcpp::StringVector imzMLoutFnames)
{
   MTPreProcessing myPreProcessing(rMSIObj_list, numOfThreads, memoryPerThreadMB,
                                  preProcessingParams, reference, 
                                  uuid, outputDataPath, imzMLoutFnames);
   return myPreProcessing.Run();
}
