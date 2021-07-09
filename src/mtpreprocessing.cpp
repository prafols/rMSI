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
                                 Rcpp::Reference preProcessingParams, Rcpp::NumericVector mass, Rcpp::NumericVector reference) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB)
{
  bool bContinuousMode = false; //TODO this is very imporant! get it from the data and set it properly here!!!! also check all data is in the same format, otherwise abort!
  
  //TODO add smoothing and baseline params here!
  
  //Get the alignment params
  Rcpp::Reference alignmentParams = preProcessingParams.field("alignment");
  bool bilinear = alignmentParams.field("bilinear");
  int alignIterations = alignmentParams.field("iterations");
  double maxShiftppm = alignmentParams.field("maxShiftppm");
  double lagRefLow = alignmentParams.field("refLow");
  double lagRefMid = alignmentParams.field("refMid");
  double lagRefHigh = alignmentParams.field("refHigh");
  int fftOverSampling = alignmentParams.field("overSampling");
  double winSizeRelative = alignmentParams.field("winSizeRelative");
  
  alngObj = new LabelFreeAlign*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    alngObj[i] = new LabelFreeAlign(mass.begin(), reference.begin(), mass.length(), bilinear, 
                                    &fftSharedMutex, alignIterations, 
                                    lagRefLow, lagRefMid, lagRefHigh,
                                    maxShiftppm,  fftOverSampling, winSizeRelative);
  }

  mLags  =  new LabelFreeAlign::TLags[numPixels]; 
}

MTPreProcessing::~MTPreProcessing()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete alngObj[i];
  }
  delete[] alngObj;
  delete[] mLags;
}

List MTPreProcessing::Run()
{
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
  
  return List::create( Named("LagLow") = LagsLow, Named("LagHigh") = LagsHigh);
}

void MTPreProcessing::ProcessingFunction(int threadSlot)
{
  //Perform alignment of each spectrum in the current loaded cube
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    mLags[cubes[threadSlot]->dataOriginal[j].pixelID] = alngObj[threadSlot]->AlignSpectrum( cubes[threadSlot]->dataInterpolated[j], 
                                                                                            cubes[threadSlot]->dataOriginal[j].imzMLmass.data(),
                                                                                            cubes[threadSlot]->dataOriginal[j].imzMLmass.size());
  }
}

// [[Rcpp::export]]
List FullImageAlign( Rcpp::List rMSIObj_list,int numOfThreads, double memoryPerThreadMB, 
                     Rcpp::Reference preProcessingParams, Rcpp::NumericVector mass, Rcpp::NumericVector reference)
{
  MTPreProcessing myPreProcessing(rMSIObj_list, numOfThreads, memoryPerThreadMB,
                                  preProcessingParams, mass, reference);

  return myPreProcessing.Run();
}






