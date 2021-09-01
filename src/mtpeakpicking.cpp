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
#include "mtpeakpicking.h"
using namespace Rcpp;

MTPeakPicking::MTPeakPicking(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB,
                             Rcpp::Reference preProcessingParams,
                             Rcpp::StringVector uuid, Rcpp::String outputImzMLPath, Rcpp::StringVector outputImzMLfnames, 
                             Rcpp::NumericVector commonMassAxis) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis, DataCubeIOMode::PEAKLIST_STORE, uuid, outputImzMLPath, outputImzMLfnames)
{
  //Get the peak-picking params
  Rcpp::Reference peakPickingParams = preProcessingParams.field("peakpicking");
  minSNR = peakPickingParams.field("SNR");
  int peakWinSize = peakPickingParams.field("WinSize");
  int peakInterpolationUpSampling = peakPickingParams.field("overSampling");
  
  peakObj = new PeakPicking*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    peakObj[i] = new PeakPicking(peakWinSize, massAxis.begin(), massAxis.length(), peakInterpolationUpSampling );  
  }
}

MTPeakPicking::~MTPeakPicking()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete peakObj[i];
  }
  delete[] peakObj;
}

Rcpp::List MTPeakPicking::Run()
{
  Rcpp::Rcout<<"\n\nPicking peaks...\n";
  
  //Run peak-picking in mutli-threading
  runMSIProcessingCpp();
  
  //Get the imzMLwriters offset and average/base spectrums
  Rcpp::List offsetLst;
  for( unsigned int i = 0; i < ioObj->get_images_count(); i++)
  {
    offsetLst.push_back(ioObj->get_OffsetsLengths(i));
  }
  
  return offsetLst;
}

void MTPeakPicking::ProcessingFunction(int threadSlot)
{
  //Perform peak-picking of each spectrum in the current loaded cube
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
     cubes[threadSlot]->peakLists[j] = peakObj[threadSlot]->peakPicking( cubes[threadSlot]->dataInterpolated[j], minSNR ); 
  }
}

// [[Rcpp::export]]
Rcpp::List CRunPeakPicking(   Rcpp::List rMSIObj_list,int numOfThreads, double memoryPerThreadMB, 
                        Rcpp::Reference preProcessingParams, 
                        Rcpp::StringVector uuid, Rcpp::String outputDataPath, Rcpp::StringVector imzMLoutFnames,
                        Rcpp::NumericVector commonMassAxis)
{
  Rcpp::List out;
  try
  {
    MTPeakPicking myPeakPicking (rMSIObj_list, numOfThreads, memoryPerThreadMB,
                               preProcessingParams,
                               uuid, outputDataPath, imzMLoutFnames, 
                               commonMassAxis);
  
    out = myPeakPicking.Run();
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  
  return out;
}
