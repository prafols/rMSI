/*************************************************************************
 *     rMSIproc - R package for MSI data processing
 *     Copyright (C) 2017 Pere Rafols Soler
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
#include <memory>
#include "mtfillpeaks.h"
using namespace Rcpp;

MTFillPeaks::MTFillPeaks(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, 
                         Rcpp::Reference preProcessingParams,
                         Rcpp::NumericVector commonMassAxis,
                         Rcpp::List peakMatrix):
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis)
{
  replacedZerosCounters = new unsigned int[numOfThreadsDouble];
  for( int i = 0; i < numOfThreadsDouble; i++)
  {
    replacedZerosCounters[i] = 0;
  }
  
  //Get the peak-picking params
  Rcpp::Reference peakPickingParams = preProcessingParams.field("peakpicking");
  int peakWinSize = peakPickingParams.field("WinSize");
  int peakInterpolationUpSampling = peakPickingParams.field("overSampling");
  
  peakObj = new PeakPicking*[numOfThreadsDouble];
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    peakObj[i] = new PeakPicking(peakWinSize, massAxis.begin(), massAxis.length(), peakInterpolationUpSampling );  
  }
  
  //Get the peak matrices parts, in theory this is passed by reference accoring to Rcpp documentation
  pkMatintensity = as<NumericMatrix>(peakMatrix["intensity"]);
  pkMatarea = as<NumericMatrix>(peakMatrix["area"]);
  pkMatsnr = as<NumericMatrix>(peakMatrix["SNR"]);
  pkMatmass = as<NumericVector>(peakMatrix["mass"]);
  
  //Build the mass index vector
  Rcout<<"Creating the mass index vector...\n";
  mass_index = new unsigned int[pkMatmass.length()];
  int lastMassId = 1;
  double mass_distance_ant, mass_distance;
  for(int i = 0; i < pkMatmass.length(); i++)
  {
    mass_distance_ant = fabs( pkMatmass[i] -  massAxis[lastMassId-1]);
    for( int imass = lastMassId; imass <  massAxis.length(); imass++)
    {
      mass_distance = fabs( pkMatmass[i] - massAxis[imass]);
      if( mass_distance > mass_distance_ant )
      {
        mass_index[i] = imass - 1;
        //Rcpp::Rcout<<"DBG mass_index[ " << i << " ] = " << mass_index[i] << "\n"; 
        lastMassId = imass;
        break;
      }
      mass_distance_ant = mass_distance;
    }
  }
}

MTFillPeaks::~MTFillPeaks()
{
  for(int i = 0; i < numOfThreadsDouble; i++)
  {
    delete peakObj[i];
  }
  delete[] peakObj;
  delete[] replacedZerosCounters;
  delete[] mass_index;
}

void MTFillPeaks::Run()
{
  Rcout<<"Filling peaks below SNR threshold...\n";
  runMSIProcessingCpp();
  int replacedZerosSum = 0;
  for( int i = 0; i < numOfThreadsDouble; i++)
  {
    replacedZerosSum += replacedZerosCounters[i];
  }
  Rcout << "A total of " << replacedZerosSum << " low intensity peaks were retrieved\n";
}

void MTFillPeaks::ProcessingFunction(int threadSlot)
{
  unsigned int peakMat_row_index;
  for( int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    //Get the row in the peak matrix of the current pixel
    peakMat_row_index = ioObj->getPeakMatrixRow(cubes[threadSlot]->cubeID, j);
    
    for( int imass = 0; imass < pkMatintensity.ncol(); imass++)
    {
      //Locate the zeros.....
      if( pkMatintensity(peakMat_row_index, imass) == 0.0 )
      {
        //Count the replaced zero
        replacedZerosCounters[threadSlot]++;
        
        //Fill matrix position with proper intensity
        pkMatintensity(peakMat_row_index, imass) = cubes[threadSlot]->dataInterpolated[j][mass_index[imass]]; 
        pkMatarea(peakMat_row_index, imass) = peakObj[threadSlot]->predictPeakArea(cubes[threadSlot]->dataInterpolated[j], mass_index[imass]);
      }
    }
  }
}

//Returning void since in theory I  can directely modify the peak matrix without returning anything
// [[Rcpp::export]]
void CRunFillPeaks( Rcpp::List rMSIObj_list,int numOfThreads, double memoryPerThreadMB, 
                    Rcpp::Reference preProcessingParams, 
                    Rcpp::NumericVector commonMassAxis,
                    Rcpp::List peakMatrix)
{
  try
  {
    MTFillPeaks myFillPeaks(rMSIObj_list, numOfThreads, memoryPerThreadMB,
                            preProcessingParams,
                            commonMassAxis,
                            peakMatrix);
    
    myFillPeaks.Run();
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
}
