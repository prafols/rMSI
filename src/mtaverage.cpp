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
#include "mtaverage.h"
using namespace Rcpp;

MTAverage::MTAverage(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, double minTIC, double maxTIC) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB),
  TICmin(minTIC), 
  TICmax(maxTIC)
{
  sm = new double*[ioObj->getNumberOfCubes()];
  validPixelCount = new unsigned int[ioObj->getNumberOfCubes()];  
  
  for(int i = 0; i < ioObj->getNumberOfCubes() ; i++)
  {
    sm[i] = new double[ioObj->getMassAxisLength()];
    validPixelCount[i] = 0;
  }
  
  for (int i = 0; i < ioObj->getNumberOfCubes() ; i++)
  {
    for (int j = 0; j < ioObj->getMassAxisLength() ; j++)
    {
    sm[i][j] = 0;
    }
  }
  
  //Get normalizations
  for(int i = 0; i < rMSIObj_list.size(); i++)
  {
    lNormalizations.push_back(Rcpp::as<Rcpp::DataFrame>((Rcpp::as<Rcpp::List>(rMSIObj_list[i]))["normalizations"]));
  }
  
}

MTAverage::~MTAverage()
{
  for (int i = 0; i < ioObj->getNumberOfCubes(); i++)
  {
     delete[] sm[i];
  }
  
  delete[] sm;
  delete[] validPixelCount;
}

NumericVector MTAverage::Run()
{
  Rcpp::Rcout<<"Calculating overall average spectrum...\n";
  
  NumericVector avg(ioObj->getMassAxisLength());
  for (int j = 0; j < ioObj->getMassAxisLength(); j++)
  {
    avg[j] = 0;
  }
  //Run in multi-threading
  runMSIProcessingCpp();
  
  //Merging all the cube's partial average spectrum and getting the total average spectrum
  unsigned int CubeCount = 0;
  for (int i = 0; i < ioObj->getNumberOfCubes(); i++)
  {
    if(validPixelCount[i] > 0)
    {
      for (int j = 0; j < ioObj->getMassAxisLength(); j++)
      {
        avg[j] += sm[i][j];
      }
      CubeCount++;
    }
  }
  
  if(CubeCount > 0)
  {
    for (int j = 0; j < ioObj->getMassAxisLength(); j++)
    {
      avg[j] /= (double)CubeCount;
    }
  }
  
  return avg;
}


void MTAverage::ProcessingFunction(int threadSlot)
{
  //Perform the average value of each mass channel in the current loaded cube
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    int imageIndex = ioObj->getImageIndex(cubes[threadSlot]->cubeID, j);
    int pixelIndex = ioObj->getPixelId(cubes[threadSlot]->cubeID, j);
    double TICval = (Rcpp::as<Rcpp::NumericVector>((Rcpp::as<Rcpp::DataFrame>(lNormalizations[imageIndex])["TIC"])))[pixelIndex];
    
    if(TICval >= TICmin && TICval <= TICmax)
    {
      for (int k= 0; k < cubes[threadSlot]->ncols; k++)
      {
        sm[cubes[threadSlot]->cubeID][k] += (cubes[threadSlot]->dataInterpolated[j][k])/TICval; //Average with TIC Normalization
        validPixelCount[cubes[threadSlot]->cubeID]++;
      }
    }
  }
  
  //Partial average
  if(validPixelCount[cubes[threadSlot]->cubeID] > 0)
  {
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      sm[cubes[threadSlot]->cubeID][k] /= (double)(validPixelCount[cubes[threadSlot]->cubeID]);
    }
  }
}

// Calculate the average spectrum from a list of rMSI objects.
// [[Rcpp::export]]
NumericVector COverallAverageSpectrum(Rcpp::List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB,
                               double minTIC, double maxTic)
{
  MTAverage myAverage(rMSIObj_list, 
                      numOfThreads, 
                      memoryPerThreadMB, 
                      minTIC,
                      maxTic);
 
  return myAverage.Run();
}
