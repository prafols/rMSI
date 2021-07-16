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

MTAverage::MTAverage(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB)
{
  Rcpp::Rcout<<"DBG trap1\n";
  sm = new double*[ioObj->getNumberOfCubes()];
  Rcpp::Rcout<<"DBG trap2\n";
  for(int i = 0; i < ioObj->getNumberOfCubes() ; i++)
  {
    sm[i] = new double[ioObj->getMassAxisLength()];
  }
  Rcpp::Rcout<<"DBG trap3\n";
  for (int i = 0; i < ioObj->getNumberOfCubes() ; i++)
  {
    for (int j = 0; j < ioObj->getMassAxisLength() ; j++)
    {
    sm[i][j] = 0;
    }
  }
  Rcpp::Rcout<<"DBG trap4\n";
}

MTAverage::~MTAverage()
{
  for (int i = 0; i < ioObj->getNumberOfCubes(); i++)
  {
     delete[] sm[i];
  }
  
  delete[] sm;
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
  for (int i = 0; i < ioObj->getNumberOfCubes(); i++)
  {
    for (int j = 0; j < ioObj->getMassAxisLength(); j++)
    {
      avg[j] += sm[i][j];
    }
  }
  
  for (int j = 0; j < ioObj->getMassAxisLength(); j++)
  {
    avg[j] /= numPixels;    //Dividir per el nombre total de pixels.
  }
  
  return avg;
}


void MTAverage::ProcessingFunction(int threadSlot)
{
  //Perform the average value of each mass channel in the current loaded cube
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
       sm[cubes[threadSlot]->cubeID][k] += cubes[threadSlot]->dataInterpolated[j][k];
    }
  }
}

// Calculate the average spectrum from a list of rMSI objects.
// [[Rcpp::export]]
NumericVector COverallAverageSpectrum(Rcpp::List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB)
{
  
  Rcpp::Rcout<<"DBG trap10\n";
  MTAverage myAverage(rMSIObj_list, 
                      numOfThreads, 
                      memoryPerThreadMB);
  Rcpp::Rcout<<"DBG trap11\n";
  
  return myAverage.Run();
}
