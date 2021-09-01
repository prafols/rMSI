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
#include "mtnormalization.h"
using namespace Rcpp;

MTNormalization::MTNormalization(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, Rcpp::NumericVector commonMassAxis) : 
  ThreadingMsiProc(rMSIObj_list, numberOfThreads, memoryPerThreadMB, commonMassAxis)
{
  for( int i = 0; i < rMSIObj_list.length(); i++)
  {
    unsigned int imagepixels = (as<NumericMatrix>(as<List>(rMSIObj_list[i])["pos"])).nrow();
    NumericVector vTIC(imagepixels);
    NumericVector vRMS(imagepixels);
    DataFrame dfnorms = DataFrame::create(Named("TIC") = vTIC, Named("RMS") = vRMS);
    lstNorms.push_back(dfnorms);
  }
}

MTNormalization::~MTNormalization()
{
  //Empty destructor
}

List MTNormalization::Run()
{
  Rcpp::Rcout<<"Calculating normalizations...\n";
  
  //Run in multi-threading
  runMSIProcessingCpp();
  
  return lstNorms;
}


void MTNormalization::ProcessingFunction(int threadSlot)
{
  //Perform the average value of each mass channel in the current loaded cube
  double TIC, RMS;
  for (int j = 0; j < cubes[threadSlot]->nrows; j++)
  {
    TIC = 0.0;
    RMS = 0.0;
    
    for (int k= 0; k < cubes[threadSlot]->ncols; k++)
    {
      TIC += cubes[threadSlot]->dataInterpolated[j][k];
      RMS += (cubes[threadSlot]->dataInterpolated[j][k] * cubes[threadSlot]->dataInterpolated[j][k]);
    }
    RMS = sqrt(RMS);
    
    int imgID = ioObj->getImageIndex(cubes[threadSlot]->cubeID, j);
    int pixelID = ioObj->getPixelId(cubes[threadSlot]->cubeID, j);
    
    mutex_copyData.lock();
    as<NumericVector>((as<DataFrame>(lstNorms[imgID]))["TIC"])[pixelID] = TIC;
    as<NumericVector>((as<DataFrame>(lstNorms[imgID]))["RMS"])[pixelID] = RMS;
    mutex_copyData.unlock();
  }
}

// Calculate the average spectrum from a list of rMSI objects.
// [[Rcpp::export]]
List CNormalizations(Rcpp::List rMSIObj_list, 
                               int numOfThreads, 
                               double memoryPerThreadMB,
                               Rcpp::NumericVector commonMassAxis)
{
  List out;
  
  try
    {
    MTNormalization myNorms(rMSIObj_list, 
                      numOfThreads, 
                      memoryPerThreadMB,
                      commonMassAxis);
    out = myNorms.Run();
    }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}
