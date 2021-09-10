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

#ifndef MT_REPLACEZEROPEAKS_H
#define MT_REPLACEZEROPEAKS_H
#include <Rcpp.h>
#include "peakpicking.h"
#include "threadingmsiproc.h"

class MTFillPeaks : public ThreadingMsiProc 
{
  public:
    
    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    // preProcessingParams: An R reference class with the pre-processing parameters.
    // mass: a numeric vector with the common mass axis
    // commonMassAxis: The common mass axis used to process and interpolate multiple datasets.
    // peakMatrix: a peak matrix as it is returned by the peak binning algorithm
    MTFillPeaks(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB, 
                Rcpp::Reference preProcessingParams,
                Rcpp::NumericVector commonMassAxis,
                Rcpp::List peakMatrix);
    
    ~MTFillPeaks();
    
    //Exectur a full imatge processing using threaded methods, nothing is returned becasue the peakmatrix is directely modified
    void Run(); 
    
  private:
    unsigned int *replacedZerosCounters;
    unsigned int *mass_index;
    
    //In theory, Rcpp passes arguments like Rcpp::Lists an matrices as reference... so it should be as simple as follows
    Rcpp::NumericVector pkMatmass;
    Rcpp::NumericMatrix pkMatintensity;
    Rcpp::NumericMatrix pkMatarea;
    Rcpp::NumericMatrix pkMatsnr;
    PeakPicking **peakObj;
  
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
};
#endif
