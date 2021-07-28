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

#ifndef MT_NORMALIZATION_H
  #define MT_NORMALIZATION_H
#include <Rcpp.h>
#include "threadingmsiproc.h"

class MTNormalization : public ThreadingMsiProc 
{
  public:

    //Constructor arguments:
    // rMSIObj_list: A list of rMSI objects to process
    // numberOfThreads: Total number of threads to use during processing
    // memoryPerThreadMB: Maximum memory allocated by each thread in MB. The total allocated memory will be: 2*numberOfThreads*memoryPerThreadMB
    MTNormalization(Rcpp::List rMSIObj_list, int numberOfThreads, double memoryPerThreadMB);
    ~MTNormalization();
    
    //Execute a full imatge processing using threaded methods
    Rcpp::List Run();
    
  private:
    Rcpp::List lstNorms; //List with multiple normalization data frames
    
    //Thread Processing function definition
    void ProcessingFunction(int threadSlot);
};
#endif
