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

using namespace Rcpp;

// Function to reduce the number of data points in part of spectrum to plot it faster.
// This is implemented in C to provid an efficient implementation.
// [[Rcpp::export]]
List ReduceDataPointsC(NumericVector mass, NumericVector intensity, double massMin, double massMax, int npoints)
{
  //Locate mass range index min
  int imin;
  for( imin = 0; imin < mass.length(); imin++)
  {
    if( mass[imin] >= massMin )
    {
      break;
    }
  }
  
  //Locate mass range index max
  int imax;
  for( imax = mass.length() - 1; imax >= 0; imax--)
  {
    if( mass[imax] <= massMax )
    {
      break;
    }
  }
  
  //Adjust extrems to handles properly N < npoints
  imin = imin == 0 ? imin : (imin - 1);
  imax = imax == ( mass.length()  - 1) ? ( mass.length()  - 1) : (imax + 1);

  int N = npoints > (imax - imin + 1) ? (imax - imin + 1) : npoints; 
  NumericVector mzData(N);
  NumericVector inData(N);
  
  //Indexing converting
  const double m = ((double)(N - 1))/((double)(imax - imin));
  const double n = - m*((double)imin);
  
  //Process the selected mass range
  int binID = 0;
  int binID_ant = 0;
  for(int i = imin; i <= imax; i++)
  {
    //Data reduction
    binID = (int)(m * (double)i + n);
    binID = binID < 0 ? 0 : binID;
    binID = binID > (N - 1) ? (N - 1) : binID;
      
    if( intensity[i] >= inData[binID]  || binID != binID_ant)
    {
      inData[binID] = intensity[i];
      mzData[binID] = mass[i];
      binID_ant = binID;
    }
  }
  
  return List::create( Named("mass") = mzData, Named("intensity") = inData);
}
