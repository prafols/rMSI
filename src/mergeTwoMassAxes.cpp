/*************************************************************************
 *     rMSI - R package for MSI data processing
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
#include <limits>

using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib rMSI, .registration = TRUE

//' MergeMassAxis.
//' 
//' Merges two mass axis in a single one using an apropiate bin size.
//' The resulting mass axis will display a bin size equal to the minimum of two supplied vectors. 
//' The bin size is calculated relative to the m/z for better accuracy.
//' The resulting mass axis range is calculated using the common range between the two mass axis.
//' If there is no overlao between the two mass axis range an error will be raised.
//' 
//' @param mz1 the first mass axis to merge.
//' @param mz2 the second mass axis to merge.
//' 
//' @return a list containing the common mass axis that represents mz1 and mz1 accurately and a boolean indicating if and error was raised.
//' @export
//' 
// [[Rcpp::export]]
List MergeMassAxis(NumericVector mz1, NumericVector mz2)
{
  int i1 = 0; //Iterator over mz1
  int i2 = 0; //Iterator over mz2
  int iN = 0; //Iterator over mzNew
  double dist, binSize, dist1 = 0.0, dist2 = 0.0;
  double *mzNew = new double[mz1.length() + mz2.length()]; //Allocate memory for the worst case
  
  while( i1 < mz1.length() || i2 < mz2.length())
  {
    if(i2 >= mz2.length())
    {
      //We run out of mz2
      mzNew[iN] = mz1[i1];
      i1++;
    }
    else if(i1 >= mz1.length())
    {
      //We run out of mz1
      mzNew[iN] = mz2[i2];
      i2++;
    }
    else
    {
      //Calculate the binsize in current region
      binSize = std::numeric_limits<double>::max(); //Start with something really high to allow min to work.
      if( i1 > 0 )
      {
        binSize = (mz1[i1] - mz1[i1-1]) < binSize ? (mz1[i1] - mz1[i1-1]) : binSize;
      }
      
      if( i1+1 < mz1.length() )
      {
        binSize = (mz1[i1+1] - mz1[i1]) < binSize ? (mz1[i1+1] - mz1[i1]) : binSize;
      }
      
      if( i2 > 0 )
      {
        binSize = (mz2[i2] - mz2[i2-1]) < binSize ? (mz2[i2] - mz2[i2-1]) : binSize;
      }
      
      if( i2+1 < mz2.length() )
      {
        binSize = (mz2[i2+1] - mz2[i2]) < binSize ? (mz2[i2+1] - mz2[i2]) : binSize;
      }
      
      dist = fabs(mz1[i1] -  mz2[i2]);
      if( iN > 0 )
      {
        dist1 = mz1[i1] - mzNew[iN-1];
        dist2 = mz2[i2] - mzNew[iN-1];
      }
      
      //Resampling
      if( mz1[i1] < mz2[i1] )
      {
        if( dist1 < binSize && iN > 0)
        {
          mzNew[iN] = mzNew[iN-1] + binSize; 
        }
        else
        {
          mzNew[iN] = mz1[i1];
        }
      }
      else
      {
        if( dist2 < binSize && iN > 0)
        {
          mzNew[iN] = mzNew[iN-1] + binSize; 
        }
        else
        {
          mzNew[iN] = mz2[i1];
        }
      }
        
      if(dist < binSize)
      {
        //Merging elements
        i1++;
        i2++;
      }
      else if(mz1[i1] < mz2[i2])
      {
        //Inserting mz1
        i1++;
      }
      else
      {
        //Inserting mz2
        i2++;
      }
    }
    iN++;
  }
  
  //if the new mass axis is empty then mz1 and mz2 are non-overlaping vectors, so an error is raised.
  bool bError = false;
  if( iN == 0)
  {
    bError = true;
    iN = 1;
    mzNew[0] = -1; //Mark as error using -1 (not possible to have an m/z value of -1)
  }
  
  //Copy to a NumericVector
  NumericVector resultMass(iN);
  memcpy(resultMass.begin(), mzNew, sizeof(double)*iN);
  delete[] mzNew;
  return List::create( Named("mass") = resultMass, Named("error") = bError );
}
