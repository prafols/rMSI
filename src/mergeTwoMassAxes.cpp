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

//Simple peak detection with bin size calculation
List simplePeakDetect(NumericVector mz, NumericVector intensity)
{
  double slope_ant = 0, slope;
  NumericVector pkMass;
  NumericVector pkIntensity;
  NumericVector pkIndex;
  NumericVector pkBinSize;
  NumericVector pkIndexLeft;
  NumericVector pkIndexRight;
  for( int i = 0; i < (mz.length()-1); i++ )
  {
    if( i < (mz.length()-1) )
    {
      slope = intensity[i+1] - intensity[i];
    }
    if( (slope_ant >  0 ) && (slope < 0 ) && i > 0 ) 
    {
      pkMass.insert(pkMass.end(), mz[i] );
      pkIntensity.insert(pkIntensity.end(), intensity[i] );
      pkIndex.insert(pkIndex.end(), i ); //inserting C index
    }
    slope_ant = slope;
  }
  
  //Calculate the mass bin size for each peak
  for( int i = 0; i < pkIndex.length(); i++)
  {
    double localBinSize = 0;
    int count = 0;
    //Left region
    if(pkIndex[i] > 0)
    {
      for( int j = pkIndex[i]; j >= 0; j--)
      {
        if( (intensity[j] <= intensity[j-1])  )
        {
          pkIndexLeft.insert(pkIndexLeft.end(), j ); //inserting C index
          break; //End of peak
        }
        localBinSize += (mz[j] - mz[j-1]);
        count++;
      }
    }
    
    //Right region
    if(pkIndex[i] < (mz.length()-1))
    {
      for( int j = pkIndex[i]; j < (mz.length()-1); j++)
      {
        if( (intensity[j] <= intensity[j+1])  )
        {
          pkIndexRight.insert(pkIndexRight.end(), j ); //inserting C index
          break; //End of peak
        }
        localBinSize += (mz[j+1] - mz[j]);
        count++;
      }
    }
    if(count >0 )
    {
      localBinSize /= (double)count;
      pkBinSize.insert(pkBinSize.end(), localBinSize);
    }
  }
  
  return List::create( Named("mass") = pkMass,
                       Named("intensity") = pkIntensity, 
                       Named("index") = pkIndex,
                       Named("indexLeft") = pkIndexLeft,
                       Named("indexRight") = pkIndexRight,
                       Named("BinSize") = pkBinSize );
}


//' CalcMassAxisBinSize.
//' 
//' Calc the bin size of a mass axis at each mass channels using simple peak-picking information.
//' 
//' @param mass the mass axis.
//' @param intensity the intensity of a given spectrum.
//' 
//' @return the bin size of each m/z channel.
//' @export
//' 
// [[Rcpp::export]]
NumericVector CalcMassAxisBinSize(NumericVector mass, NumericVector intensity)
{
  //Calc bin size at each input vector position
  List pks2 = simplePeakDetect(mass, intensity);
  NumericVector bins2(mass.length());
  int currInd = 0;
  NumericVector indeces = pks2["indexRight"];
  NumericVector binss = pks2["BinSize"];
  for(int i=0; i < indeces.length(); i++ )
  {
    for( int  j = currInd; j <= indeces[i]; j++)
    {
      bins2[j] = binss[i];
    }
    currInd = indeces[i]; 
  }
  while(currInd < bins2.length() )
  {
    bins2[currInd] = binss[binss.length()-1];
    currInd++;
  }
  
  return bins2;
}

//' MergeMassAxis.
//' 
//' Merges two mass axis in a single one using an apropiate bin size.
//' The resulting mass axis will display a bin size equal to the minimum of two supplied vectors. 
//' The bin size must be supplied along each input mass axis.
//' The first mass axis (mz1) can be a zero-length vector.
//' 
//' @param mz1 the first mass axis to merge.
//' @param bins1 the bins size for the first mass axis.
//' @param mz2 the second mass axis to merge.
//' @param intensity2 the spectral intensities corresponding to the second mass axis.
//' 
//' @return a list containing the common mass axis that represents mz1 and mz1 accurately and a boolean indicating if and error was raised.
//' @export
//' 
// [[Rcpp::export]]
List MergeMassAxis(NumericVector mz1, NumericVector bins1, NumericVector mz2, NumericVector bins2)
{
  //Error check
  if(mz1.length() != bins1.length())
  {
    stop("Error: mz1 and bins1 must be the same length.\n");
  }
  if(mz2.length() != bins2.length())
  {
    stop("Error: mz2 and intensity2 must be the same length.\n");
  }
  if(mz2.length() == 0)
  {
    stop("Error: mz2 does not contain any element.\n");
  }
  
  //Vectors concatenation and sorting
  double *newMz = new double[mz1.length() + mz2.length()];
  double *newBins = new double[mz1.length() + mz2.length()];
  int i1 = 0;
  int i2 = 0;
  int inew = 0;
  if( mz1.length() == 0 )
  {
    //mz1 is empty, so just copy mz2
    for(int i = 0; i < mz2.length(); i++)
    {
      newMz[inew] = mz2[i];
      newBins[inew] = bins2[i];
      inew++;
    }
  }
  else
  {
    //Merge mz1 and mz2
    while(true)
    {
      if(i1 >= mz1.length() && i2 >= mz2.length())
      {
        //We run out of elements on both mass axes
        break;
      }
      else if(i1 >= mz1.length())
      {
        //We run out of elements in mz1
        if(inew == 0)
        {
          newMz[inew] = mz2[i2];
          inew++;
        }
        else if( (mz2[i2] -  newMz[inew - 1]) >= bins2[i2] )
        {
          newMz[inew] = mz2[i2];
          newBins[inew] = bins2[i2];
          inew++;
        }
        i2++;
      }
      else if(i2 >= mz2.length())
      {
        //We run out of elements in mz2
        if(inew == 0)
        {
          newMz[inew] = mz1[i1];
          inew++;
        }
        else if( (mz1[i1] -  newMz[inew - 1]) >= bins1[i1] )
        {
          newMz[inew] = mz1[i1];
          newBins[inew] = bins1[i1];
          inew++;
        }
        i1++;
      }
      else
      {
        //There are remaining elements in both mass axes
        if(mz1[i1] < mz2[i2])
        {
          if(inew == 0)
          {
            newMz[inew] = mz1[i1];
            inew++;
          }
          else if( (mz1[i1] -  newMz[inew - 1]) >= bins1[i1] )
          {
            newMz[inew] = mz1[i1];
            newBins[inew] = bins1[i1];
            inew++;
          }
          i1++;
        }
        else
        {
          if(inew == 0)
          {
            newMz[inew] = mz2[i2];
            inew++;
          }
          else if( (mz2[i2] -  newMz[inew - 1]) >= bins2[i2] )
          {
            newMz[inew] = mz2[i2];
            newBins[inew] = bins2[i2];
            inew++;
          }
          i2++;
        }
      }
    }
  }
  
  //copy used elements
  NumericVector resMass(inew);
  NumericVector resBins(inew);
  memcpy(resMass.begin(), newMz,  sizeof(double)*inew);
  memcpy(resBins.begin(), newBins,  sizeof(double)*inew);
  delete[] newMz;
  delete[] newBins;
  
  return List::create( Named("mass") = resMass, Named("bins") = resBins );
}
