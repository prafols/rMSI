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
#include <vector>
#include "smoothing.h"
#include "mlinterp.hpp" //Used for linear interpolation
using namespace Rcpp;

Smoothing::Smoothing(int kernelSize )
{
  switch(kernelSize)
  {
  case 5:
    sgC.push_back(-3.0/35.0);
    sgC.push_back(12.0/35.0);
    sgC.push_back(17.0/35.0);
    sgC.push_back(12.0/35.0);
    sgC.push_back(-3.0/35.0);
    break;
  case 7:
    sgC.push_back(-2.0/21.0);
    sgC.push_back(3.0/21.0);
    sgC.push_back(6.0/21.0);
    sgC.push_back(7.0/21.0);
    sgC.push_back(6.0/21.0);
    sgC.push_back(3.0/21.0);
    sgC.push_back(-2.0/21.0);
    break;
  case 9:
    sgC.push_back(-21.0/231.0);
    sgC.push_back(14.0/231.0);
    sgC.push_back(39.0/231.0);
    sgC.push_back(54.0/231.0);
    sgC.push_back(59.0/231.0);
    sgC.push_back(54.0/231.0);
    sgC.push_back(39.0/231.0);
    sgC.push_back(14.0/231.0);
    sgC.push_back(-21.0/231.0);
    break;
  case 11:
    sgC.push_back(-36.0/429.0);
    sgC.push_back(9.0/429.0);
    sgC.push_back(44.0/429.0);
    sgC.push_back(69.0/429.0);
    sgC.push_back(84.0/429.0);
    sgC.push_back(89.0/429.0);
    sgC.push_back(84.0/429.0);
    sgC.push_back(69.0/429.0);
    sgC.push_back(44.0/429.0);
    sgC.push_back(9.0/429.0);
    sgC.push_back(-36.0/429.0);
    break;
  case 13:
    sgC.push_back(-11.0/143.0);
    sgC.push_back(0.0/143.0);
    sgC.push_back(9.0/143.0);
    sgC.push_back(16.0/143.0);
    sgC.push_back(21.0/143.0);
    sgC.push_back(24.0/143.0);
    sgC.push_back(25.0/143.0);
    sgC.push_back(24.0/143.0);
    sgC.push_back(21.0/143.0);
    sgC.push_back(16.0/143.0);
    sgC.push_back(9.0/143.0);
    sgC.push_back(0.0/143.0);
    sgC.push_back(-11.0/143.0);
    break;
  case 15:
    sgC.push_back(-78.0/1105.0);
    sgC.push_back(-13.0/1105.0);
    sgC.push_back(42.0/1105.0);
    sgC.push_back(87.0/1105.0);
    sgC.push_back(122.0/1105.0);
    sgC.push_back(147.0/1105.0);
    sgC.push_back(162.0/1105.0);
    sgC.push_back(167.0/1105.0);
    sgC.push_back(162.0/1105.0);
    sgC.push_back(147.0/1105.0);
    sgC.push_back(122.0/1105.0);
    sgC.push_back(87.0/1105.0);
    sgC.push_back(42.0/1105.0);
    sgC.push_back(-13.0/1105.0);
    sgC.push_back(-78.0/1105.0);
    break;
  default:
    stop("Error, not valid SavitzkyGolay kernel size, valid values are: 5, 7, 9, 11, 13, 15");
  }
  
}

Smoothing::~Smoothing()
{
  
}

NumericVector Smoothing::smoothSavitzkyGolay(NumericVector x)
{
  double *xC = new double[x.length()];
  memcpy(xC, x.begin(), sizeof(double)*x.length());
  smoothSavitzkyGolay(xC, x.length());
  NumericVector y(x.length());
  memcpy(y.begin(), xC, sizeof(double)*y.length());
  delete[] xC;
  return y;
}

//Performs the SavitzkyGolay smoothing and overwites the original data with smoothed data
void Smoothing::smoothSavitzkyGolay(double *x, int length)
{
  //Convolution with SavitzkyGolay kernel
  double *y = new double[length];
  for( int i = 0; i < length; i++)
  {
    y[i] = 0; //Init a zero
  }
  for( int i = (sgC.size() - 1)/2; i < length - ((sgC.size() - 1)/2); i++)
  {
    for( int j = 0; j < sgC.size(); j++)
    {
      y[i] = y[i] + x[i + j -((sgC.size() - 1)/2) ] * sgC[j];
    }
  }
  
  //Overwrite input pointer
  memcpy(x, y, sizeof(double)*length);
  delete[] y;
}

//Smooth a given spectrum and interpolate it to imzML processed mode. 
// Arguments:
// - commonMassAxis: a pointer to the common mass axis used in the interpolated spectrum.
// - intensityDataInterpolated: a pointer to the intensity spectrum. If data is in processed mode this argument corresponds to the interpolated spectrum.
// - length: length of the interpolated intensity vector.
// -massData: a pointer to the mass axis independently if data is in processed or continuous mode (as is in the imzML file)
// -intensityData: a pointer to the mass axis independently if data is in processed or continuous mode (as is in the imzML file)
// - N: number of mass channels in the spectrum massData. Set to zero for continous data mode.
void Smoothing::smoothSavitzkyGolay(double *commonMassAxis, double *intensityDataInterpolated, int length, double *massData, double *intensityData, int N)
{
  //Smooth the interpolated intensity spectrum
  smoothSavitzkyGolay(intensityDataInterpolated, length);
  
  //Interpolate to the common mass axis
  if(N > 0)
  {
    mlinterp::interp(
      &length, N, // Number of points (interpolated --> imzML original)
      intensityDataInterpolated, intensityData, // Y axis  (interpolated --> imzML original)
      commonMassAxis, massData // X axis  (interpolated --> imzML original)
    );
  }
}

//' Smoothing_SavitzkyGolay.
//' 
//' Computes the Savitzky-Golay smoothing of a vector x using a filter size of sgSize.
//' @param x the data vector to smooth.
//' @param sgSize valid values are: 5, 7, 9, 11, 13, 15.
//' @return the smoothed data vector.
// [[Rcpp::export]]
NumericVector Smoothing_SavitzkyGolay(NumericVector x, int sgSize = 5)
{
  Smoothing smObj(sgSize);
  return smObj.smoothSavitzkyGolay(x);
}


