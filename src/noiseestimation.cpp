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
#include "noiseestimation.h"
using namespace Rcpp;


NoiseEstimation::NoiseEstimation(int dataLength)
{
  //Init FFT objects according dataLength
  FFT_Size = (int)pow(2.0, std::ceil(log2(dataLength)));
  fft_in = fftw_alloc_real(FFT_Size);
  fft_out = fftw_alloc_real(FFT_Size);
  fft_pdirect = fftw_plan_r2r_1d(FFT_Size, fft_in, fft_out, FFTW_R2HC, FFTW_ESTIMATE);
  fft_pinvers = fftw_plan_r2r_1d(FFT_Size, fft_out, fft_in, FFTW_HC2R, FFTW_ESTIMATE);
  filWin = new double[1+FFT_Size/2];
  filWinMode = none;
  filWinSize = 0;
}

NoiseEstimation::~NoiseEstimation()
{
  fftw_destroy_plan(fft_pdirect);
  fftw_destroy_plan(fft_pinvers);
  fftw_free(fft_in); 
  fftw_free(fft_out);
  delete[] filWin;
}

void NoiseEstimation::NoiseEstimationFFTCosWin( double *data, int dataLength, int WinSize)
{
  if( filWinSize != WinSize ||  filWinMode != cos )
  {
    ComputeCosWin(WinSize);
  }
  
  NoiseEstimationFFT(data, dataLength);
}

void NoiseEstimation::NoiseEstimationFFTExpWin( double *data, int dataLength, int WinSize)
{
  if( filWinSize != WinSize ||  filWinMode != exp )
  {
    ComputeExpWin(WinSize);
  }
  
  NoiseEstimationFFT(data, dataLength);
}

void NoiseEstimation::NoiseEstimationFFT(double *data, int dataLength)
{
  
  //Copy data adding padding zeros
  double *xc_data = new double[FFT_Size];
  for( int i = 0; i < FFT_Size; i++)
  {
    xc_data[i] = i < dataLength? data[i] : 0.0;
  }

  if( filWinSize == 0 || filWinMode == none )
  {
    stop("Error: Filtering Windows has not been calculated yet");
    delete[] xc_data;
    return; 
  }
  
  //Copy data to a FFT objects
  memcpy(fft_in, xc_data, sizeof(double)*FFT_Size);

  //FFT data
  fftw_execute(fft_pdirect);
  
  //Apply the window function
  for( int i = 0; i <= FFT_Size/2; i++)
  {
    fft_out[i] *= filWin[i]; //The real part
    if(i > 0 && i < FFT_Size/2)
    {
      fft_out[ FFT_Size - i] *= filWin[i]; //The imaginary part
    } 
  }
  
  //The invers FFT
  fftw_execute(fft_pinvers);
  
  //Copy data from FFT object respecting original data size
  for( int i = 0; i < dataLength; i++)
  {
    data[i] = fft_in[i] / (double)FFT_Size; //Amplitude scaling to fit original range
  }
  
  delete[] xc_data;
}

NumericVector NoiseEstimation::NoiseEstimationFFTCosWin( NumericVector data, int WinSize )
{
  double *xc_data = new double[data.length()];
  memcpy(xc_data, data.begin(), sizeof(double)*data.length());
  NoiseEstimationFFTCosWin(xc_data, data.length(), WinSize);
  NumericVector y(data.length());
  memcpy(y.begin(), xc_data, sizeof(double)*data.length());
  delete[] xc_data;
  return y;
}

NumericVector NoiseEstimation::NoiseEstimationFFTExpWin( NumericVector data, int WinSize )
{
  double *xc_data = new double[data.length()];
  memcpy(xc_data, data.begin(), sizeof(double)*data.length());
  NoiseEstimationFFTExpWin(xc_data, data.length(), WinSize);
  NumericVector y(data.length());
  memcpy(y.begin(), xc_data, sizeof(double)*data.length());
  delete[] xc_data;
  return y;
}

void NoiseEstimation::ComputeCosWin(int WinSize)
{
  filWinSize = WinSize;
  for( int i = 0; i <= FFT_Size/2; i++)
  {
    filWin[i] = i < filWinSize ? (0.5 + 0.5*std::cos((2.0*M_PI*((double)i + 1.0))/(2.0 * (double)filWinSize))) : 0.0;
  }
  filWinMode = cos;
}

void NoiseEstimation::ComputeExpWin(int WinSize)
{
  filWinSize = WinSize;
  for( int i = 0; i <= FFT_Size/2; i++)
  {
    filWin[i] = std::exp(-5.0*i/(filWinSize));
  }
  filWinMode = exp;
}

int NoiseEstimation::getFFTSize()
{
  return FFT_Size;
}


////// Rcpp Exported methods //////////////////////////////////////////////////////////
//' NoiseEstimationFFTCosWin.
//' 
//' Estimate the noise of a spectrum using a FFT filter and a cosinus window in frequency domain.
//' 
//' @param x an Rcpp::NumericVector containing the spectrum intensities.
//' @param filWinSize an integer specified the cosinus win size in samples.
//' 
//' @return an Rcpp::NumericVector containing the estimated noise.
//' @export
// [[Rcpp::export]]
NumericVector NoiseEstimationFFTCosWin(NumericVector x,  int filWinSize = 40)
{
  NoiseEstimation neObj(x.length());
  return neObj.NoiseEstimationFFTCosWin(x, filWinSize);
}

//' NoiseEstimationFFTExpWin.
//' 
//' Estimate the noise of a spectrum using a FFT filter and a decay exponential window in frequency domain.
//' 
//' @param x an Rcpp::NumericVector containing the spectrum intensities.
//' @param filWinSize an integer specified the cosinus win size in samples.
//' 
//' @return an Rcpp::NumericVector containing the estimated noise.
//' @export
// [[Rcpp::export]]
NumericVector NoiseEstimationFFTExpWin(NumericVector x,  int filWinSize = 40)
{
  NoiseEstimation neObj(x.length());
  return neObj.NoiseEstimationFFTExpWin(x, filWinSize);
}

//' NoiseEstimationFFTCosWinMat.
//' 
//' Estimate the noise of some spectra using a FFT filter and a cosinus window in frequency domain.
//' 
//' @param x an Rcpp::NumericMatrix containing the spectra intensities. Each spectrum in a row.
//' @param filWinSize an integer specified the cosinus win size in samples.
//' 
//' @return an Rcpp::NumericMatrix containing the estimated noise in a matrix where each spectrum is a row.
//' @export
// [[Rcpp::export]]
NumericMatrix NoiseEstimationFFTCosWinMat(NumericMatrix x,  int filWinSize = 40)
{
  NoiseEstimation neObj(x.cols());
  NumericMatrix y(x.rows(), x.cols());
  for( int i = 0; i < x.rows(); i++)
  {
    y.row(i) = neObj.NoiseEstimationFFTCosWin(x.row(i), filWinSize);
  }
  return y;
}

//' NoiseEstimationFFTExpWinMat.
//' 
//' Estimate the noise of some spectra using a FFT filter and a decay exponential window in frequency domain.
//' 
//' @param x an Rcpp::NumericMatrix containing the spectra intensities. Each spectrum in a row.
//' @param filWinSize an integer specified the cosinus win size in samples.
//' 
//' @return an Rcpp::NumericMatrix containing the estimated noise in a matrix where each spectrum is a row.
//' @export
// [[Rcpp::export]]
NumericMatrix NoiseEstimationFFTExpWinMat(NumericMatrix x,  int filWinSize = 40)
{
  NoiseEstimation neObj(x.cols());
  NumericMatrix y(x.rows(), x.cols());
  for( int i = 0; i < x.rows(); i++)
  {
    y.row(i) = neObj.NoiseEstimationFFTExpWin(x.row(i), filWinSize);
  }
  return y;
}
