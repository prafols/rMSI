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

#ifndef NOISE_ESTIMATION_H
  #define NOISE_ESTIMATION_H
  
#include <Rcpp.h>
#include <fftw3.h>

class NoiseEstimation
{
  public:
    NoiseEstimation(int dataLength);
    ~NoiseEstimation();
    void NoiseEstimationFFTCosWin( double *data, int dataLength, int WinSize );
    void NoiseEstimationFFTExpWin( double *data, int dataLength, int WinSize );
    Rcpp::NumericVector NoiseEstimationFFTCosWin( Rcpp::NumericVector data, int WinSize );
    Rcpp::NumericVector NoiseEstimationFFTExpWin( Rcpp::NumericVector data, int WinSize );
    int getFFTSize();
    
    
  private:
    int FFT_Size;
    double *fft_in;
    double *fft_out;
    double *filWin; //Windows used in FFT space for noise estimation
    enum FilWinType {none, cos, exp};
    FilWinType filWinMode;
    int filWinSize;
    fftw_plan fft_pdirect;
    fftw_plan fft_pinvers;
    
    void ComputeCosWin(int WinSize);
    void ComputeExpWin(int WinSize);
    void NoiseEstimationFFT(double *data, int dataLength);
};
  
#endif