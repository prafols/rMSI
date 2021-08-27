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
#ifndef PEAK_PICKING_H
  #define PEAK_PICKING_H

#include <Rcpp.h>
#include <fftw3.h>
#include <vector>
#include "noiseestimation.h"

class PeakPicking
{
  public:
    PeakPicking(int WinSize, double *massAxis, int numOfDataPoints, int UpSampling = 10 );
    ~PeakPicking();
    Rcpp::NumericVector getHannWin();

    //Peak data structure
    typedef struct
    {
      std::vector<double> mass;
      std::vector<double> intensity;
      std::vector<double> SNR;
      std::vector<double> area;
      std::vector<double> binSize;
    } Peaks;
    
    //Detect all local maximums and Filter peaks using SNR min value in a sliding window
    //Returns a pointer to a Peaks structur, is programer responsability to free memory of returned data pointer.
    Peaks *peakPicking(double *spectrum, double SNR = 5);
    
    Rcpp::List PeakObj2List(PeakPicking::Peaks *pks);
    
    //Function to test interpolations
    Rcpp::NumericVector getInterpolatedPeak(Rcpp::NumericVector data, int iPeak, bool ApplyHanning);
    
    //Function to test Hanning Window implementation
    Rcpp::NumericVector getHanningWindow();
    Rcpp::NumericVector getAreaWindow();
    double *mass;
    
    double predictPeakArea( double *spectrum, int iPeakMass ); //Public because is used by zeroremover
    
  private:
    int FFT_Size; //First FFT windows size
    int FFTInter_Size; //Second FFT Windows size, used for peak interpolation
    double *HanningWin; //Hanning function computed in constructor for acurate peak mass prediction.
    double *AreaWin; //Wide window function to aboid Gibbs phenomenon on area interpolation.
    int dataLength;
    
    NoiseEstimation *neObj;
    
    //Data for FFT interpolation
    double *fft_in1; //Used for first fft input buffer with a length off FFT_Size
    double *fft_out1; //Used for first fft output buffer with a length off FFT_Size
    double *fft_in2; //Used for second fft input buffer with a length off FFTInter_Size
    double *fft_out2; //Used for second fft output buffer with a length off FFTInter_Size
    fftw_plan fft_pdirect;
    fftw_plan fft_pinvers;
    
    double predictPeakMass( double *spectrum, int iPeakMass );
    int interpolateFFT(double *spectrum, int iPeakMass, bool ApplyHanning); //Interpolatea m/z peak and return its location in interpolated space (fft_out2)
    Peaks *detectPeaks( double *spectrum, double *noise, double SNR );
};
  
#endif  