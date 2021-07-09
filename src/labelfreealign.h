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

#ifndef LABEL_FREE_ALIGN_H
  #define LABEL_FREE_ALIGN_H

#include <Rcpp.h>
#include <fftw3.h>

#define BOTTOM_SPECTRUM 0
#define CENTER_SPECTRUM 1
#define TOP_SPECTRUM 2
//#define EXTRA_DEBUG_INFO //Used for debuggin, comment out!


class LabelFreeAlign
{
  public:
    //spectraSplit the low/high part of spectra to keep (the resting points to unit will be removed).
    LabelFreeAlign(double *mass, double *ref_spectrum, int numOfPoints,
                   bool bilinear, int iterations = 3, 
                   double lagRefLow = 0.1, double lagRefMid = 0.5, double lagRefHigh = 0.9,
                   double lagLimitppm = 200, int fftOverSampling = 10, double winSizeRelative = 0.6);
    ~LabelFreeAlign();
    
    //Data accessors to internal vars to test this class
    Rcpp::NumericVector getHannWindow();
    Rcpp::NumericVector getHannWindowCenter();
    Rcpp::NumericVector getRefLowFFT();
    Rcpp::NumericVector getRefCenterFFT();
    Rcpp::NumericVector getRefHighFFT();
#ifdef EXTRA_DEBUG_INFO
    Rcpp::NumericVector DBG_getRefLow_fftinbuffer();
    Rcpp::NumericVector DBG_getRefMid_fftinbuffer();
    Rcpp::NumericVector DBG_getRefHIGH_fftinbuffer();
    Rcpp::NumericVector DBG_getSignalLow_fftinbuffer();
    Rcpp::NumericVector DBG_getSignalMid_fftinbuffer();
    Rcpp::NumericVector DBG_getSignalHigh_fftinbuffer();
#endif
    
    typedef struct
    {
      double lagHigh;
      double lagMid;
      double lagLow;
    }TLags;
    
    //Algin a given spectrum to the reference specified in the constructor. 
    //The used High and Low lags are returned as a TLags structure.
    // Arguments:
    // - intensityDataInterpolated: a pointer to the intensity spectrum. If data is in processed mode this argument corresponds to the interpolated spectrum so its length is not N.
    // -massData: a pointer to the mass axis independently if data is in processed or continuous mode (as is in the imzML file)
    // -intensityData: a pointer to the mass axis independently if data is in processed or continuous mode (as is in the imzML file)
    // - N: number of mass channels in the spectrum massData. Set to zero for continous data mode.
    TLags AlignSpectrum(double *intensityDataInterpolated, double *massData, double *intensityData, int N);

  private:
    void ComputeRef(double *data_ref, int spectrumPart);
    void ZeroPadding(double *data,bool reverse, int targetSize, int dataSize);
    void CopyData2Window(double *data_int, double *data_out,  int spectrumPart);
    void TimeWindow(double *data, int spectrumPart);
    double FourierBestCor(double *data, double *ref);
    
    int dataLength; //Number of points used in each spectrum
    double *commonMassAxis; //The common mass axis for the whole dataset
    int WinLength; //Number of points of spectrum retained in hanning window
    int FFT_Size_direct; //Number of points used for fft direct
    int FFT_Size_inverse; //Number of points used for fft inverse (which is diferent than direct to allow interpolation for lag values)
    
    fftw_plan fft_pdirect;
    fftw_plan fft_pinvers;
    
    double *fft_direct_in;
    double *fft_direct_out;
    double *fft_inverse_in;
    double *fft_inverse_out;
    
    //Mem space to store pre-computed reference FFT space values
    double *fft_ref_low;
    double *fft_ref_center;
    double *fft_ref_high;
    double *HannWindow;
    double *HannWindowCenter;
    int refMassLowIndex;
    int refMassMidIndex;
    int refMassHighIndex;
    double lagMaxppm;
    bool bBilinear;
    int AlignIterations;
    int FFTInterpolationOverSampling;
    
#ifdef EXTRA_DEBUG_INFO
   Rcpp::NumericVector DBG_refLow_fftInBuffer;
   Rcpp::NumericVector DBG_refMid_fftInBuffer;
   Rcpp::NumericVector DBG_refHigh_fftInBuffer;
   Rcpp::NumericVector DBG_signalLow_fftInBuffer;
   Rcpp::NumericVector DBG_signalMid_fftInBuffer;
   Rcpp::NumericVector DBG_signalHigh_fftInBuffer;
#endif
    
};

#endif
