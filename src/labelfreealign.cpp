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
#include "labelfreealign.h"
#include "mlinterp.hpp" //Used for linear interpolation

using namespace Rcpp;

LabelFreeAlign::LabelFreeAlign(double *mass, double *ref_spectrum, int numOfPoints,
                               bool bilinear, int iterations, 
                               double lagRefLow, double lagRefMid, double lagRefHigh,
                               double lagLimitppm, int fftOverSampling, double winSizeRelative ):
dataLength(numOfPoints),
commonMassAxis(mass),
bBilinear(bilinear),
AlignIterations(iterations),
lagMaxppm(lagLimitppm),
FFTInterpolationOverSampling(fftOverSampling)

{
  WinLength = (int)round(winSizeRelative * (double)dataLength);
  
  //Pre-compute the Hanning Windows
  HannWindow = new double[WinLength];
  HannWindowCenter = new double[WinLength];
  for(int i = 0; i < WinLength; i++)
  {
    HannWindow[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( 2.0 * WinLength ) ));
    HannWindowCenter[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( WinLength ) ));
  }
  
  FFT_Size_direct = (int)pow(2.0, std::ceil(log2(WinLength)));
  FFT_Size_inverse = (int)pow(2.0, std::ceil(log2(WinLength*FFTInterpolationOverSampling)));
  fft_direct_in = fftw_alloc_real(FFT_Size_direct);
  fft_direct_out = fftw_alloc_real(FFT_Size_direct);
  fft_inverse_in = fftw_alloc_real(FFT_Size_inverse);
  fft_inverse_out = fftw_alloc_real(FFT_Size_inverse);
  
  fft_pdirect = fftw_plan_r2r_1d(FFT_Size_direct, fft_direct_in, fft_direct_out, FFTW_R2HC, FFTW_ESTIMATE);
  fft_pinvers = fftw_plan_r2r_1d(FFT_Size_inverse, fft_inverse_in, fft_inverse_out, FFTW_HC2R, FFTW_ESTIMATE);
  
  fft_ref_low = new  double[FFT_Size_inverse];
  fft_ref_center = new  double[FFT_Size_inverse];
  fft_ref_high = new  double[FFT_Size_inverse];

  //Lag limits realtive to mass references
  refMassLowIndex = (int)round(lagRefLow * (double)dataLength);
  refMassMidIndex = (int)round(lagRefMid * (double)dataLength);
  refMassHighIndex = (int)round(lagRefHigh * (double)dataLength);
    
  ComputeRef(ref_spectrum, BOTTOM_SPECTRUM);
  ComputeRef(ref_spectrum, CENTER_SPECTRUM);
  ComputeRef(ref_spectrum, TOP_SPECTRUM);

#ifdef EXTRA_DEBUG_INFO 
  Rcpp::Rcout<<"\n===============================================\n";
  Rcpp::Rcout<<"DBG: refMassLowIndex = "<<refMassLowIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassMidIndex = "<<refMassMidIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassHighIndex = "<<refMassHighIndex<<"\n";
  Rcpp::Rcout<<"===============================================\n";
#endif
  
  refMassLowIndex = refMassLowIndex < 0 ? 0 : refMassLowIndex;
  refMassLowIndex = refMassLowIndex > (dataLength - 1) ? (dataLength - 1) : refMassLowIndex;
  
  refMassMidIndex = refMassMidIndex < 0 ? 0 : refMassMidIndex;
  refMassMidIndex = refMassMidIndex > (dataLength - 1) ? (dataLength - 1) : refMassMidIndex;
  
  refMassHighIndex = refMassHighIndex < 0 ? 0 : refMassHighIndex;
  refMassHighIndex = refMassHighIndex > (dataLength - 1) ? (dataLength - 1) : refMassHighIndex;
  
#ifdef EXTRA_DEBUG_INFO 
  Rcpp::Rcout<<"\n===============================================\n";
  Rcpp::Rcout<<"DBG: refMassLowIndex = "<<refMassLowIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassMidIndex = "<<refMassMidIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassHighIndex = "<<refMassHighIndex<<"\n";
  Rcpp::Rcout<<"===============================================\n";
#endif
  
}

LabelFreeAlign::~LabelFreeAlign()
{
  fftw_destroy_plan(fft_pdirect);
  fftw_destroy_plan(fft_pinvers);
  fftw_free(fft_direct_in); 
  fftw_free(fft_direct_out);
  fftw_free(fft_inverse_in); 
  fftw_free(fft_inverse_out);
  delete[] fft_ref_low;
  delete[] fft_ref_center;
  delete[] fft_ref_high;
  delete[] HannWindow;
  delete[] HannWindowCenter;
}

void LabelFreeAlign::ComputeRef(double *data_ref, int spectrumPart)
{
  double *data_ptr;
  switch(spectrumPart)
  {
    case BOTTOM_SPECTRUM:
      data_ptr = fft_ref_low; 
      break;
      
    case CENTER_SPECTRUM:
      data_ptr = fft_ref_center; 
      break;
      
    case TOP_SPECTRUM:
      data_ptr = fft_ref_high;
      break;
  }
  
  CopyData2Window(data_ref, fft_direct_in, spectrumPart);
  TimeWindow(fft_direct_in,  spectrumPart);
  ZeroPadding(fft_direct_in, spectrumPart == TOP_SPECTRUM, FFT_Size_direct, WinLength);
  
#ifdef EXTRA_DEBUG_INFO
  switch(spectrumPart)
  {
  case BOTTOM_SPECTRUM:
    DBG_refLow_fftInBuffer = Rcpp::NumericVector(FFT_Size_direct); 
    memcpy(DBG_getRefLow_fftinbuffer().begin(), fft_direct_in, sizeof(double)*FFT_Size_direct);
    break;
    
  case CENTER_SPECTRUM:
    DBG_refMid_fftInBuffer = Rcpp::NumericVector(FFT_Size_direct); 
    memcpy(DBG_refMid_fftInBuffer.begin(), fft_direct_in, sizeof(double)*FFT_Size_direct);
    break;
    
  case TOP_SPECTRUM:
    DBG_refHigh_fftInBuffer = Rcpp::NumericVector(FFT_Size_direct); 
    memcpy(DBG_refHigh_fftInBuffer.begin(), fft_direct_in, sizeof(double)*FFT_Size_direct);
    break;
  }
#endif  
  
  fftw_execute(fft_pdirect);
  
  //FFT domain Interpolation
  if(FFT_Size_direct < FFT_Size_inverse)
  {
    memcpy(data_ptr, fft_direct_out, sizeof(double)*(1+(FFT_Size_direct/2)));  //Fill the real part
    memcpy(data_ptr + (FFT_Size_inverse - (FFT_Size_direct/2) + 1), fft_direct_out + (1+(FFT_Size_direct/2)), sizeof(double)*((FFT_Size_direct/2)-1)); //Fill the imaginary part
    for(int i = (1+(FFT_Size_direct/2)); i < (FFT_Size_inverse - (FFT_Size_direct/2) + 1); i++) data_ptr[i] = 0.0; //Zero padding
  }
  else
  {
    //No interpolation... just copy
    memcpy(data_ptr, fft_direct_out, sizeof(double)*FFT_Size_direct);
  }
  
  //Compute Conj
  for( int i = ((FFT_Size_inverse/2) + 1) ; i < FFT_Size_inverse; i++)
  {
    data_ptr[i] *= (-1.0); 
  }
}

void LabelFreeAlign::ZeroPadding(double *data,bool reverse, int targetSize, int dataSize)
{
  if(dataSize == targetSize)
  {
    //No padding needed
    return;
  }
  if(reverse)
  {
    for( int i = 0; i < (targetSize - dataSize); i++)
    {
      data[i] = 0.0;
    }
  }
  else
  {
    for( int i = dataSize; i < targetSize; i++)
    {
      data[i] = 0.0;
    }
  }
}

void LabelFreeAlign::CopyData2Window(double *data_int, double *data_out, int spectrumPart)
{
  int offset_in;
  int offset_out;
  switch(spectrumPart)
  {
    case BOTTOM_SPECTRUM:
      offset_in = refMassLowIndex;
      offset_out = 0;
      break;
      
    case CENTER_SPECTRUM:
      offset_in = refMassMidIndex - WinLength/2;
      offset_out = 0;
      break;
      
    case TOP_SPECTRUM:
      offset_in = refMassHighIndex - WinLength;
      offset_out = FFT_Size_direct - WinLength;
      break;
  }
  
  //Move pointer offsets to avoid a segfault accesing outside vector
  if( (offset_in + WinLength) >  dataLength)
  {
    offset_in -= (offset_in + WinLength) - dataLength;
  }
  if( offset_in < 0 )
  {
    offset_in = 0;
  }
  
  memcpy(data_out + offset_out, data_int + offset_in, sizeof(double)*WinLength);
}

void LabelFreeAlign::TimeWindow(double *data,  int spectrumPart)
{
  for( int i = 0; i < WinLength; i++)
  {
    
    switch(spectrumPart)
    {
    case BOTTOM_SPECTRUM:
      data[i] *= HannWindow[WinLength - i - 1]; 
      break;
      
    case CENTER_SPECTRUM:
      data[i] *= HannWindowCenter[i];
      break;
      
    case TOP_SPECTRUM:
      data[i + (FFT_Size_direct - WinLength)] *= HannWindow[i]; 
      break;
    }
  }
}

LabelFreeAlign::TLags LabelFreeAlign::AlignSpectrum(double *intensityDataInterpolated, double *massData, double *intensityData, int N)
{
  //Hanning Windowing
  double *topWin_data = new double[FFT_Size_direct];
  double *midWin_data = new double[FFT_Size_direct];
  double *botWin_data = new double[FFT_Size_direct];
  TLags firstLag;
  
  //Prepare data pointer for continuous mode:
  double *ptrMass, *ptrIntensity; 
  bool bDataInContinuousMode;
  if(N == 0)
  {
    //imzML in continuous mode
    ptrMass = new double[dataLength];   
    ptrIntensity = new double[dataLength]; 
    memcpy(ptrMass, commonMassAxis, sizeof(double)*dataLength); //Copy the common mass axis to the realigned mass axis
    memcpy(ptrIntensity, intensityDataInterpolated, sizeof(double)*dataLength); //Copy the intensity to the realigned intensity
    N = dataLength;
    bDataInContinuousMode = true;
  }
  else
  {
    //imzML in processed mode 
    ptrMass = massData;
    ptrIntensity = intensityData;
    bDataInContinuousMode = false;
  }
  
  for(int i = 0; i < AlignIterations; i++) 
  {
    CopyData2Window(intensityDataInterpolated, topWin_data, TOP_SPECTRUM);
    CopyData2Window(intensityDataInterpolated, midWin_data, CENTER_SPECTRUM);
    CopyData2Window(intensityDataInterpolated, botWin_data, BOTTOM_SPECTRUM);
    
    TimeWindow(topWin_data, TOP_SPECTRUM);
    TimeWindow(midWin_data, CENTER_SPECTRUM);
    TimeWindow(botWin_data, BOTTOM_SPECTRUM);
    
    //Zero-padding 2 improve fft performance
    ZeroPadding(topWin_data, true, FFT_Size_direct, WinLength);
    ZeroPadding(botWin_data, false, FFT_Size_direct, WinLength);
    ZeroPadding(midWin_data, false, FFT_Size_direct, WinLength);
    
#ifdef EXTRA_DEBUG_INFO
    DBG_signalLow_fftInBuffer = Rcpp::NumericVector(FFT_Size_direct); 
    memcpy(DBG_signalLow_fftInBuffer.begin(), botWin_data, sizeof(double)*FFT_Size_direct);

    DBG_signalMid_fftInBuffer = Rcpp::NumericVector(FFT_Size_direct); 
    memcpy(DBG_signalMid_fftInBuffer.begin(), midWin_data, sizeof(double)*FFT_Size_direct);
    
    DBG_signalHigh_fftInBuffer = Rcpp::NumericVector(FFT_Size_direct); 
    memcpy(DBG_signalHigh_fftInBuffer.begin(), topWin_data, sizeof(double)*FFT_Size_direct);
#endif  
    
    //Get lags
    TLags lags;
    lags.lagLow = FourierBestCor(botWin_data, fft_ref_low);
    lags.lagMid = FourierBestCor(midWin_data, fft_ref_center);
    lags.lagHigh = FourierBestCor(topWin_data, fft_ref_high);

#ifdef EXTRA_DEBUG_INFO    
    Rcpp::Rcout<<"\n===============================================\n";
    Rcpp::Rcout<<"DBG: lags.lagLow = "<<lags.lagLow<<"\n";
    Rcpp::Rcout<<"DBG: lags.lagMid = "<<lags.lagMid<<"\n";
    Rcpp::Rcout<<"DBG: lags.lagHigh = "<<lags.lagHigh<<"\n";
    Rcpp::Rcout<<"===============================================\n";
#endif
    
    //Limit lag Low
    int targetMassLowIndex = refMassLowIndex + lags.lagLow;
    targetMassLowIndex = targetMassLowIndex >= dataLength ? dataLength - 1 : targetMassLowIndex;
    targetMassLowIndex = targetMassLowIndex < 0 ? 0 : targetMassLowIndex;
    double lagLowppm = 1e6*fabs(commonMassAxis[refMassLowIndex] - commonMassAxis[targetMassLowIndex])/commonMassAxis[refMassLowIndex];
    lags.lagLow = lagLowppm > lagMaxppm ? 0.0 : lags.lagLow;
    
    //Limit lag Mid
    int targetMassMidIndex = refMassMidIndex + lags.lagMid;
    targetMassMidIndex = targetMassMidIndex >= dataLength ? dataLength - 1 : targetMassMidIndex;
    targetMassMidIndex = targetMassMidIndex < 0 ? 0 : targetMassMidIndex;
    double lagMidppm = 1e6*fabs(commonMassAxis[refMassMidIndex] - commonMassAxis[targetMassMidIndex])/commonMassAxis[refMassMidIndex];
    lags.lagMid = lagMidppm > lagMaxppm ? 0.0 : lags.lagMid;
    
    //Limit lag High
    int targetMassHighIndex = refMassHighIndex + lags.lagHigh;
    targetMassHighIndex = targetMassHighIndex >= dataLength ? dataLength - 1 : targetMassHighIndex;
    targetMassHighIndex = targetMassHighIndex < 0 ? 0 : targetMassHighIndex;
    double lagHighppm = 1e6*fabs(commonMassAxis[refMassHighIndex] - commonMassAxis[targetMassHighIndex])/commonMassAxis[refMassHighIndex];
    lags.lagHigh = lagHighppm > lagMaxppm ? 0.0 : lags.lagHigh;

#ifdef EXTRA_DEBUG_INFO    
    Rcpp::Rcout<<"\n===============================================\n";
    Rcpp::Rcout<<"DBG: lagLowppm = "<<lagLowppm<<"\n";
    Rcpp::Rcout<<"DBG: lagMidppm = "<<lagMidppm<<"\n";
    Rcpp::Rcout<<"DBG: lagHighppm = "<<lagHighppm<<"\n";
    Rcpp::Rcout<<"===============================================\n";
    
    Rcpp::Rcout<<"\n===============================================\n";
    Rcpp::Rcout<<"DBG: lags.lagLow = "<<lags.lagLow<<"\n";
    Rcpp::Rcout<<"DBG: lags.lagMid = "<<lags.lagMid<<"\n";
    Rcpp::Rcout<<"DBG: lags.lagHigh = "<<lags.lagHigh<<"\n";
    Rcpp::Rcout<<"===============================================\n";
#endif
    
    if(i == 0)
    {
      firstLag = lags;
    }
  
    //Spectra warping constants
    double K1, K2, Sh1, Sh2;

    if(bBilinear)
    {
      //Shifting mass indexes in the common mass axis
      int indexLowMassShifted =(int)round(refMassLowIndex + lags.lagLow);
      int indexMidMassShifted =(int)round(refMassMidIndex + lags.lagMid);
      int indexHighMassShifted =(int)round(refMassHighIndex + lags.lagHigh);
      
      //Saturate to a valid range
      indexLowMassShifted = indexLowMassShifted >= dataLength ? (dataLength - 1) : indexLowMassShifted;
      indexLowMassShifted = indexLowMassShifted < 0 ? 0 : indexLowMassShifted;
      
      indexMidMassShifted = indexMidMassShifted >= dataLength ? (dataLength - 1) : indexMidMassShifted;
      indexMidMassShifted = indexMidMassShifted < 0 ? 0 : indexMidMassShifted;
      
      indexHighMassShifted = indexHighMassShifted >= dataLength ? (dataLength - 1) : indexHighMassShifted;
      indexHighMassShifted = indexHighMassShifted < 0 ? 0 : indexHighMassShifted;
      
      //Calculate mass shift and scaling constants
      K1 = (commonMassAxis[indexMidMassShifted] - commonMassAxis[indexLowMassShifted])/(commonMassAxis[refMassMidIndex] - commonMassAxis[refMassLowIndex]); 
      Sh1 = commonMassAxis[indexLowMassShifted] - commonMassAxis[refMassLowIndex]*K1; // y = mx + n --> n = y - mx
      
      K2 = (commonMassAxis[indexHighMassShifted] - commonMassAxis[indexMidMassShifted])/(commonMassAxis[refMassHighIndex] - commonMassAxis[refMassMidIndex]); 
      Sh2 = commonMassAxis[indexMidMassShifted] - commonMassAxis[refMassMidIndex]*K2; // y = mx + n --> n = y - mx
      
      //Apply the alignment to the imzML mass axis
      for(int i = 0; i < N; i++) 
      {
        if(ptrMass[i] < commonMassAxis[refMassMidIndex])
        {
          ptrMass[i] = ptrMass[i]*K1 + Sh1;
        }
        else
        {
          ptrMass[i] = ptrMass[i]*K2 + Sh2;
        }
      }
    }
    else
    {
        //Shifting mass indexes in the common mass axis
        int indexLowMassShifted = refMassLowIndex + lags.lagLow;
        int indexHighMassShifted = refMassHighIndex + lags.lagHigh;
        
        //Saturate to a valid range
        indexLowMassShifted = indexLowMassShifted >= dataLength ? (dataLength - 1) : indexLowMassShifted;
        indexLowMassShifted = indexLowMassShifted < 0 ? 0 : indexLowMassShifted;
        
        indexHighMassShifted = indexHighMassShifted >= dataLength ? (dataLength - 1) : indexHighMassShifted;
        indexHighMassShifted = indexHighMassShifted < 0 ? 0 : indexHighMassShifted;
        
        //Calculate mass shift and scaling constants
        K1 = (commonMassAxis[indexHighMassShifted] - commonMassAxis[indexLowMassShifted])/(commonMassAxis[refMassHighIndex] - commonMassAxis[refMassLowIndex]); 
        Sh1 = commonMassAxis[indexLowMassShifted] - commonMassAxis[refMassLowIndex]*K1; // y = mx + n --> n = y - mx
        
        //Apply the alignment to the imzML mass axis
        for(int i = 0; i < N; i++)
        {
          ptrMass[i] = ptrMass[i]*K1 + Sh1;
        }
    }
    
    //Interpolate to the common mass axis
    mlinterp::interp(
      &N, (int)dataLength, // Number of points (imzML original, interpolated )
      ptrIntensity, intensityDataInterpolated, // Y axis  (imzML original, interpolated )
      ptrMass, commonMassAxis // X axis  (imzML original, interpolated )
    );
  }
  
  delete[] topWin_data;
  delete[] botWin_data;
  delete[] midWin_data;
  if(bDataInContinuousMode)
  {
    delete[] ptrMass;
    delete[] ptrIntensity;
  }
    
  return firstLag;
}

double LabelFreeAlign::FourierBestCor(double *data, double *ref)
{
  double *fft_direct_out_interpolated = new double[FFT_Size_inverse];
  memcpy(fft_direct_in, data, sizeof(double)*FFT_Size_direct);
  fftw_execute(fft_pdirect);
  
  //FFT domain Interpolation
  if(FFT_Size_direct < FFT_Size_inverse)
  {
    memcpy(fft_direct_out_interpolated, fft_direct_out, sizeof(double)*(1+(FFT_Size_direct/2)));  //Fill the real part
    memcpy(fft_direct_out_interpolated + (FFT_Size_inverse - (FFT_Size_direct/2) + 1), fft_direct_out + (1+(FFT_Size_direct/2)), sizeof(double)*((FFT_Size_direct/2)-1)); //Fill the imaginary part
    for(int i = (1+(FFT_Size_direct/2)); i < (FFT_Size_inverse - (FFT_Size_direct/2) + 1); i++) fft_direct_out_interpolated[i] = 0.0; //Zero padding
  }
  else
  {
    //No interpolation... just copy
    memcpy(fft_direct_out_interpolated, fft_direct_out, sizeof(double)*FFT_Size_direct);
  }

  //Mult fft complex values, the ref is assumed already Conj (this is automatically done by ComputeRef method)
  for( int i = 0; i <= FFT_Size_inverse/2; i++)
  {
    if( i > 0 && i < FFT_Size_inverse/2)
    {
      fft_inverse_in[i] = ref[i] * fft_direct_out_interpolated[i] - ref[FFT_Size_inverse - i] * fft_direct_out_interpolated[FFT_Size_inverse - i];
      fft_inverse_in[FFT_Size_inverse - i] = ref[i] * fft_direct_out_interpolated[FFT_Size_inverse - i] + ref[FFT_Size_inverse - i] * fft_direct_out_interpolated[i];
    }
    else
    {
      fft_inverse_in[i] = ref[i] * fft_direct_out_interpolated[i];
    }
  }
  delete[] fft_direct_out_interpolated;
  
  fftw_execute(fft_pinvers);
  
  //Locate the max correlation
  double dMax = 0.0;
  int lag = 0;
  for( int i = 0; i < FFT_Size_inverse; i++)
  {
    if(fft_inverse_out[i] > dMax)
    {
      dMax = fft_inverse_out[i];
      lag = i;
    }
  }

#ifdef EXTRA_DEBUG_INFO   
  Rcpp::Rcout<<"\n===============================================\n";
  Rcpp::Rcout<<"DBG: lag RAW = "<<lag<<"\n";
  Rcpp::Rcout<<"DBG: FFT_Size_direct = "<<FFT_Size_direct<<"\n";
  Rcpp::Rcout<<"DBG: FFT_Size_inverse = "<<FFT_Size_inverse<<"\n";
  Rcpp::Rcout<<"===============================================\n";
#endif
  
  if( lag >= FFT_Size_inverse/2)
  {
    lag = FFT_Size_inverse - lag; 
  }
  else
  {
    lag = -lag;
  }
  

  return ((double)lag) * (((double)(FFT_Size_direct)) / ((double)(FFT_Size_inverse)) );
}

NumericVector LabelFreeAlign::getHannWindow()
{
  NumericVector hannWin(WinLength);
  memcpy(hannWin.begin(), HannWindow, sizeof(double)*WinLength);
  return hannWin;
}

NumericVector LabelFreeAlign::getHannWindowCenter()
{
  NumericVector hannWin(WinLength);
  memcpy(hannWin.begin(), HannWindowCenter, sizeof(double)*WinLength);
  return hannWin;
}

NumericVector LabelFreeAlign::getRefLowFFT()
{
  NumericVector refFft(FFT_Size_inverse);
  memcpy(refFft.begin(), fft_ref_low, sizeof(double)*FFT_Size_inverse);
  return refFft;
}

NumericVector LabelFreeAlign::getRefCenterFFT()
{
  NumericVector refFft(FFT_Size_inverse);
  memcpy(refFft.begin(), fft_ref_center, sizeof(double)*FFT_Size_inverse);
  return refFft;
}

NumericVector LabelFreeAlign::getRefHighFFT()
{
  NumericVector refFft(FFT_Size_inverse);
  memcpy(refFft.begin(), fft_ref_high, sizeof(double)*FFT_Size_inverse);
  return refFft;
}

#ifdef EXTRA_DEBUG_INFO
Rcpp::NumericVector LabelFreeAlign::DBG_getRefLow_fftinbuffer()
{
  return DBG_refLow_fftInBuffer;
}

Rcpp::NumericVector LabelFreeAlign::DBG_getRefMid_fftinbuffer()
{
  return DBG_refMid_fftInBuffer;
}

Rcpp::NumericVector LabelFreeAlign::DBG_getRefHIGH_fftinbuffer()
{
  return DBG_refHigh_fftInBuffer;
}

Rcpp::NumericVector LabelFreeAlign::DBG_getSignalLow_fftinbuffer()
{
  return DBG_signalLow_fftInBuffer;
}

Rcpp::NumericVector LabelFreeAlign::DBG_getSignalMid_fftinbuffer()
{
  return DBG_signalMid_fftInBuffer;
}

Rcpp::NumericVector LabelFreeAlign::DBG_getSignalHigh_fftinbuffer()
{
  return DBG_signalHigh_fftInBuffer;
}
#endif

//TESTS//////////////////////////////////////////////////////////////////////////////////////////////
/*
// [[Rcpp::export]]
List TestComputeRefAndHannWin(NumericVector refSpectrum)
{
  boost::mutex mtx;
  double refC[refSpectrum.length()];
  memcpy(refC, refSpectrum.begin(), sizeof(double)*refSpectrum.length());
  LabelFreeAlign alngObj(refC, refSpectrum.length(), false, &mtx);
  return List::create( Named("HannWin") = alngObj.getHannWindow(), Named("HannWinCenter") =  alngObj.getHannWindowCenter(), Named("RefLow") = alngObj.getRefLowFFT(), Named("RefCenter") = alngObj.getRefCenterFFT(),  Named("RefHigh") = alngObj.getRefHighFFT());
}
*/

/*
//To run this debug function the ZeroPadding method must be set as public
// [[Rcpp::export]]
NumericVector TestZeroPadding(NumericVector x, bool rev)
{
  boost::mutex mtx;
  const int NewLength = (int)pow(2.0, std::ceil(log2((double)x.length())));
  double refC[NewLength];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length(), false, &mtx);
  
  alngObj.ZeroPadding(refC, rev,  NewLength, x.length());
    
  NumericVector y(NewLength);
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}
*/

/*
//To run this debug function the TestTimeWindow method must be set as public
// [[Rcpp::export]]
NumericVector TestTimeWindow(NumericVector x, bool bHigh)
{
  boost::mutex mtx;
  double refC[x.length()];
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length(), false, &mtx);
  
  alngObj.TimeWindow(refC, bHigh);
  
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y; 
}
*/

/*
//To run this debug function the FourierLinerScaleShift method must be set as public
// [[Rcpp::export]]
NumericVector TestFourierLinerScaleShift(NumericVector x, double scaling,  double shift)
{
  double refC[x.length()];
  boost::mutex mtx;
  memcpy(refC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, x.length(), false, &mtx);
  
  alngObj.FourierLinerScaleShift(refC, scaling, shift);
  
  NumericVector y(x.length());
  memcpy(y.begin(), refC, sizeof(double)*x.length());
  return y;
}
*/

/*
//To run this debug function theFourierBestCormethod must be set as public
// [[Rcpp::export]]
double TestFourierBestCor(NumericVector ref, NumericVector x, bool bRefLow)
{
  boost::mutex mtx;
  double refC[ref.length()];
  double xC[x.length()];
  memcpy(refC, ref.begin(), sizeof(double)*ref.length());
  memcpy(xC, x.begin(), sizeof(double)*x.length());
  LabelFreeAlign alngObj(refC, ref.length(), false, &mtx);
  
  NumericVector ref_ptr;
  if(bRefLow)
  {
    ref_ptr = alngObj.getRefLowFFT();
  }
  else
  {
    ref_ptr = alngObj.getRefHighFFT();
  }
  
  return (double)alngObj.FourierBestCor(xC, ref_ptr.begin());
}
*/

//To run for debug purposes
// [[Rcpp::export]]
Rcpp::List AlignSpectrumToReference( NumericVector mass, NumericVector ref, NumericVector spectrumInterpolated,
                                        NumericVector massProcessedMode, NumericVector intensityProcessedMode, bool bilinear = false, 
                                        double lagRefLow = 0.1, double lagRefMid = 0.5, double lagRefHigh = 0.9,
                                        int iterations = 1, double lagLimitppm = 200, int fftOverSampling = 10, double winSizeRelative = 0.6 )
{
  NumericVector y(spectrumInterpolated.length());
  memcpy(y.begin(), spectrumInterpolated.begin(), sizeof(double)*spectrumInterpolated.length());
  
  LabelFreeAlign alngObj(mass.begin(), ref.begin(), ref.length(), 
                         bilinear, iterations, 
                         lagRefLow, lagRefMid, lagRefHigh, lagLimitppm, fftOverSampling, winSizeRelative);
  
  Rcpp::NumericVector refLow = alngObj.getRefLowFFT();
  Rcpp::NumericVector refMid = alngObj.getRefCenterFFT();
  Rcpp::NumericVector refHigh = alngObj.getRefHighFFT();
  
  LabelFreeAlign::TLags lags = alngObj.AlignSpectrum(y.begin(), massProcessedMode.begin(), intensityProcessedMode.begin(), massProcessedMode.length());
  Rcout<<"Lag low = "<<lags.lagLow<<" Lag center = "<<lags.lagMid<<" Lag high = "<<lags.lagHigh<<"\n";
  
#ifdef EXTRA_DEBUG_INFO
  return List::create( Named("InterpolatedSpectrum") = y, 
                       Named("ProcessedModeMass") = massProcessedMode, 
                       Named("RefLowFFT") = refLow,
                       Named("RefMidFFT") = refMid,
                       Named("RefHighFFT") = refHigh,
                       Named("RefLowFFTinbuffer") =alngObj.DBG_getRefLow_fftinbuffer(),
                       Named("RefMidFFTinbuffer") =alngObj.DBG_getRefMid_fftinbuffer(),
                       Named("RefHighFFTinbuffer") =alngObj.DBG_getRefHIGH_fftinbuffer(),
                       Named("SignalLowFFTinbuffer") =alngObj.DBG_getSignalLow_fftinbuffer(),
                       Named("SignalMidFFTinbuffer") =alngObj.DBG_getSignalMid_fftinbuffer(),
                       Named("SignalHighFFTinbuffer") =alngObj.DBG_getSignalHigh_fftinbuffer()
                       );
#else
  return List::create( Named("InterpolatedSpectrum") = y, 
                       Named("ProcessedModeMass") = massProcessedMode, 
                       Named("RefLowFFT") = refLow,
                       Named("RefMidFFT") = refMid,
                       Named("RefHighFFT") = refHigh
                       );
#endif
}


