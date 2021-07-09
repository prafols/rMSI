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

using namespace Rcpp;

LabelFreeAlign::LabelFreeAlign(double *mass, double *ref_spectrum, int numOfPoints,
                               bool bilinear, std::mutex *sharedMutex, int iterations, 
                               double lagRefLow, double lagRefMid, double lagRefHigh,
                               double lagLimitppm, int fftOverSampling, double winSizeRelative ):
dataLength(numOfPoints),
commonMassAxis(mass),
bBilinear(bilinear),
AlignIterations(iterations),
lagMaxppm(lagLimitppm),
FFTScaleShiftOverSampling(fftOverSampling)

{
  WinLength = (int)round(winSizeRelative * (double)dataLength);
  FFT_Size_SH = (int)pow(2.0, std::ceil(log2(dataLength)));

  //Pre-compute the Hanning Windows
  HannWindow = new double[WinLength];
  HannWindowCenter = new double[WinLength];
  for(int i = 0; i < WinLength; i++)
  {
    HannWindow[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( 2.0 * WinLength ) ));
    HannWindowCenter[i] = 0.5 * (1.0 - cos( 2.0 * M_PI * ((double)i + 1.0) / ( WinLength ) ));
  }
  
  FFT_Size_direct = (int)pow(2.0, std::ceil(log2(WinLength)));
  FFT_Size_inverse = (int)pow(2.0, std::ceil(log2(WinLength*FFTScaleShiftOverSampling)));
  fft_direct_in = fftw_alloc_real(FFT_Size_direct);
  fft_direct_out = fftw_alloc_real(FFT_Size_direct);
  fft_inverse_in = fftw_alloc_real(FFT_Size_inverse);
  fft_inverse_out = fftw_alloc_real(FFT_Size_inverse);
  fft_direct_shiftScale_in = fftw_alloc_complex(FFT_Size_SH);
  fft_direct_shiftScale_out = fftw_alloc_complex(FFT_Size_SH);
  
  fft_pdirect = fftw_plan_r2r_1d(FFT_Size_direct, fft_direct_in, fft_direct_out, FFTW_R2HC, FFTW_ESTIMATE);
  fft_pinvers = fftw_plan_r2r_1d(FFT_Size_inverse, fft_inverse_in, fft_inverse_out, FFTW_HC2R, FFTW_ESTIMATE);
  fft_pdshiftScale = fftw_plan_dft_1d(FFT_Size_SH, fft_direct_shiftScale_in, fft_direct_shiftScale_out, FFTW_FORWARD, FFTW_ESTIMATE);
  
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

  Rcpp::Rcout<<"\n===============================================\n";
  Rcpp::Rcout<<"DBG: refMassLowIndex = "<<refMassLowIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassMidIndex = "<<refMassMidIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassHighIndex = "<<refMassHighIndex<<"\n";
  Rcpp::Rcout<<"===============================================\n";
  
  refMassLowIndex = refMassLowIndex < 0 ? 0 : refMassLowIndex;
  refMassLowIndex = refMassLowIndex > (dataLength - 1) ? (dataLength - 1) : refMassLowIndex;
  
  refMassMidIndex = refMassMidIndex < 0 ? 0 : refMassMidIndex;
  refMassMidIndex = refMassMidIndex > (dataLength - 1) ? (dataLength - 1) : refMassMidIndex;
  
  refMassHighIndex = refMassHighIndex < 0 ? 0 : refMassHighIndex;
  refMassHighIndex = refMassHighIndex > (dataLength - 1) ? (dataLength - 1) : refMassHighIndex;
  
  Rcpp::Rcout<<"\n===============================================\n";
  Rcpp::Rcout<<"DBG: refMassLowIndex = "<<refMassLowIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassMidIndex = "<<refMassMidIndex<<"\n";
  Rcpp::Rcout<<"DBG: refMassHighIndex = "<<refMassHighIndex<<"\n";
  Rcpp::Rcout<<"===============================================\n";
  
  //The mutext shared for all threads
  fftwMtx = sharedMutex;
}

LabelFreeAlign::~LabelFreeAlign()
{
  fftw_destroy_plan(fft_pdirect);
  fftw_destroy_plan(fft_pinvers);
  fftw_destroy_plan(fft_pdshiftScale);
  fftw_free(fft_direct_in); 
  fftw_free(fft_direct_out);
  fftw_free(fft_inverse_in); 
  fftw_free(fft_inverse_out);
  fftw_free(fft_direct_shiftScale_in);
  fftw_free(fft_direct_shiftScale_out);
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
  
  
  //TODO simple test to check is the issue is realted with interpolation...
  /*for(int i = 0; i < FFT_Size_inverse; i++) data_ptr[i]=0.0;
  memcpy(data_ptr, data_ref, sizeof(double)*dataLength);
  return;*/
  
  CopyData2Window(data_ref, fft_direct_in, spectrumPart);

  //TODO simple test to check is the issue is realted with interpolation...
  /*for(int i = 0; i < FFT_Size_inverse; i++) data_ptr[i]=0.0;
  memcpy(data_ptr, fft_direct_in, sizeof(double)*FFT_Size_direct);
  return;*/
  
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
  
  //TODO simple test to check is the issue is realted with interpolation...
  /* for(int i = 0; i < FFT_Size_inverse; i++) data_ptr[i]=0.0;
   memcpy(data_ptr, fft_direct_out, sizeof(double)*FFT_Size_direct);
   return; */
  
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

LabelFreeAlign::TLags LabelFreeAlign::AlignSpectrum(double *interpolatedIntensitySpectrum, double *imzMLmass, int imzMLmassLength)
{
  //Hanning Windowing
  double *topWin_data = new double[FFT_Size_direct];
  double *midWin_data = new double[FFT_Size_direct];
  double *botWin_data = new double[FFT_Size_direct];
  TLags firstLag;
  
  for(int i = 0; i < AlignIterations; i++) //TODO with multiple iterations the processed mode will not work!
  {
    CopyData2Window(interpolatedIntensitySpectrum, topWin_data, TOP_SPECTRUM); //TODO warning! with multiple iterations the windowed spectrum could be trimmed multiple times!
    CopyData2Window(interpolatedIntensitySpectrum, midWin_data, CENTER_SPECTRUM);
    CopyData2Window(interpolatedIntensitySpectrum, botWin_data, BOTTOM_SPECTRUM);
    
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
    
    Rcpp::Rcout<<"\n===============================================\n";
    Rcpp::Rcout<<"DBG: lags.lagLow = "<<lags.lagLow<<"\n";
    Rcpp::Rcout<<"DBG: lags.lagMid = "<<lags.lagMid<<"\n";
    Rcpp::Rcout<<"DBG: lags.lagHigh = "<<lags.lagHigh<<"\n";
    Rcpp::Rcout<<"===============================================\n";
    
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
    
    if(i == 0)
    {
      firstLag = lags;
    }
  
    //Spectra warping constants
    double K1, K2, Sh1, Sh2;
    double Rl = (double) refMassLowIndex;
    double Rm = (double) refMassMidIndex;
    double Rh = (double) refMassHighIndex;
    
    if(bBilinear)
    {
      if(imzMLmassLength > 0)
      {
        //Align data in processed mode
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
        for(int i = 0; i < imzMLmassLength; i++)
        {
          if(imzMLmass[i] < commonMassAxis[refMassMidIndex])
          {
            imzMLmass[i] = imzMLmass[i]*K1 + Sh1;
          }
          else
          {
            imzMLmass[i] = imzMLmass[i]*K2 + Sh2;
          }
        }
      }
      else
      {
        //Align data in continuous mode
        K1 = (Rm + lags.lagMid - Rl - lags.lagLow)/(Rm - Rl); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
        Sh1 = (Rm* lags.lagLow - Rl*lags.lagMid)/(Rm - Rl); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
        
        K2 = (Rh + lags.lagHigh - Rm - lags.lagMid)/(Rh - Rm); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
        Sh2 = (Rh* lags.lagMid - Rm*lags.lagHigh)/(Rh - Rm); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
        
        //Apply the scaling and shift to original data, after that the pointer to data contains the aligned sptectrum, so copy data before
        double *data2 = new double[dataLength];
        memcpy(data2, interpolatedIntensitySpectrum, sizeof(double)*dataLength);
        
        FourierLinerScaleShift(interpolatedIntensitySpectrum, K1, Sh1); //Now data will contain the aligned left part of spectrum
        FourierLinerScaleShift(data2, K2, Sh2);
        
        //Merge the two vectors in a single one using the Rm as a center...
        const int iVectorCenterOffset =  (int)round(Rm); //0.5 is the central reference
        memcpy(interpolatedIntensitySpectrum+iVectorCenterOffset, data2+iVectorCenterOffset, sizeof(double)*(dataLength - iVectorCenterOffset));
        delete[] data2;
      }
    }
    else
    {
      if(imzMLmassLength > 0)
      {
        //Align data in processed mode
        
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
        for(int i = 0; i < imzMLmassLength; i++)
        {
          imzMLmass[i] = imzMLmass[i]*K1 + Sh1;
        }
      }
      else
      {
        //Align data in continuous mode
        K1 = (Rh + lags.lagHigh - Rl - lags.lagLow)/(Rh - Rl); //New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
        Sh1 = (Rh* lags.lagLow - Rl*lags.lagHigh)/(Rh - Rl); //If scaling is performed before shift. New implementation supporting Other refs diferents than 0 and N. Extremes are Rl and Rh
        
#ifdef EXTRA_DEBUG_INFO
        Rcpp::Rcout << "===============================================\n";
        Rcpp::Rcout << "FFT Shift & Scaling constants:\n";
        Rcpp::Rcout << "K1 = " << K1 << "\n";
        Rcpp::Rcout << "Sh1 = " << Sh1 << "\n";
        Rcpp::Rcout << "===============================================\n\n";
#endif
        
        //Apply the scaling and shift to original data, after that the pointer to data contains the aligned sptectrum
        FourierLinerScaleShift(interpolatedIntensitySpectrum, K1, Sh1);
      }
    }
  }
  
  delete[] topWin_data;
  delete[] botWin_data;
  delete[] midWin_data;
    
  return firstLag;
}

void LabelFreeAlign::FourierLinerScaleShift(double *data, double scaling, double shift)
{
  const int NewFFT_Size = FFT_Size_SH + (int)round( (scaling - 1.0) * (double)FFT_Size_SH);
  const int NewFFT_SizeOversampled = NewFFT_Size * FFTScaleShiftOverSampling;
  fftwMtx->lock();
  fftw_complex *fft_odd_in = fftw_alloc_complex(NewFFT_SizeOversampled);
  fftw_complex *fft_odd_out = fftw_alloc_complex(NewFFT_SizeOversampled);
  fftw_plan fft_pOddInvers = fftw_plan_dft_1d(NewFFT_SizeOversampled, fft_odd_in, fft_odd_out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftwMtx->unlock();
  
  //Copy data as complex number and add padding zeros and calculate spectrum Max
  double maxInit = 0.0;
  for(int i = 0; i < FFT_Size_SH; i++)
  {
    if(i < dataLength)
    {
      fft_direct_shiftScale_in[i][0] = data[i];
    }
    else
    {
      fft_direct_shiftScale_in[i][0] = 0.0;
    }
    fft_direct_shiftScale_in[i][1] = 0.0; //The imaginary part is always zero because input es real
    
    maxInit = fft_direct_shiftScale_in[i][0] > maxInit ? fft_direct_shiftScale_in[i][0] : maxInit;
  }
  fftw_execute(fft_pdshiftScale);
  
  //Copy data to fft_odd_in and apply scaling and shift
  double arg, Re, Im, auxRe, auxIm;
  for( int i = 0; i < NewFFT_SizeOversampled; i++ )
  {
    //Scaling in FFT domain
    if( i < FFT_Size_SH/2 && i < NewFFT_SizeOversampled/2)
    {
      fft_odd_in[i][0] = fft_direct_shiftScale_out[i][0]; //Copy real part
      fft_odd_in[i][1] = fft_direct_shiftScale_out[i][1]; //Copy imaginary part
    }
    else if( i >= (NewFFT_SizeOversampled - FFT_Size_SH/2))
    {
      fft_odd_in[i][0] = fft_direct_shiftScale_out[i - NewFFT_SizeOversampled + FFT_Size_SH][0]; //Copy real part
      fft_odd_in[i][1] = fft_direct_shiftScale_out[i - NewFFT_SizeOversampled + FFT_Size_SH][1]; //Copy imaginary part
    }
    else
    {
      //Zero padding
      fft_odd_in[i][0] = 0.0;
      fft_odd_in[i][1] = 0.0;
    }
    
    //Shift in FFT domain
    arg = -2.0*M_PI*(1.0/((double)NewFFT_SizeOversampled))*((double)i + 1.0)*round( ((double)FFTScaleShiftOverSampling) * shift);
    Re = cos(arg);
    Im = sin(arg);
    auxRe = fft_odd_in[i][0] * Re - fft_odd_in[i][1]*Im;
    auxIm = fft_odd_in[i][0] * Im + fft_odd_in[i][1]*Re;
    fft_odd_in[i][0] = auxRe;
    fft_odd_in[i][1] = auxIm;
  }
  fftw_execute(fft_pOddInvers);

  //Copy data to output vector downsampling it by averaging values
  double maxEnd = 0.0;
  int lastiUp = -FFTScaleShiftOverSampling/2;
  int iup = 0;
  double *arg_curr = new double[dataLength];
  double arg_aux;
  for( int i = 0; i < dataLength; i++)
  {
    data[i] = 0.0;
    arg_curr[i] = 0.0;
    while( iup <  ( lastiUp + FFTScaleShiftOverSampling ))
    {
      arg_aux = fabs(atan2(fft_odd_out[iup][1],fft_odd_out[iup][0])); //Abs of phase
      arg_curr[i] = arg_aux >  arg_curr[i] ? arg_aux : arg_curr[i]; //Keep the more intense
      data[i] += sqrt(fft_odd_out[iup][0]*fft_odd_out[iup][0] + fft_odd_out[iup][1]*fft_odd_out[iup][1]); //Modul
      iup++;
    }

    //Prepare lastiUp for next samples
    lastiUp = iup;
    
    //It is not necessary to divide data[i] by FFTScaleShiftOverSampling since I'll compensate for the max.
    maxEnd = data[i] > maxEnd ? data[i] : maxEnd;
  }
  
  //Compensate gain
  double argAverage3 = 0.0;
  int iWinC;
  if( maxEnd > 0.0 && maxInit > 0.0 )
  {
    for(int i = 0; i < dataLength; i++)
    {
      iWinC = i == 0 ? 1 : i;
      iWinC = i == (dataLength - 1)  ? (dataLength - 2) : i;
      argAverage3 =  arg_curr[iWinC - 1] +  arg_curr[iWinC] +  arg_curr[iWinC + 1];
      argAverage3 /= 3.0; 
      if(argAverage3 < M_PI/3.0 )
      {
        data[i] /= maxEnd/maxInit;
      }
      else
      {
        data[i] = 0.0;
      }
    }
  }
  
  fftwMtx->lock();
  fftw_destroy_plan(fft_pOddInvers);
  fftw_free(fft_odd_in);
  fftw_free(fft_odd_out);
  fftwMtx->unlock();
  delete[] arg_curr;
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
  
  Rcpp::Rcout<<"\n===============================================\n";
  Rcpp::Rcout<<"DBG: lag RAW = "<<lag<<"\n";
  Rcpp::Rcout<<"DBG: FFT_Size_direct = "<<FFT_Size_direct<<"\n";
  Rcpp::Rcout<<"DBG: FFT_Size_inverse = "<<FFT_Size_inverse<<"\n";
  Rcpp::Rcout<<"===============================================\n";
  
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
                                        NumericVector massProcessedMode, bool bilinear = false, 
                                        double lagRefLow = 0.1, double lagRefMid = 0.5, double lagRefHigh = 0.9,
                                        int iterations = 1, double lagLimitppm = 200, int fftOverSampling = 10, double winSizeRelative = 0.6 )
{
  std::mutex mtx;
  NumericVector y(spectrumInterpolated.length());
  memcpy(y.begin(), spectrumInterpolated.begin(), sizeof(double)*spectrumInterpolated.length());
  
  LabelFreeAlign alngObj(mass.begin(), ref.begin(), ref.length(), 
                         bilinear, &mtx, iterations, 
                         lagRefLow, lagRefMid, lagRefHigh, lagLimitppm, fftOverSampling, winSizeRelative);
  
  Rcpp::NumericVector refLow = alngObj.getRefLowFFT();
  Rcpp::NumericVector refMid = alngObj.getRefCenterFFT();
  Rcpp::NumericVector refHigh = alngObj.getRefHighFFT();
  
  LabelFreeAlign::TLags lags = alngObj.AlignSpectrum(y.begin(), massProcessedMode.begin(), massProcessedMode.length());

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
#elif
  return List::create( Named("InterpolatedSpectrum") = y, 
                       Named("ProcessedModeMass") = massProcessedMode, 
                       Named("RefLowFFT") = refLow,
                       Named("RefMidFFT") = refMid,
                       Named("RefHighFFT") = refHigh
  );
#endif
}


