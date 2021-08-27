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
  #include "peakpicking.h"
using namespace Rcpp;

#define AREA_WINDOW_SIDE_WIDTH 3

PeakPicking::PeakPicking(int WinSize, double *massAxis, int numOfDataPoints, int UpSampling )
{
  FFT_Size = (int)pow(2.0, std::ceil(log2(WinSize)));
  FFT_Size = FFT_Size < 16 ? 16 : FFT_Size; //Minimum allowed windows size is 16 points.
  FFTInter_Size = (int)pow(2.0, std::ceil(log2(UpSampling*FFT_Size))); //FFT interpolation buffer
  
  dataLength = numOfDataPoints;
  mass = new double[dataLength];
  memcpy(mass, massAxis, sizeof(double) * dataLength);
  
  //Compute Hanning and Area windows
  HanningWin = new double[FFT_Size];
  AreaWin = new double[FFT_Size];
  
  const int leftHannLength = 2*((FFT_Size/2) - 1) + 1;
  const int rightHannLength = 2*(FFT_Size - (1+FFT_Size/2)) - 1;
  
  for( int i = 0; i < FFT_Size; i++)
  {
    if(i < (FFT_Size/2) )
    {
      //Left part hanning
      HanningWin[i] = 0.5*(1.0-cos((2.0*M_PI*(double)(i))/(double)(leftHannLength - 1))); 
    }
    else if( i > (FFT_Size/2) )
    {
      //Right part hanning
      HanningWin[i] = 0.5*(1.0-cos((2.0*M_PI*(double)(i - 3))/(double)(rightHannLength - 1)));
    }
    else
    {
      //The tree sample at center (peak will be always at FFT_Size/2)
      HanningWin[i] = 1.0;
    }
    
    if(i < AREA_WINDOW_SIDE_WIDTH)
    {
      AreaWin[i] = 0.5*(1.0-cos((2.0*M_PI*(double)(i))/(double)(2*AREA_WINDOW_SIDE_WIDTH - 1)));
    }
    else if( i >= (FFT_Size - AREA_WINDOW_SIDE_WIDTH) )
    {
      AreaWin[i] = 0.5*(1.0-cos((2.0*M_PI*(double)(i - (FFT_Size - 2*AREA_WINDOW_SIDE_WIDTH )))/(double)(2*AREA_WINDOW_SIDE_WIDTH - 1)));
    }
    else
    {
      //Area central part
      AreaWin[i] = 1.0;
    }
  }
  
  //Prepare FFT objects
  fft_in1 = fftw_alloc_real(FFT_Size);
  fft_out1 = fftw_alloc_real(FFT_Size);
  fft_in2 = fftw_alloc_real(FFTInter_Size);
  fft_out2 = fftw_alloc_real(FFTInter_Size);
  fft_pdirect = fftw_plan_r2r_1d(FFT_Size, fft_in1, fft_out1, FFTW_R2HC, FFTW_ESTIMATE);
  fft_pinvers = fftw_plan_r2r_1d(FFTInter_Size, fft_in2, fft_out2, FFTW_HC2R, FFTW_ESTIMATE);
  
  //Prepare NoiseEstimation oject
  neObj = new NoiseEstimation(dataLength);
}

PeakPicking::~PeakPicking()
{
  delete[] mass;
  delete[] HanningWin;
  delete[] AreaWin;
  fftw_destroy_plan(fft_pdirect);
  fftw_destroy_plan(fft_pinvers);
  fftw_free(fft_in1); 
  fftw_free(fft_out1);
  fftw_free(fft_in2); 
  fftw_free(fft_out2);
  delete neObj;
}

PeakPicking::Peaks *PeakPicking::peakPicking(double *spectrum, double SNR )
{
  //Calculate noise
  double *noise = new double[dataLength];
  memcpy(noise, spectrum, sizeof(double)*dataLength);
  neObj->NoiseEstimationFFTExpWin(noise, dataLength, FFT_Size*2);
  
  //Detect peaks
  PeakPicking::Peaks *pks = detectPeaks(spectrum, noise, SNR);
  delete[] noise;
  return pks;
}

//Detect all local maximums and Filter peaks using SNR min value in a sliding window
PeakPicking::Peaks *PeakPicking::detectPeaks( double *spectrum, double *noise, double SNR )
{
  const int HalfWinSize = FFT_Size/2;
  double slope = 0.0;
  double slope_ant = 0.0;
  PeakPicking::Peaks *m_peaks = new PeakPicking::Peaks();
  for( int i=0; i < (dataLength - 1); i++)
  {
    slope = spectrum[i + 1] - spectrum[i]; //Compute 1st derivative
    //Look for a zero crossing at first derivate (negative sign of product) and negative 2nd derivate value (local maxim)
    if(slope*slope_ant <= 0.0 && (slope - slope_ant) < 0.0)
    {
      if(spectrum[i]/noise[i] >= SNR)
      {
        m_peaks->mass.push_back( predictPeakMass(spectrum, i)); //Compute peak accurately using FFT interpolation); 
        m_peaks->intensity.push_back(spectrum[i]);
        m_peaks->SNR.push_back(spectrum[i]/noise[i]); 
        m_peaks->area.push_back(predictPeakArea(spectrum, i)); //Normalized to non-FFT sapce
        m_peaks->binSize.push_back( fabs(mass[i + 1] - mass[i]) );
      }
    }
    slope_ant = slope;
  }
  return m_peaks;
}

int PeakPicking::interpolateFFT(double *spectrum, int iPeakMass, bool ApplyHanning)
{
  //Fill fft data vector taking care of extrems
  int idata = 0; //Indeix of data
  for( int i = 0; i < FFT_Size; i++)
  {
    idata = iPeakMass - FFT_Size/2 + i;
    fft_in1[i] = (idata >= 0 && idata < dataLength) ? spectrum[idata] : 0.0;
    if(ApplyHanning)
    {
      fft_in1[i] *= HanningWin[i]; 
    }
    else
    {
      //No Hanning here so apply area windows to avoid Gibbs
      fft_in1[i] *= AreaWin[i];
    }
  }
  
  //Compute FFT
  fftw_execute(fft_pdirect);
  
  //Zero padding in fft space
  int i1 = 0;
  for( int i = 0; i < FFTInter_Size; i++)
  {
    if( i <= (FFT_Size/2) || (FFTInter_Size - i) <= ((FFT_Size + 1)/2 - 1 )  )
    {
      fft_in2[i] = fft_out1[i1];
      i1++;
    }
    else
    {
      fft_in2[i] = 0.0;  
    }
  }
  
  //Compute Inverse FFT, fft_out2 contains the interpolated peak
  fftw_execute(fft_pinvers);
  
  //Peak is around FFT_Size/2 so just look there for it
  int imax = -1;
  double dmax = -1.0;
  for( int i = (int)std::floor((((double)FFTInter_Size)/((double)FFT_Size))*(0.5*(double)FFT_Size - 1.0)); 
       i <= (int)std::ceil((((double)FFTInter_Size)/((double)FFT_Size))*(0.5*(double)FFT_Size + 1.0)); 
       i++)
  {
    if( fft_out2[i] > dmax )
    {
      imax = i;
      dmax = fft_out2[i];
    }
  }
  
  return imax;
}

double PeakPicking::predictPeakMass( double *spectrum, int iPeakMass )
{
  double pMass; 
  int iPeak = interpolateFFT(spectrum, iPeakMass, true);
  
  //Compute the original mass indexing space
  double imass = (double)iPeakMass - (0.5*(double)FFT_Size) + ((double)iPeak) * ((double)FFT_Size)/((double)FFTInter_Size);

  //Convert peak position indexes to mass values
  double peakPosL = std::floor(imass);
  double peakPosR = std::ceil(imass);
  if(peakPosL == peakPosR)
  {
    //Peak exactly on mass channel, no decimal part
    pMass =  mass[(int)imass];
  }
  else
  {
    //Calculate peak position as a compositon of both neighbours mass channels
    pMass = (peakPosR - imass)*mass[(int)peakPosL] + (imass - peakPosL)*mass[(int)peakPosR];
  }

  return pMass;
}

double PeakPicking::predictPeakArea( double *spectrum, int iPeakMass )
{
  //Calculate integration range befor interpolation to avoid FFT Gibbs issues (left part).
  int integrationLimitLeft = FFT_Size/2;
  for( int i = (iPeakMass - 1); i >= (iPeakMass - FFT_Size/2); i-- )
  {
    if( i < 0)
    {
      //Break if out of spectrum
      break;
    }
    if(spectrum[i] >= spectrum[i + 1] ) 
    {
      //Break loop if slope increases too much at some point
      break; 
    }
    integrationLimitLeft--;
  }
  
  //Calculate integration range befor interpolation to avoid FFT Gibbs issues (rigth part).
  int integrationLimitRight = FFT_Size/2;
  for( int i = (iPeakMass + 1); i < (iPeakMass + FFT_Size/2); i++)
  {
    if( i >= dataLength)
    {
      //Break if out of spectrum
      break;
    }
    if(spectrum[i] >= spectrum[i - 1] ) 
    {
      //Break loop if slope increases too much at some point
      break;
    }
    integrationLimitRight++;
  }
  
  //Calculate mass step size linear approx on interpolated space
  int massIdL = iPeakMass - FFT_Size/2;
  int massIdR = iPeakMass + FFT_Size/2 - 1;
  massIdL = massIdL < 0 ? 0 : massIdL;
  massIdR = massIdR >= dataLength ? (dataLength - 1) : massIdR;
  double massStep = (mass[massIdR] - mass[massIdL])/((double)FFTInter_Size);
  
  //Transform integration range to interpolated indexes
  const double rangeConverter = ((double)FFTInter_Size)/((double)FFT_Size);
  integrationLimitLeft = (int)round(((double)integrationLimitLeft) * rangeConverter);
  integrationLimitRight = (int)round(((double)integrationLimitRight) * rangeConverter);
  integrationLimitLeft = integrationLimitLeft < 0 ? 0 : integrationLimitLeft;
  integrationLimitRight = integrationLimitRight > (FFTInter_Size - 1) ? (FFTInter_Size - 1) : integrationLimitRight;

  //Peak interpolation
  int iPeak = interpolateFFT(spectrum, iPeakMass, false);
  
  //Area integration
  double pArea = 0.0;
  for( int i = integrationLimitLeft; i <= integrationLimitRight; i++)
  {
    pArea+=fft_out2[i];
  }
 
  if( fft_out2[iPeak] == 0.0 )
  {
    pArea = 0.0; //Avoid zero division
  }
  else
  {
    pArea *= (spectrum[iPeakMass] / fft_out2[iPeak]) * massStep;
  }
  return pArea; //Return the area de-normalizing the FFT space
}

//Returns the internaly used Hanning windows (only for test purposes)
NumericVector PeakPicking::getHannWin()
{
  NumericVector hannR(FFT_Size);
  memcpy(hannR.begin(), HanningWin, sizeof(double)*FFT_Size);
  return hannR;
}

//Convert a pointer to a Peaks structure to an R list
List PeakPicking::PeakObj2List(PeakPicking::Peaks *pks)
{
  NumericVector mass(pks->mass.size());
  NumericVector intensity(pks->intensity.size());
  NumericVector SNR(pks->SNR.size());
  NumericVector area(pks->area.size());
  for(int i = 0; i < mass.length(); i++ )
  {
    mass(i) = pks->mass[i];
    intensity(i) = pks->intensity[i];
    SNR(i) = pks->SNR[i];
    area(i) = pks->area[i];
  }
  
  return List::create( Named("mass") = mass, Named("intensity") = intensity, Named("SNR") = SNR);
}

NumericVector PeakPicking::getInterpolatedPeak(NumericVector data, int iPeak, bool ApplyHanning)
{
  double *Cdata = new double[data.length()];
  memcpy(Cdata, data.begin(), sizeof(double)*data.length());
  interpolateFFT(Cdata, iPeak, ApplyHanning);
    
  NumericVector interpData(FFTInter_Size);
  memcpy(interpData.begin(), fft_out2, sizeof(double)*FFTInter_Size);
  delete[] Cdata;
  return(interpData);
}

NumericVector PeakPicking::getHanningWindow()
{
  NumericVector RHann(FFT_Size);
  memcpy(RHann.begin(), HanningWin, sizeof(double)*FFT_Size);
  return RHann;
}

NumericVector PeakPicking::getAreaWindow()
{
  NumericVector RArea(FFT_Size);
  memcpy(RArea.begin(), AreaWin, sizeof(double)*FFT_Size);
  return RArea;
}

////// Rcpp Exported methods //////////////////////////////////////////////////////////
//' DetectPeaks_C.
//' 
//' Detect peaks from a Rcpp::NumericVector object and returns data in a R matrix.
//' This method is only exported to be use by R function DetectPeaks which is an actual R function.
//' The returned peak positions follows C indexing style, this is starts with zero.
//' 
//' @param mass a NumericVector containing the mass axis of the spectrum.
//' @param intensity a NumericVector where peaks must be detected.
//' @param SNR Only peaks with an equal or higher SNR are retained.
//' @param WinSize The windows used to detect peaks and caculate noise.
//' @param UpSampling the oversampling used for acurate mass detection and area integration.
//' 
//' @return a NumerixMatrix of 5 rows corresponding to: mass, intensity of the peak, SNR, area and binSize.
//' 
// [[Rcpp::export]]
NumericMatrix DetectPeaks_C(NumericVector mass, NumericVector intensity, double SNR = 5, int WinSize = 20, int UpSampling = 10)
{
  if(mass.length() != intensity.length())
  {
    Rcpp::stop("Error in DetectPeaks_C() function: mass and intensity length must be equal.");
    return NumericMatrix(0,0);
  }
    
  double *massC = new double[mass.length()];
  double *spectrum = new double[mass.length()];
  
  //Copy R data to C arrays
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  memcpy(spectrum, intensity.begin(), sizeof(double)*intensity.length());
  
  PeakPicking ppObj(WinSize, massC, mass.length(), UpSampling);
  PeakPicking::Peaks *peaks = ppObj.peakPicking(spectrum, SNR);
  
  //Convert peaks to R matrix like object
  NumericMatrix mp(5, peaks->intensity.size());
  for(int i = 0; i < peaks->mass.size(); i++)
  {
    mp(0, i) = peaks->mass[i];
    mp(1, i) = peaks->intensity[i];
    mp(2, i) = peaks->SNR[i];
    mp(3, i) = peaks->area[i];
    mp(4, i) = peaks->binSize[i];
  }
  
  delete peaks;
  delete[] massC;
  delete[] spectrum;
  rownames(mp) =  CharacterVector::create("mass", "intensity", "SNR", "area", "binSize");
  return mp;
}
//' TestPeakInterpolation_C.
//' 
//' 
//' @param mass a NumericVector containing the mass axis of the spectrum.
//' @param intensity a NumericVector where peaks must be detected.
//' @param peakIndex the location of the peak to interpolate in the spectrum.  
//' @param WinSize The windows used to detect peaks and caculate noise.
//' @param UpSampling the oversampling used for acurate mass detection and area integration.
//' @param useHanning if hanning windowing must be used befor interpolation.
//' @param Iterations number of iterations to perform. This is just for testing interpolation efficiency
//' 
//' @return a NumerixVector with the FFT interpolated peak shape.
//' 
// [[Rcpp::export]]
NumericVector TestPeakInterpolation_C(NumericVector mass, NumericVector intensity, int peakIndex, int WinSize = 20, int UpSampling = 10, bool useHanning = false, int Iterations = 1)
{
  if(mass.length() != intensity.length())
  {
    Rcpp::stop("Error in DetectPeaks_C() function: mass and intensity length must be equal.");
    return NumericVector(0);
  }
  
  double *massC = new double[mass.length()];
  double *spectrum = new double[mass.length()];
  
  //Copy R data to C arrays
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  memcpy(spectrum, intensity.begin(), sizeof(double)*intensity.length());
  
  PeakPicking ppObj(WinSize, massC, mass.length(), UpSampling);
  
  NumericVector interpolated;
  for( int i = 0; i < Iterations; i++)
  {
    interpolated = ppObj.getInterpolatedPeak(intensity, peakIndex, useHanning);
  }
  
  delete[] massC;
  delete[] spectrum;
  
  return(interpolated);
}

//' TestHanningWindow.
//' 
//' Method to test the implementation of Hanning window in R session.
//' @param mass a NumericVector containing the mass axis of the spectrum.
//' @param WinSize The windows used to detect peaks and caculate noise.
//' @param UpSampling the oversampling used for acurate mass detection and area integration.
//' 
//' @return a NumericVector containing the Hanning Window.
//' 
// [[Rcpp::export]]
NumericVector TestHanningWindow(NumericVector mass, int WinSize = 20, int UpSampling = 10)
{
  double *massC = new double[mass.length()];
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  PeakPicking ppObj(WinSize, massC, mass.length(), UpSampling);
  NumericVector hannWin = ppObj.getHanningWindow();
  delete[] massC;
  return hannWin;
}

//' TestAreaWindow.
//' 
//' Method to test the implementation of Area window in R session.
//' @param mass a NumericVector containing the mass axis of the spectrum.
//' @param WinSize The windows used to detect peaks and caculate noise.
//' @param UpSampling the oversampling used for acurate mass detection and area integration.
//' 
//' @return a NumericVector containing the Area Window.
//' 
// [[Rcpp::export]]
NumericVector TestAreaWindow(NumericVector mass, int WinSize = 20, int UpSampling = 10)
{
  double *massC = new double[mass.length()];
  memcpy(massC, mass.begin(), sizeof(double)*mass.length());
  PeakPicking ppObj(WinSize, massC, mass.length(), UpSampling);
  NumericVector areaWin = ppObj.getAreaWindow();
  delete[] massC;
  return areaWin;
}
