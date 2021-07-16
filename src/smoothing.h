/*************************************************************************
 * rMSIproc - R package for MSI data processing
 * Copyright (C) 2014 Pere Rafols Soler
 *
 * This file is part of rMSIproc.
 *
 * rMSIproc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * rMSIproc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rMSIproc.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/
#ifndef SMOOTHING_H
  #define SMOOTHING_H

#include <Rcpp.h>
  
class Smoothing
{
  public:
    Smoothing(int kernelSize = 5);
    ~Smoothing();
    Rcpp::NumericVector smoothSavitzkyGolay(Rcpp::NumericVector x);
    void smoothSavitzkyGolay(double *x, int length);
    
    //Smooth a given spectrum and interpolate it to imzML processed mode. 
    // Arguments:
    // - commonMassAxis: a pointer to the common mass axis used in the interpolated spectrum.
    // - intensityDataInterpolated: a pointer to the intensity spectrum. If data is in processed mode this argument corresponds to the interpolated spectrum.
    // - length: length of the interpolated intensity vector.
    // -massData: a pointer to the mass axis independently if data is in processed or continuous mode (as is in the imzML file)
    // -intensityData: a pointer to the mass axis independently if data is in processed or continuous mode (as is in the imzML file)
    // - N: number of mass channels in the spectrum massData. Set to zero for continous data mode.
    void smoothSavitzkyGolay(double *commonMassAxis, double *intensityDataInterpolated, int length, double *massData, double *intensityData, int N); 
    
  private:
    std::vector<double> sgC; //This is the SavitzkyGolay kernel
};


#endif