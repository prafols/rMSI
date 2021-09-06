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
#include "peakbinning.h"
#include "progressbar.h"
#include <limits>       // std::numeric_limits
#include <stdexcept>
using namespace Rcpp;

PeakBinning::PeakBinning(Rcpp::Reference preProcessingParams, int numberOfThreads) //TODO currently the numberOfThreads is not used since the implementation is not multithreaded.. think about it
{
  //Start with total number of pixels set to zero, it will be increased when data is appended
  totalNumOfPixels = 0;
  
  //Get the parameters
  Rcpp::Reference binningParams = preProcessingParams.field("peakbinning");
  tolerance = binningParams.field("tolerance");
  tolerance_in_ppm = binningParams.field("tolerance_in_ppm");
  binFilter = binningParams.field("binFilter"); 
}

PeakBinning::~PeakBinning()
{
  for( int i = 0; i < imzMLReaders.size(); i++)
  {
    delete imzMLReaders[i];
  }
}

void PeakBinning::appedImageData(Rcpp::List imzMLDescriptor)
{
  //Set the imzML reader
  DataFrame imzMLrun = as<DataFrame>(imzMLDescriptor["run_data"]);
  std::string sFilePath = as<std::string>(imzMLDescriptor["path"]);
  std::string sFnameImzML = as<std::string>(imzMLDescriptor["file"]);
  std::string sFnameImzMLInput = sFilePath + "/" + sFnameImzML + ".ibd";
  bool bPeakListrMSIformat = as<bool>(imzMLDescriptor["rMSIpeakList"]);
  
  if(!bPeakListrMSIformat && !tolerance_in_ppm)
  {
    throw std::runtime_error("ERROR: a peak list is not in rMSI imzML format, please set the bining tolerance in ppm to use non-rMSI peak lists.\n");
  }
  
  //Append the imzML Reader
  imzMLReaders.push_back(new ImzMLBinRead(sFnameImzMLInput.c_str(), 
                                          imzMLrun.nrows(), 
                                          as<String>(imzMLDescriptor["mz_dataType"]),
                                          as<String>(imzMLDescriptor["int_dataType"]) ,
                                          false, //A peak list is allways processed mode
                                          false, //Do not call the file open() on constructor
                                          bPeakListrMSIformat
                                          )); 
  
  NumericVector imzML_mzLength = imzMLrun["mzLength"];
  NumericVector imzML_mzOffsets = imzMLrun["mzOffset"];
  NumericVector imzML_intLength = imzMLrun["intLength"];
  NumericVector imzML_intOffsets = imzMLrun["intOffset"];
  imzMLReaders.back()->set_mzLength(&imzML_mzLength);  
  imzMLReaders.back()->set_mzOffset(&imzML_mzOffsets);
  imzMLReaders.back()->set_intLength(&imzML_intLength);
  imzMLReaders.back()->set_intOffset(&imzML_intOffsets);
  
  totalNumOfPixels += imzMLrun.nrows();
}

List PeakBinning::BinPeaks()
{
  Rcout<<"Binning peaks...\n";
  unsigned int iCount = 0; //Used for the progress bar
  PeakPicking::Peaks *mpeaks; //Pointer to the current peaklist
  
  //The mass name for each matrix column
  std::vector<double> binMass; 
  
  //Store the number of peaks added to each column to apply the bin filter at the end
  std::vector<unsigned int> columnsPeakCounters;
  
  //Peak matrices starting with zero columns
  std::vector<std::vector<TBin> > binMat; //A matrix containing binning values
  
  for(int iimg = 0; iimg < imzMLReaders.size(); iimg++)
  {
    imzMLReaders[iimg]->open();
    for(int ipixel = 0; ipixel < imzMLReaders[iimg]->get_number_of_pixels(); ipixel++)
    {
      progressBar(iCount, totalNumOfPixels, "=", " ");
      iCount++;
      mpeaks = imzMLReaders[iimg]->ReadPeakList(ipixel);
      for(int imass = 0; imass < mpeaks->mass.size(); imass++)
      {
        //Search the closest mass bin
        double minMassDistance = std::numeric_limits<double>::max();
        int minDistanceIndex = -1;
        double currentMassDistance;
        for(int icol = 0; icol < binMass.size(); icol++)
        {
          currentMassDistance = fabs(mpeaks->mass[imass] - binMass[icol]);
          if(currentMassDistance < minMassDistance)
          {
            minMassDistance = currentMassDistance;
            minDistanceIndex = icol;
          }
        }
        
        //Select the kind of binning tolerance
        double compTolerance;
        if(tolerance_in_ppm)
        {
          minMassDistance = 1e6*(minMassDistance/mpeaks->mass[imass]); //Compute distance in ppm
          compTolerance = tolerance;
        }
        else
        {
          compTolerance = tolerance * mpeaks->binSize[imass];
        } 
        
        /**** DEBUG Pints
        Rcout << "ipixel = " << ipixel << " of " << totalNumOfPixels << "\n"; 
        Rcout << "mpeaks->mass[ " << imass << " ] = " << mpeaks->mass[imass] << "\n"; 
        Rcout << "minMassDistance = " << minMassDistance << "\n";
        Rcout << "minDistanceIndex = " << minDistanceIndex << "\n";
        Rcout << "compTolerance = " << compTolerance << "\n";
        Rcout << "binMat ncols = " << binMat.size() << "\n";
        if(binMat.size() > 20000)
        {
          stop("ABORTING!!!! too much columns!!!\n");
        }
        */
        
        if( (minMassDistance < compTolerance) && (minDistanceIndex >= 0) )
        {
          //The peak must be binned with the minDistanceIndex column of the peak matrix
          columnsPeakCounters[minDistanceIndex]++;
          binMass[minDistanceIndex] = 0.5*(binMass[minDistanceIndex] + mpeaks->mass[imass]); //Recompute the peak matrix mass by simply averaging it
          binMat[minDistanceIndex][ipixel].intensity = mpeaks->intensity[imass] > binMat[minDistanceIndex][ipixel].intensity ? mpeaks->intensity[imass] : binMat[minDistanceIndex][ipixel].intensity;
          if(imzMLReaders[iimg]->get_rMSIPeakListFormat())
          {
            binMat[minDistanceIndex][ipixel].SNR = mpeaks->SNR[imass] > binMat[minDistanceIndex][ipixel].SNR ? mpeaks->SNR[imass] : binMat[minDistanceIndex][ipixel].SNR;
            binMat[minDistanceIndex][ipixel].area = mpeaks->area[imass] > binMat[minDistanceIndex][ipixel].area ? binMat[minDistanceIndex][ipixel].area : binMat[minDistanceIndex][ipixel].area;
          }
        }
        else
        {
          //A new column must be added to the peak matrix
          binMass.push_back(mpeaks->mass[imass]); //Append element to the mass vector (names of bin Matrix)
          columnsPeakCounters.push_back(1);
          binMat.resize(binMat.size() + 1); //Append new column
          binMat[binMat.size() - 1].resize(totalNumOfPixels); //Extenend all new column elements
          binMat[binMat.size() - 1][ipixel].intensity = mpeaks->intensity[imass];
          if(imzMLReaders[iimg]->get_rMSIPeakListFormat())
          {
            binMat[binMat.size() - 1][ipixel].SNR = mpeaks->SNR[imass]; 
            binMat[binMat.size() - 1][ipixel].area = mpeaks->area[imass]; 
          }
        }
      }
      
      delete mpeaks;
    }
    imzMLReaders[iimg]->close();
  }
  
  //Apply binFilter
  std::vector<bool> keepColumns(columnsPeakCounters.size());
  unsigned int number_of_columns_to_remove = 0;
  for( int i=0; i < keepColumns.size(); i++)
  {
    keepColumns[i] = (double)columnsPeakCounters[i] > binFilter*(double)totalNumOfPixels;
    if(!keepColumns[i])
    {
      number_of_columns_to_remove++;
    }
  }
  if(number_of_columns_to_remove >= binMass.size())
  {
    throw std::runtime_error("ERROR: all peaks were removed by the bin filter. Consider decresin the bin filter or SNR.\n");
  }

  Rcout<<"Bining complete with a total number of "<<binMass.size() - number_of_columns_to_remove<<" bins\n"; 
  
  //Copy data to R matrices
  Rcout<<"Coping data to R object...\n";
  NumericMatrix binMatIntensity(totalNumOfPixels, binMass.size() - number_of_columns_to_remove);
  NumericMatrix binMatSNR(totalNumOfPixels, binMass.size() - number_of_columns_to_remove);
  NumericMatrix binMatArea(totalNumOfPixels, binMass.size()- number_of_columns_to_remove);
  
  //Sort columns by mass
  Rcout<<"Sorting columns by mass...\n";
  NumericVector massCopy(binMass.size());
  NumericVector massSorted(binMass.size() - number_of_columns_to_remove);
  memcpy(massCopy.begin(), binMass.data(), sizeof(double)*binMass.size());
  int sortedInds[binMass.size()];
  double minVal;
  for(int i = 0; i < binMass.size(); i++)
  {
    minVal = std::numeric_limits<double>::max();
    for( int j = 0; j < binMass.size(); j++ )
    {
      if( massCopy[j] < minVal && massCopy[j] > 0 )
      {
        minVal = massCopy[j];
        sortedInds[i] = j;
      }  
    }
    massCopy[ sortedInds[i]  ] = -1; //Mark as sorted
  }
  
  //Copy the mass axis sorting it
  unsigned int icopy = 0;
  for(int i = 0; i < binMass.size(); i++)
  {
    if(keepColumns[i])
    {
      massSorted[icopy] = binMass[sortedInds[i]];
      icopy++;
    }
  }
  
  //Copy the matrix sorting it
  for( int ir = 0;  ir < totalNumOfPixels; ir++)
  {
    icopy = 0;
    for(int ic = 0; ic < binMass.size(); ic++)
    {
      if(keepColumns[ic])
      {
        binMatIntensity(ir, icopy ) = binMat[sortedInds[ic]][ir].intensity;
        binMatSNR      (ir, icopy ) = binMat[sortedInds[ic]][ir].SNR;
        binMatArea     (ir, icopy ) = binMat[sortedInds[ic]][ir].area;
        icopy++;
      }
    }
  }
  
  return List::create( Named("mass") = massSorted, Named("intensity") = binMatIntensity, Named("SNR") = binMatSNR, Named("area") = binMatArea );
}

// [[Rcpp::export]]
List CRunPeakBinning(Rcpp::List imzMLDescriptor, Rcpp::Reference preProcessingParams, int numOfThreads)
{
  List out;
  try
  {
    PeakBinning myPeakBinning(preProcessingParams, numOfThreads);
    
    for(int i = 0; i < imzMLDescriptor.length(); i++)
    {
      myPeakBinning.appedImageData(imzMLDescriptor[i]);
    }
    
    out = myPeakBinning.BinPeaks(); 
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return out;
}

