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
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h> 
#include <stdexcept>
#include "rmsicdatacubeio.h"
using namespace Rcpp;


CrMSIDataCubeIO::CrMSIDataCubeIO(Rcpp::NumericVector massAxis, double cubeMemoryLimitMB, bool storeDataInImzML, Rcpp::String imzMLOutputPath)
  :mass(massAxis), storeData(storeDataInImzML), dataOutputPath(imzMLOutputPath.get_cstring())
{
  cubeMaxNumRows = std::ceil((1024*1024*cubeMemoryLimitMB)/(double)(8*mass.length()));
}

CrMSIDataCubeIO::~CrMSIDataCubeIO()
{
  for( int i = 0; i < imzMLReaders.size(); i++)
  {
    delete imzMLReaders[i];
  }
  
  for( int i = 0; i < imzMLWriters.size(); i++)
  {
    delete imzMLWriters[i];
  }
}

void CrMSIDataCubeIO::appedImageData(Rcpp::List rMSIobj, std::string outputImzMLuuid, std::string outputImzMLfname)
{
  if(outputImzMLuuid.empty() && storeData)
  {
    throw std::runtime_error("Error: A valid UUID must be provided for the stored imzML file.\n");
  }
  
  if(outputImzMLfname.empty() && storeData)
  {
    throw std::runtime_error("Error: A filename suffix must be provided for the stored imzML file.\n");
  }
  
  if(outputImzMLuuid.length() != 32 && storeData)
  {
    throw std::runtime_error("Error: The UUID must contain 16 bytes.\n");
  }
  
  //Set the imzML reader
  List data = rMSIobj["data"];
  List imzML = data["imzML"];
  DataFrame imzMLrun = as<DataFrame>(imzML["run"]);
  std::string sFilePath = as<std::string>(data["path"]);
  std::string sFnameImzML = as<std::string>(imzML["file"]);
  std::string sFnameImzMLInput = sFilePath + "/" + sFnameImzML + ".ibd";
  
  //Rcpp::Rcout << "CrMSIDataCubeIO::appedImageData()--> sFnameImzML = "<< sFnameImzML << std::endl; //DEBUG line

  //Append the imzML Reader
  imzMLReaders.push_back(new ImzMLBinRead(sFnameImzMLInput.c_str(), 
                                     imzMLrun.nrows(), 
                                     as<String>(imzML["mz_dataType"]),
                                     as<String>(imzML["int_dataType"]) ,
                                     as<bool>(imzML["continuous_mode"]), 
                                     false)); //Do not call the file open() on constructor
  
  NumericVector imzML_mzLength = imzMLrun["mzLength"];
  NumericVector imzML_mzOffsets = imzMLrun["mzOffset"];
  NumericVector imzML_intLength = imzMLrun["intLength"];
  NumericVector imzML_intOffsets = imzMLrun["intOffset"];
  imzMLReaders.back()->set_mzLength(&imzML_mzLength);  
  imzMLReaders.back()->set_mzOffset(&imzML_mzOffsets);
  imzMLReaders.back()->set_intLength(&imzML_intLength);
  imzMLReaders.back()->set_intOffset(&imzML_intOffsets);
  
  
  //Create the imzMLwriter corresponding to the current imzMLreader
  if(storeData)
  {
    std::string sFnameImzMLOutput= dataOutputPath + "/" + outputImzMLfname + ".ibd"; 
    imzMLWriters.push_back(new ImzMLBinWrite(sFnameImzMLOutput.c_str(),
                                            imzMLrun.nrows(), 
                                            as<String>(imzML["mz_dataType"]),
                                            as<String>(imzML["int_dataType"]) ,
                                            as<bool>(imzML["continuous_mode"]),
                                            true, //Data cubes will be sotred using the sequential mode
                                            false)); //Dont call the file open() on constructor to allow directely writing the uuid, then I'll close it!
    imzMLWriters.back()->open(true); //Open and truncate the file
    imzMLWriters.back()->writeUUID(outputImzMLuuid);
    imzMLWriters.back()->close(); 
    
    acumulatedSpectrum.push_back(Rcpp::NumericVector(mass.length()));
    baseSpectrum.push_back(Rcpp::NumericVector(mass.length()));
  }
  
  //Initialize the cube description if this is the first call to appedImageData()
  if(dataCubesDesc.size() == 0)
  {
    dataCubesDesc.push_back(DataCubeDescription());
  }
  
  for( unsigned int i = 0; i < imzMLrun.nrows(); i++) //For each pixel in imzML file
  {
    if(dataCubesDesc.back().size() == cubeMaxNumRows)
    {
      //No more space left in the cube, so create a new cube
      dataCubesDesc.push_back(DataCubeDescription());
    }
    
    dataCubesDesc.back().push_back(PixelDescription());
    dataCubesDesc.back().back().pixel_ID = i;
    dataCubesDesc.back().back().imzML_ID = imzMLReaders.size()-1;
  }
}

CrMSIDataCubeIO::DataCube *CrMSIDataCubeIO::loadDataCube(int iCube)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  
  //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> iCube=" << iCube << std::endl; //DEBUG line
  
  DataCube *data_ptr = new DataCube();
  data_ptr->cubeID = iCube;
  data_ptr->nrows = dataCubesDesc[iCube].size();
  data_ptr->ncols = mass.length();
  data_ptr->dataOriginal = new imzMLSpectrum[data_ptr->nrows];
  data_ptr->dataInterpolated = new double*[data_ptr->nrows];

  //Data reading
  int current_imzML_id;
  int previous_imzML_id = -1; //Start previous as -1 to indicate an unallocated imzML
  for(unsigned int i = 0; i < data_ptr->nrows; i++) //For each spectrum belonging to the selected datacube
  {
    current_imzML_id = dataCubesDesc[iCube][i].imzML_ID;
    
    //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> current_imzML_id=" << current_imzML_id << std::endl; //DEBUG line!
    //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> mzLength(0)=" << imzMLReaders[current_imzML_id]->get_mzLength(0)  << std::endl; //DEBUG line
    //Rcpp::Rcout << "CrMSIDataCubeIO::loadDataCube()--> ibd file=" << imzMLReaders[current_imzML_id]->getIbdFilePath()  << std::endl; //DEBUG line
    
    
    if(current_imzML_id != previous_imzML_id)
    {
      if(previous_imzML_id != -1)
      {
        imzMLReaders[previous_imzML_id]->close();
      }
      imzMLReaders[current_imzML_id]->open();
    }
    
    data_ptr->dataInterpolated[i] = new double[data_ptr->ncols];
    data_ptr->dataOriginal[i] = imzMLReaders[current_imzML_id]->ReadSpectrum(dataCubesDesc[iCube][i].pixel_ID, //pixel id to read
                                                0, //unsigned int ionIndex
                                                mass.length(),//unsigned int ionCount
                                                data_ptr->dataInterpolated[i], //Store data directely at the datacube mem
                                                mass.length(),  //commonMassLength,
                                                mass.begin()); //double *commonMass
                                  
    previous_imzML_id = current_imzML_id;
  }
  imzMLReaders[current_imzML_id]->close(); //Force to close the last opened imzML
  
  return data_ptr;
}

void CrMSIDataCubeIO::freeDataCube(DataCube *data_ptr)
{
  for( int i = 0; i < data_ptr->nrows; i++ )
  {
    delete[] data_ptr->dataInterpolated[i];
  }
  delete[] data_ptr->dataOriginal;
  delete[] data_ptr->dataInterpolated;
  delete data_ptr;
}

void CrMSIDataCubeIO::storeDataCube(int iCube, DataCube *data_ptr) 
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  
  int current_imzML_id;
  int previous_imzML_id = -1; //Start previous as -1 to indicate an unallocated imzML
  for(unsigned int i = 0; i < data_ptr->nrows; i++) //For each spectrum belonging to the selected datacube
  {
    current_imzML_id = dataCubesDesc[iCube][i].imzML_ID;
    
    //Rcpp::Rcout << "CrMSIDataCubeIO::storeDataCube()--> current_imzML_id=" << current_imzML_id << std::endl; //DEBUG line!
    //Rcpp::Rcout << "CrMSIDataCubeIO::storeDataCube()--> mzLength(0)=" << imzMLReaders[current_imzML_id]->get_mzLength(0)  << std::endl; //DEBUG line
    //Rcpp::Rcout << "CrMSIDataCubeIO::storeDataCube()--> ibd file=" << imzMLReaders[current_imzML_id]->getIbdFilePath()  << std::endl; //DEBUG line
    
    if(current_imzML_id != previous_imzML_id)
    {
      if(previous_imzML_id != -1)
      {
        imzMLWriters[previous_imzML_id]->close();
      }
      imzMLWriters[current_imzML_id]->open();
    }
    
    //Calc average and base spectrum
    for(unsigned int j = 0; j < mass.length(); j++)
    {
      acumulatedSpectrum[current_imzML_id][j] += data_ptr->dataInterpolated[i][j];  
      baseSpectrum[current_imzML_id][j] = data_ptr->dataInterpolated[i][j] > baseSpectrum[current_imzML_id][j] ? data_ptr->dataInterpolated[i][j] : baseSpectrum[current_imzML_id][j];
    }
    
    if(imzMLWriters[current_imzML_id]->get_continuous())
    {
      //Continuous mode write
      imzMLWriters[current_imzML_id]->writeMzData(mass.length(), mass.begin()); //The mass axis will be writtem only once
      imzMLWriters[current_imzML_id]->writeIntData(mass.length(), data_ptr->dataInterpolated[i]);
    }
    else
    {
      //Processed mode write
      imzMLWriters[current_imzML_id]->writeMzData(data_ptr->dataOriginal[i].imzMLmass.size(), data_ptr->dataOriginal[i].imzMLmass.data());
      imzMLWriters[current_imzML_id]->writeIntData(data_ptr->dataOriginal[i].imzMLintensity.size(), data_ptr->dataOriginal[i].imzMLintensity.data());
    }
    
    previous_imzML_id = current_imzML_id;
  }
  imzMLReaders[current_imzML_id]->close(); //Force to close the last opened imzML
}

int CrMSIDataCubeIO::getNumberOfCubes()
{
  return dataCubesDesc.size();
}

int CrMSIDataCubeIO::getMassAxisLength()
{
  return mass.length();
}

int CrMSIDataCubeIO::getNumberOfPixelsInCube(int iCube)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  
  return dataCubesDesc[iCube].size();
}

Rcpp::DataFrame CrMSIDataCubeIO::get_OffsetsLengths(unsigned int index)
{
  if(index >= imzMLWriters.size())
  {
    throw std::runtime_error("Error: imzMLWriter index out of range\n");
  }
  return imzMLWriters[index]->get_OffsetsLengths();
}

unsigned int CrMSIDataCubeIO::get_images_count()
{
 return imzMLReaders.size();
}

int CrMSIDataCubeIO::getImageIndex(int iCube, int cubeRow)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  if( cubeRow >= dataCubesDesc[iCube].size())
  {
    throw std::runtime_error("Error: DataCube row out of range\n");
  }
  return dataCubesDesc[iCube][cubeRow].imzML_ID;
}

int CrMSIDataCubeIO::getPixelId(int iCube, int cubeRow)
{
  if(iCube >= dataCubesDesc.size())
  {
    throw std::runtime_error("Error: DataCube index out of range\n");
  }
  if( cubeRow >= dataCubesDesc[iCube].size())
  {
    throw std::runtime_error("Error: DataCube row out of range\n");
  }
  return dataCubesDesc[iCube][cubeRow].pixel_ID;
}

Rcpp::NumericVector CrMSIDataCubeIO::get_AverageSpectrum(unsigned int index)
{
  if(index >= acumulatedSpectrum.size())
  {
    throw std::runtime_error("Error: acumulatedSpectrum index out of range\n");
  }
  
  Rcpp::NumericVector AverageSpectrum(mass.length());
  for( unsigned int i = 0; i < mass.length(); i++)
  {
    AverageSpectrum[i] = acumulatedSpectrum[index][i] / ((double) (imzMLWriters[index]->get_number_of_pixels()));
  }
  
  return AverageSpectrum;
}

Rcpp::NumericVector CrMSIDataCubeIO::get_BaseSpectrum(unsigned int index)
{
  if(index >= baseSpectrum.size())
  {
    throw std::runtime_error("Error: baseSpectrum index out of range\n");
  }
  
  return baseSpectrum[index];
}