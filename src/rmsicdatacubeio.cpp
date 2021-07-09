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


CrMSIDataCubeIO::CrMSIDataCubeIO(NumericVector massAxis, double cubeMemoryLimitMB)
:mass(massAxis)
{
  storeData = false;
  cubeMaxNumRows = std::ceil((1024*1024*cubeMemoryLimitMB)/(double)(8*mass.length()));
}

CrMSIDataCubeIO::~CrMSIDataCubeIO()
{
  for( int i = 0; i < imzMLReaders.size(); i++)
  {
    delete imzMLReaders[i];
  }
}

void CrMSIDataCubeIO::setDataOutputPath(const char* out_path)
{
  storeData = true;
  //TODO create the imzML writers... well maybe I don't have the info yet since data is not appended... then just store the path? think about it
  throw std::runtime_error("Error: TODO setDataOutputPath() not implemented yet!\n");
}

void CrMSIDataCubeIO::appedImageData(Rcpp::List rMSIobj)
{
  //Set the imzML reader
  List data = rMSIobj["data"];
  List imzML = data["imzML"];
  DataFrame imzMLrun = as<DataFrame>(imzML["run"]);
  std::string sFilePath = as<std::string>(data["path"]);
  std::string sFnameImzML = as<std::string>(imzML["file"]);
  sFnameImzML= sFilePath + "/" + sFnameImzML + ".ibd";
  
  //Rcpp::Rcout << "CrMSIDataCubeIO::appedImageData()--> sFnameImzML = "<< sFnameImzML << std::endl; //DEBUG line

  //Append the imzML Reader
  imzMLReaders.push_back(new ImzMLBinRead(sFnameImzML.c_str(), 
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
  
  //TODO implement me!
  throw std::runtime_error("Error: storeDataCube() is not implemented yet!\n");
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
