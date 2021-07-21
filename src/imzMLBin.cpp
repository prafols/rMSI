/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2019 Pere Rafols Soler
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

#include "imzMLBin.h"
#include "mlinterp.hpp" //Used for linear interpolation
#include <stdexcept>
#include <Rcpp.h>

//#define __DEBUG__
#define MAX_MEMORY_MB 250 //Maxim usable memory

ImzMLBin::ImzMLBin(const char* ibd_fname,  unsigned int num_of_pixels,Rcpp::String Str_mzType, Rcpp::String Str_intType, bool continuous, Mode mode):
  ibdFname(ibd_fname), Npixels(num_of_pixels), bContinuous(continuous), fileMode(mode)
{
  
  mzDataType = string2imzMLDatatype(Str_mzType);
  intDataType = string2imzMLDatatype(Str_intType);
  
#ifdef __DEBUG__
  Rcpp::Rcout << "ImzMLBin::ImzMLBin() entring... with " << Npixels << " num. of pixels...";
#endif
  
  switch(mzDataType)
  {
    case int32: 
    case float32:  
      mzDataPointBytes = 4;
      break;
      
    case int64:
    case float64:
      mzDataPointBytes = 8;
      break;

    default:
      throw std::runtime_error("ERROR: ImzMLBin class constructor: invalid mzType ");
      break;
  }
  
  switch(intDataType)
  {
    case int32: 
    case float32:  
      intDataPointBytes = 4;
      break;
      
    case int64:
    case float64:
      intDataPointBytes = 8;
      break;
      
    default:
      throw std::runtime_error("ERROR: ImzMLBin class constructor: invalid intType.\n");
      break;
  }
  
  imzLength = new unsigned int[Npixels];
  lmzOffset  = new unsigned long long[Npixels];  
  iintLength = new unsigned int[Npixels];  
  lintOffset = new unsigned long long[Npixels];  
    
#ifdef __DEBUG__
Rcpp::Rcout << "ImzMLBin() constructor end successfuly\n";
Rcpp::Rcout << "imzLength addr: " << imzLength << "\n";
Rcpp::Rcout << "lmzOffset addr: " << lmzOffset << "\n";
Rcpp::Rcout << "iintLength addr: " << iintLength << "\n";
Rcpp::Rcout << "lintOffset addr: " << lintOffset << "\n";
#endif
}

ImzMLBin::~ImzMLBin()
{
  close();
  
#ifdef __DEBUG__
  Rcpp::Rcout << "imzLength addr: " << imzLength << "\n";
  Rcpp::Rcout << "lmzOffset addr: " << lmzOffset << "\n";
  Rcpp::Rcout << "iintLength addr: " << iintLength << "\n";
  Rcpp::Rcout << "lintOffset addr: " << lintOffset << "\n";
#endif
  
  delete[] imzLength;
  delete[] lmzOffset;
  delete[] iintLength;
  delete[] lintOffset;
  
#ifdef __DEBUG__
  Rcpp::Rcout << "ImzMLBin() destructor end successfuly\n";
#endif
}

const char* ImzMLBin::getIbdFilePath()
{
  return ibdFname.get_cstring();
}

bool ImzMLBin::get_continuous()
{
  return bContinuous;
}

unsigned int ImzMLBin::get_mzLength(unsigned int index)
{
  if(index >= Npixels)
  {
    throw std::runtime_error("ERROR: ImzMLBin class get_mzLength(): out of range.\n");
  }
  return imzLength[index];
}

unsigned long long ImzMLBin::get_mzOffset(unsigned int index)
{
  if(index >= Npixels)
  {
    throw std::runtime_error("ERROR: ImzMLBin class get_mzOffset(): out of range.\n");
  }
  return lmzOffset[index];
}

unsigned int ImzMLBin::get_intLength(unsigned int index)
{
  if(index >= Npixels)
  {
    throw std::runtime_error("ERROR: ImzMLBin class get_intLength(): out of range.\n");
  }
  return iintLength[index];
}

unsigned long long ImzMLBin::get_intOffset(unsigned int index)
{
  if(index >= Npixels)
  {
    throw std::runtime_error("ERROR: ImzMLBin class get_intOffset(): out of range.\n");
  }
  return lintOffset[index];
}

Rcpp::DataFrame ImzMLBin::get_OffsetsLengths()
{
  Rcpp::NumericVector RmzLengths(Npixels);
  Rcpp::NumericVector RmzOffsets(Npixels);
  Rcpp::NumericVector RintLengths(Npixels);
  Rcpp::NumericVector RintOffsets(Npixels);
  
  for(unsigned int i = 0; i< Npixels; i++)
  {
    RmzLengths[i] = imzLength[i];
    RmzOffsets[i] = lmzOffset[i];
    RintLengths[i] = iintLength[i];
    RintOffsets[i] = lintOffset[i];
  }
  
  return Rcpp::DataFrame::create( Rcpp::Named("mzLength") = RmzLengths,
                                  Rcpp::Named("mzOffset") = RmzOffsets,
                                  Rcpp::Named("intLength") = RintLengths,
                                  Rcpp::Named("intOffset") = RintOffsets
                                  );
}

void ImzMLBin::set_mzLength(Rcpp::NumericVector* mzLength_vector)
{
  if(fileMode == Mode::SequentialWriteFile)
  {
    throw std::runtime_error("ERROR: mzLength_vector cannot be set in sequential write mode.\n");
  }
  
  for(unsigned int i = 0; i < Npixels; i++)
  {
    imzLength[i] = (unsigned int)(*mzLength_vector)[i];
  }
}

void ImzMLBin::set_mzOffset(Rcpp::NumericVector* mzOffset_vector)
{
  if(fileMode == Mode::SequentialWriteFile)
  {
    throw std::runtime_error("ERROR: mzOffset_vector cannot be set in sequential write mode.\n");
  }
  
  for(unsigned int i = 0; i < Npixels; i++)
  {
    lmzOffset[i] = (unsigned long)(*mzOffset_vector)[i];
  }
}

void ImzMLBin::set_intLength(Rcpp::NumericVector* intLength_vector)
{
  if(fileMode == Mode::SequentialWriteFile)
  {
    throw std::runtime_error("ERROR: intLength_vector cannot be set in sequential write mode.\n");
  }
  
  for(unsigned int i = 0; i < Npixels; i++)
  {
    iintLength[i] = (unsigned int)(*intLength_vector)[i];
  }
}

void ImzMLBin::set_intOffset(Rcpp::NumericVector* intOffset_vector)
{
  if(fileMode == Mode::SequentialWriteFile)
  {
    throw std::runtime_error("ERROR: intOffset_vector cannot be set in sequential write mode.\n");
  }
  
  for(unsigned int i = 0; i < Npixels; i++)
  {
    lintOffset[i] = (unsigned long)(*intOffset_vector)[i];
  }
}

unsigned int ImzMLBin::get_mzEncodingBytes()
{
  return mzDataPointBytes;
}

unsigned int ImzMLBin::get_intEncodingBytes()
{
  return intDataPointBytes;
}

unsigned int ImzMLBin::get_number_of_pixels()
{
  return Npixels;  
}

void ImzMLBin::close()
{
#ifdef __DEBUG__
  Rcpp::Rcout << "ImzMLBin() destructor start...";
#endif
  if(ibdFile.is_open())
  {
    ibdFile.close();
#ifdef __DEBUG__
    Rcpp::Rcout << "ibdFile closed...\n";
#endif
  }
}

ImzMLBin::imzMLDataType ImzMLBin::string2imzMLDatatype(Rcpp::String data_type)
{
  imzMLDataType dataType;
  if( data_type == "float" )
  {
    dataType = imzMLDataType::float32;
  }
  else if( data_type == "double" )
  {
    dataType = imzMLDataType::float64;
  }
  else if( data_type == "int" )
  {
    dataType = imzMLDataType::int32;
  }
  else if( data_type == "long" )
  {
    dataType = imzMLDataType::int64;
  }
  else
  {
    throw std::runtime_error("ERROR: string2imzMLDatatype() invalid imzML datatype.\n");
  }
  
  return dataType;
}

template<typename T> 
void ImzMLBin::convertBytes2Double(char* inBytes, double* outPtr, unsigned int N)
{
  //First, copy the data to an intermediate vector with the desired type
  T* auxBuffer = new T[N];
  memcpy(auxBuffer, inBytes, sizeof(T)*N);
  
  //Finally, move the data to the double pointer
  for(int i = 0; i < N; i++)
  {
    outPtr[i] = (double)auxBuffer[i];
  }
  
  delete[] auxBuffer;
}

template<typename T> 
void ImzMLBin::convertDouble2Bytes(double* inPtr, char* outBytes, unsigned int N)
{
  //First, copy the data to an intermediate vector with the desired type
  T* auxBuffer = new T[N];
  for(int i = 0; i < N; i++)
  {
    auxBuffer[i] = (T)inPtr[i]; //Conversion from double to T type
  }
  
  //Finally, move the data to the output bytes buffer
  memcpy(outBytes, auxBuffer, sizeof(T)*N);
  
  delete[] auxBuffer;
}

ImzMLBinRead::ImzMLBinRead(const char* ibd_fname, unsigned int num_of_pixels, Rcpp::String Str_mzType, Rcpp::String Str_intType, bool continuous, bool openIbd):
  ImzMLBin(ibd_fname, num_of_pixels, Str_mzType, Str_intType, continuous, Mode::Read)
{
  if(openIbd)
  {
    open();
  }
#ifdef __DEBUG__
  Rcpp::Rcout << "ImzMLBinRead() constructor end successfuly\n";
#endif
}

ImzMLBinRead::~ImzMLBinRead()
{
  //Empty desctructor
}

void ImzMLBinRead::open()
{
#ifdef __DEBUG__
  Rcpp::Rcout << "ImzMLBinRead() open start...\nibdfile is:"<<  ibdFname.get_cstring() << "\n";
#endif
  ibdFile.open(ibdFname.get_cstring(), std::fstream::in | std::ios::binary);
  if(!ibdFile.is_open())
  {
    throw std::runtime_error("ERROR: ImzMLBinRead could not open the imzML ibd file.\n"); 
  }
}

void ImzMLBinRead::readDataCommon(unsigned long offset, unsigned int N, double* ptr, unsigned int dataPointBytes, imzMLDataType dataType)
{
  unsigned int byteCount = N*dataPointBytes;
  char* buffer = new char [byteCount];
  
  ibdFile.seekg(offset);
  if(ibdFile.eof())
  {
    throw std::runtime_error("ERROR: ImzMLBinRead reached EOF seeking the imzML ibd file.\n"); 
  }
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBinRead got fail or bad bit condition seeking the imzML ibd file.\n"); 
  }
  
  ibdFile.read (buffer, byteCount);
  if(ibdFile.eof())
  {
    throw std::runtime_error("ERROR: ImzMLBinRead reached EOF reading the imzML ibd file.\n"); 
  }
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBinRead got fail or bad bit condition reading the imzML ibd file.\n"); 
  }
  
  switch(dataType)
  {
  case int32:
    convertBytes2Double<int>(buffer, ptr, N);
    break;
    
  case float32:  
    convertBytes2Double<float>(buffer, ptr, N);
    break;
    
  case int64:
    convertBytes2Double<long>(buffer, ptr, N);
    break;
    
  case float64:
    //convertBytes2Double<double>(buffer, ptr, N); If double there is no need of intermediate conversion
    memcpy(ptr, buffer, sizeof(double)*N);
    break;
  }
  
  delete[] buffer;
}

void ImzMLBinRead::readUUID(char* uuid)
{
  ibdFile.seekg(0);
  if(ibdFile.eof())
  {
    throw std::runtime_error("ERROR: ImzMLBinRead reached EOF seeking the imzML ibd file.\n"); 
  }
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBinRead got fail or bad bit condition seeking the imzML ibd file.\n"); 
  }
  
  ibdFile.read (uuid, 16);
  if(ibdFile.eof())
  {
    throw std::runtime_error("ERROR: ImzMLBin reached EOF reading the imzML ibd file.\n"); 
  }
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBin got fail or bad bit condition reading the imzML ibd file.\n"); 
  }
}

void ImzMLBinRead::readMzData(unsigned long offset, unsigned int N, double* ptr )
{
  readDataCommon(offset, N, ptr, mzDataPointBytes, mzDataType);
}

void ImzMLBinRead::readIntData(unsigned long offset, unsigned int N, double* ptr )
{
  readDataCommon(offset, N, ptr, intDataPointBytes, intDataType);
}

//Read a single spectrum from the imzML data
//If data is in processed mode the spectrum will be interpolated to the common mass axis
//pixelID: the pixel ID of the spectrum to read.
//ionIndex: the ion index at which to start reading the spectrum (0 means reading from the begining).
//ionCount: the number of mass channels to read (massLength means reading the whole spectrum).
//out: a pointer where data will be stored.
//commonMassLength: number of points in the common mass axis.
//commonMass: pointer to the common mass axis
imzMLSpectrum ImzMLBinRead::ReadSpectrum(int pixelID, unsigned int ionIndex, unsigned int ionCount, double *out, unsigned int commonMassLength, double *commonMass)
{
  if( (ionIndex+ionCount) > commonMassLength )
  {
    throw std::runtime_error("Error: mass channels out of range\n"); 
  }
  
  imzMLSpectrum imzMLSpc;
  imzMLSpc.pixelID = pixelID;
  
  if(get_continuous())
  {
    //Continuous mode, just load the spectrum intensity vector
    readIntData(get_intOffset(pixelID) + ionIndex*get_intEncodingBytes(), ionCount, out);  
  }
  else
  {
    //Processed mode, interpolation needed
    
    //Intermediate buffers to load data before interpolation
    const int massLength = get_mzLength(pixelID);
    if( massLength != get_intLength(pixelID))
    {
      throw std::runtime_error("Error: different mass and intensity length in the imzML data\n"); 
    }
    imzMLSpc.imzMLmass.resize(massLength);
    imzMLSpc.imzMLintensity.resize(massLength);
    
    //Read processed mode mass and intensity
    try
    {
      readMzData(get_mzOffset(pixelID),  get_mzLength(pixelID),  imzMLSpc.imzMLmass.data());
      readIntData(get_intOffset(pixelID), get_intLength(pixelID), imzMLSpc.imzMLintensity.data());  
    }
    catch(std::runtime_error &e)
    {
      throw std::runtime_error(e.what());
    }
    
    //Linear interpolation
    mlinterp::interp(
      &massLength, (int)ionCount, // Number of points (imzML original, interpolated )
      imzMLSpc.imzMLintensity.data(), out, // Y axis  (imzML original, interpolated )
      imzMLSpc.imzMLmass.data(), commonMass + ionIndex // X axis  (imzML original, interpolated )
    );
    
  }
  
  return imzMLSpc;
}

ImzMLBinWrite::ImzMLBinWrite(const char* ibd_fname,  unsigned int num_of_pixels, Rcpp::String Str_mzType, Rcpp::String Str_intType, bool continuous, bool sequentialMode, bool openIbd) :
  ImzMLBin(ibd_fname, num_of_pixels, Str_mzType, Str_intType, continuous, sequentialMode? Mode::SequentialWriteFile : Mode::ModifyFile ),
  sequentialWriteIndex_IntData(0), 
  sequentialWriteIndex_MzData(0)
{
  if(openIbd)
  {
    open();
  }
}

ImzMLBinWrite::~ImzMLBinWrite()
{
  //Empty desctructor
}

void ImzMLBinWrite::open(bool truncate)
{
  if(fileMode == Mode::ModifyFile)
  {
    //Open for modifing registers in the ibd file
    ibdFile.open(ibdFname.get_cstring(), std::fstream::in | std::fstream::out | std::ios::binary); //TODO check the following comment and validate this works!
        ///According some discussion in stackoverflows... 
        ///std::ios::binary|std::ios::out|std::ios::in mus be used to replace contents in binary files... raro....
  }
  else
  {
    if(truncate)
    {
      //Open for serial writing, any content of the ibd file will be removed
      ibdFile.open(ibdFname.get_cstring(), std::fstream::out | std::ios::binary | std::fstream::trunc);
    }
    else
    {
      //Open for serial writing, append data to the end
      ibdFile.open(ibdFname.get_cstring(), std::fstream::out | std::ios::binary | std::fstream::app);
    }
  }
  
  if(!ibdFile.is_open())
  {
    throw std::runtime_error("Error: ImzMLBinWrite could not open the imzML ibd file.\n");
  }
}

void ImzMLBinWrite::writeUUID(const char* uuid)
{
  if(fileMode == Mode::ModifyFile)
  {
    throw std::runtime_error("ERROR: ibd file was opened in an invalid mode for sequencial writing and UUID can only be set in sequential mode");
  }
  
  if(ibdFile.tellp() != 0)
  {
    throw std::runtime_error("ERROR: the ibd writing pointer is not at zero so UUID cannot be writen");
  }
  
  ibdFile.write (uuid, 16);
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBinWrite got fail or bad bit condition writing the imzML ibd file.\n"); 
  }
}

void ImzMLBinWrite::writeUUID(std::string suuid)
{
  //Convert UUID in string format to byte array
  char uuid[16];
  for( int i = 0; i < 16; i++)
  {
    uuid[i] = strtol(suuid.substr(i*2, 2).c_str(), NULL, 16);
  }
  writeUUID(uuid);
}

void ImzMLBinWrite::writeMzData(unsigned long offset, unsigned int N, double* ptr )
{
  if(fileMode == Mode::SequentialWriteFile)
  {
    throw std::runtime_error("ERROR: ibd file was opened in an invalid mode for data modification");
  }
  
  writeDataCommon(offset, N, ptr, mzDataPointBytes, mzDataType);
}

void ImzMLBinWrite::writeMzData(unsigned int N, double* ptr )
{
  if(fileMode == Mode::ModifyFile)
  {
    throw std::runtime_error("ERROR: ibd file was opened in an invalid mode for sequencial writing");
  }
  
  //Store the mz offset and length
  if(sequentialWriteIndex_MzData >= Npixels)
  {
    throw std::runtime_error("ERROR: trying to write more spectral data than the maximum number of pixels set in the constructor");
  }
  lmzOffset[sequentialWriteIndex_MzData] = ibdFile.tellp();
  imzLength[sequentialWriteIndex_MzData] = N;
  sequentialWriteIndex_MzData++;
  
  //Write data
  writeDataCommon(N, ptr, mzDataPointBytes, mzDataType);
}

void ImzMLBinWrite::writeIntData(unsigned long offset, unsigned int N, double* ptr )
{
  if(fileMode == Mode::SequentialWriteFile)
  {
    throw std::runtime_error("ERROR: ibd file was opened in an invalid mode for data modification");
  }
  
  writeDataCommon(offset, N, ptr, intDataPointBytes, intDataType);
}

void ImzMLBinWrite::writeIntData(unsigned int N, double* ptr )
{
  if(fileMode == Mode::ModifyFile)
  {
    throw std::runtime_error("ERROR: ibd file was opened in an invalid mode for sequencial writing");
  }
  
  //Store the mz offset and length
  if(sequentialWriteIndex_IntData >= Npixels)
  {
    throw std::runtime_error("ERROR: trying to write more spectral data than the maximum number of pixels set in the constructor");
  }
  lintOffset[sequentialWriteIndex_IntData] = ibdFile.tellp();
  iintLength[sequentialWriteIndex_IntData] = N;
  sequentialWriteIndex_IntData++;
  
  //In continuos mode the offset and length for the mass axis is the first, so just replicate it
  if(bContinuous && sequentialWriteIndex_MzData > 0)
  {
    lmzOffset[sequentialWriteIndex_MzData] = lmzOffset[0];
    imzLength[sequentialWriteIndex_MzData] = imzLength[0];
    sequentialWriteIndex_MzData++;
  }
  
  //Write data
  writeDataCommon(N, ptr, intDataPointBytes, intDataType);
}

void ImzMLBinWrite::writeDataCommon(unsigned long offset, unsigned int N, double* ptr, unsigned int dataPointBytes, imzMLDataType dataType)
{
  ibdFile.seekp(offset);
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBinWrite got fail or bad bit condition seeking the imzML ibd file.\n"); 
  }
  
  writeDataCommon(N, ptr, dataPointBytes, dataType);
}

void ImzMLBinWrite::writeDataCommon(unsigned int N, double* ptr, unsigned int dataPointBytes, imzMLDataType dataType)
{
  unsigned int byteCount = N*dataPointBytes;
  char* buffer = new char [byteCount];
  
  //copy the ptr contents to the wrting buffer in the apropiate format 
  switch(dataType)
  {
    case int32:
      convertDouble2Bytes<int>(ptr, buffer, N);
      break;
      
    case float32:  
      convertDouble2Bytes<float>(ptr, buffer, N);
      break;
      
    case int64:
      convertDouble2Bytes<long>(ptr, buffer, N);
      break;
      
    case float64:
      //convertDouble2Bytes<double>(ptr, buffer, N); If double there is no need of intermediate conversion
      memcpy(buffer, ptr, sizeof(double)*N);
      break;
  }
  
  ibdFile.write (buffer, byteCount);
  if(ibdFile.fail() || ibdFile.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBinWrite got fail or bad bit condition writing the imzML ibd file.\n"); 
  }
  
  delete[] buffer;
}

///DEBUG METHODS////////////////////////////////////////////////////////////////////////

//' Testing the imzMLreader
//' testingimzMLBinRead
//' @param ibdFname: full path to the ibd file.
//' @param NPixels: Total number of pixels in the image.
//' @param N: number of elemetns (or data point to read).
//' @param offset: offset in bytes at which the reading operation is started.
//' @param read_mz: if true m/z data is readed, otherwise intensities are readed.
//' @param continuous: true if imzML data is in continuous mode
// [[Rcpp::export(name=".debug_imzMLBinReader")]]
Rcpp::NumericVector testingimzMLBinRead(const char* ibdFname, unsigned int NPixels, unsigned int N, unsigned long offset, Rcpp::String dataTypeString, bool read_mz, bool continuous)
{
  Rcpp::NumericVector x(N);
  try
  {
    ImzMLBinRead myReader(ibdFname, NPixels, dataTypeString, dataTypeString, continuous);
    if(read_mz)
    {
      myReader.readMzData(offset, N, x.begin());  
    }
    else
    {
      myReader.readIntData(offset, N, x.begin());
    }
    
  }
  catch(std::runtime_error &e)
  {
    Rcpp::Rcout << "Catch Error: "<< e.what() << "\n";
  }

  return x;
}

//' Testing the imzMLwriter in sequential mode
//' This function creates a new ibd file with the provided data descibed in the following params
//' @param ibdFname: full path to the ibd file.
//' @param mz_dataTypeString: String to specify the data format used to encode m/z values.
//' @param int_dataTypeString: String to specify the data format used to encode intensity values.
//' @param uuid: 16 bytes long UUID.
//' @param mzArray: A matrix with the m/z values for all pixels. Each pixel corresponds to a row. If there is only one row data will be saved in continuous mode
//' @param intArray: A matrix with the intensity values for all pixels. Each pixel corresponds to a row so the number of pixels is extracted from here.
// [[Rcpp::export(name=".debug_imzMLBinWriterSequential")]]
Rcpp::DataFrame testingimzMLBinWriteSequential(const char* ibdFname, Rcpp::String mz_dataTypeString, Rcpp::String int_dataTypeString,
                                    Rcpp::String str_uuid, Rcpp::NumericMatrix mzArray, Rcpp::NumericMatrix intArray)
{
  try
  {
    
    if(mzArray.ncol() != intArray.ncol())
    {
      throw std::runtime_error("FATAL ERROR: mass channels must have the same length as intensity data");
    }
    
    double *ptr_data = new double[mzArray.ncol()];
    
    ImzMLBinWrite myWriter(ibdFname, intArray.nrow(), mz_dataTypeString, int_dataTypeString, mzArray.nrow()==1, true);
    myWriter.writeUUID(str_uuid.get_cstring());
    
    for(int i = 0; i < intArray.nrow(); i++)
    {
      Rcpp::Rcout << "Storing... "<< i << " of " << intArray.nrow() << "\n";
      
      //Store m/z data only for the first iteration if continuous
      if(!myWriter.get_continuous() || i == 0)
      {
        for(int j = 0; j < mzArray.ncol(); j++)
        {
          ptr_data[j] = mzArray(i, j);
        }
        myWriter.writeMzData(mzArray.ncol(), ptr_data);
      }  
      
      //Store intensity data
      for(int j = 0; j < intArray.ncol(); j++)
      {
        ptr_data[j] = intArray(i, j);
      }
      myWriter.writeIntData(intArray.ncol(), ptr_data);
    }
    
    return myWriter.get_OffsetsLengths();
    
  }
  catch(std::runtime_error &e)
  {
    Rcpp::Rcout << "Catch Error: "<< e.what() << "\n";
  }
  
  return NULL;
}

//TODO i'm here, writeing two testing methods, one for sequential mode and another for modify mode. Then I'll have to create an interface in the rMSIXBin class for both modes.

//TODO data check in continuous vs. processed modes when sequential and modify writing, revise this!!!
