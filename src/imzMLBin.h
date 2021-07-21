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


#ifndef IMZML_BIN_H
#define IMZML_BIN_H

#include <fstream>
#include <vector>
#include <Rcpp.h>

typedef struct
{
  int pixelID;
  std::vector<double> imzMLmass; 
  std::vector<double> imzMLintensity; 
}imzMLSpectrum; //If data is in continuous mode the std::vectors will be empty

class ImzMLBin
{
  public:
    enum imzMLDataType
    {
      int32, //32 bits integer  (int)
      int64, //64 bits interger (long)
      float32, //32 bits float  (float)
      float64  //64 bits float  (double)
    } ;
    
    enum Mode { Read, SequentialWriteFile, ModifyFile }; //The mode must be spcified in the constructor
    
    ImzMLBin(const char* ibd_fname, unsigned int num_of_pixels, Rcpp::String Str_mzType, Rcpp::String Str_intType, bool continuous, Mode mode);
    ~ImzMLBin();
    
    const char* getIbdFilePath();
    bool get_continuous();
    unsigned int get_mzLength(unsigned int index);
    unsigned long long get_mzOffset(unsigned int index);
    unsigned int get_intLength(unsigned int index);
    unsigned long long get_intOffset(unsigned int index);
    Rcpp::DataFrame get_OffsetsLengths(); //Get all offsets and legnth in a R data frame
    
    void set_mzLength(Rcpp::NumericVector* mzLength_vector);
    void set_mzOffset(Rcpp::NumericVector* mzOffset_vector);
    void set_intLength(Rcpp::NumericVector* intLength_vector);
    void set_intOffset(Rcpp::NumericVector* intOffset_vector);
    
    unsigned int get_number_of_pixels();
    
    //Get the number of bytes used for encoding
    unsigned int get_mzEncodingBytes();
    unsigned int get_intEncodingBytes();
    
    //Close the file connection
    void close();
    
  protected:
    Mode fileMode; //Define the mode used to acces the binary file
    Rcpp::String ibdFname;
    std::fstream ibdFile; 
    unsigned int Npixels; //Total number of pixels in the image;
    unsigned int mzDataPointBytes; //Number of bytes used to encode a mass channel 
    unsigned int intDataPointBytes; //Number of bytes used to encode an intensity data point
    imzMLDataType mzDataType, intDataType;
    bool bContinuous;
    unsigned int* imzLength;
    unsigned long long* lmzOffset;
    unsigned int* iintLength;
    unsigned long long* lintOffset;
    
    //Get the imzMLDataType from a string
    imzMLDataType string2imzMLDatatype(Rcpp::String data_type);
    
    template<typename T> 
    void convertBytes2Double(char* inBytes, double* outPtr, unsigned int N);
    
    template<typename T> 
    void convertDouble2Bytes(double* inPtr, char* outBytes, unsigned int N);
};

class ImzMLBinRead : public ImzMLBin
{
  public: 
    ImzMLBinRead(const char* ibd_fname, unsigned int num_of_pixels, Rcpp::String Str_mzType, Rcpp::String Str_intType, bool continuous, bool openIbd = true);
    ~ImzMLBinRead();
    
    //Open the ibd file in reading mode
    void open();
    
    //Get the 16 bytes UUID from the imzML ibd file. uuid must be allocated by the user.
    void readUUID(char* uuid);
    
    //Read N elements from the ibd file and decode them as m/z data.
    //offset: offset in bytes at which the reading operation is started.
    //N: number of elements to read from the ibd file (N is elements, not bytes!)
    //ptr: Data will be stored at the ptr pointer
    void readMzData(unsigned long offset, unsigned int N, double* ptr );
    
    //Read N elements from the ibd file and decode them as intensity data.
    //offset: offset in bytes at which the reading operation is started.
    //N: number of elements to read from the ibd file (N is elements, not bytes!)
    //ptr: Data will be stored at the ptr pointer
    void readIntData(unsigned long offset, unsigned int N, double* ptr );
    
    //Read a single spectrum from the imzML data
    //If data is in processed mode the spectrum will be interpolated to the common mass axis
    //pixelID: the pixel ID of the spectrum to read.
    //ionIndex: the ion index at which to start reading the spectrum (0 means reading from the begining).
    //ionCount: the number of mass channels to read (massLength means reading the whole spectrum).
    //out: a pointer where data will be stored.
    //commonMassLength: number of points in the common mass axis.
    //commonMass: pointer to the common mass axis
    imzMLSpectrum ReadSpectrum(int pixelID, unsigned int ionIndex, unsigned int ionCount, double *out, unsigned int commonMassLength, double *commonMass);
    
  private:
    //Read N elements from the ibd file and decode them.
    //offset: offset in bytes at which the reading operation is started.
    //N: number of elements to read from the ibd file (N is elements, not bytes!)
    //ptr: Data will be stored at the ptr pointer
    //dataPointBytes: number of bytes used to encode a single data point.
    //dataType: data type used for the encoding.
    void readDataCommon(unsigned long offset, unsigned int N, double* ptr, unsigned int dataPointBytes, imzMLDataType dataType);
};

class ImzMLBinWrite : public ImzMLBin
{
  public: 
    ImzMLBinWrite(const char* ibd_fname,  unsigned int num_of_pixels, Rcpp::String Str_mzType, Rcpp::String Str_intType, bool continuous, bool sequentialMode, bool openIbd = true);
    ~ImzMLBinWrite();
    
    //Open the ibd file in writing mode
    //If trucate is set to true the file will be completelly removed and started over again. thats the case when writing the uuid
    void open(bool truncate = false);
    
    //Write the 16 bytes UUID to the imzML ibd file.
    void writeUUID(const char* uuid);
    
    //Write the 16 bytes UUID to the imzML ibd file. Uuid provided as a std::String
    void writeUUID(std::string suuid);
    
    //Write N elements to the ibd file at the given offset as m/z channels
    //Data is obtained from ptr pointer
    void writeMzData(unsigned long offset, unsigned int N, double* ptr );
    
    //Append N elements to the ibd file (write in sequential mode)
    //Data is obtained from ptr pointer
    void writeMzData( unsigned int N, double* ptr );
    
    //Write N elements to the ibd file at the given offset as spectrum intensities
    //Data is obtained from ptr pointer
    void writeIntData(unsigned long offset, unsigned int N, double* ptr );
    
    //Append N elements to the ibd file (write in sequential mode)
    //Data is obtained from ptr pointer
    void writeIntData( unsigned int N, double* ptr );
    
  private:
    unsigned int sequentialWriteIndex_MzData; //When sequentially writing data, this integers provides the index of the next pixel to store
    unsigned int sequentialWriteIndex_IntData; //When sequentially writing data, this integers provides the index of the next pixel to store
    
    //Write N elements to the ibd file encoded in the specified format.
    //This method is tailored to data modification mode
    //offset: offset in bytes at which the writing operation is started.
    //N: number of elements to write to the ibd file (N is elements, not bytes!)
    //ptr: Pointer to the data to write
    //dataPointBytes: number of bytes used to encode a single data point.
    //dataType: data type used for the encoding.
    void writeDataCommon(unsigned long offset, unsigned int N, double* ptr, unsigned int dataPointBytes, imzMLDataType dataType);
    
    //Write N elements to the ibd file encoded in the specified format.
    //This method is tailored to sequential writing mode
    //N: number of elements to write to the ibd file (N is elements, not bytes!)
    //ptr: Pointer to the data to write
    //dataPointBytes: number of bytes used to encode a single data point.
    //dataType: data type used for the encoding.
    void writeDataCommon(unsigned int N, double* ptr, unsigned int dataPointBytes, imzMLDataType dataType);
};

#endif
