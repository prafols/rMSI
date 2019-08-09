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
#include <Rcpp.h>

enum imzMLDataType
{
  int32, //32 bits integer  (int)
  int64, //64 bits interger (long)
  float32, //32 bits float  (float)
  float64  //64 bits float  (double)
};

class ImzMLBin
{
  public:
    ImzMLBin(const char* ibd_fname, unsigned int num_of_pixels, imzMLDataType mzType, imzMLDataType intType, bool continuous);
    ~ImzMLBin();
    
    bool get_continuous();
    unsigned int get_mzLength(unsigned int index);
    unsigned long get_mzOffset(unsigned int index);
    unsigned int get_intLength(unsigned int index);
    unsigned long get_intOffset(unsigned int index);
    
    void set_mzLength(Rcpp::NumericVector* mzLength_vector);
    void set_mzOffset(Rcpp::NumericVector* mzOffset_vector);
    void set_intLength(Rcpp::NumericVector* intLength_vector);
    void set_intOffset(Rcpp::NumericVector* intOffset_vector);
    
  protected:
    std::fstream ibdFile; 
    unsigned int Npixels; //Total number of pixels in the image;
    unsigned int mzDataPointBytes; //Number of bytes used to encode a mass channel 
    unsigned int intDataPointBytes; //Number of bytes used to encode an intensity data point
    imzMLDataType mzDataType, intDataType;
    bool bContinuous;
    unsigned int* imzLength;
    unsigned long* lmzOffset;
    unsigned int* iintLength;
    unsigned long* lintOffset;
    
    template<typename T> 
    void covertBytes2Double(char* inBytes, double* outPtr, unsigned int N);
};

class ImzMLBinRead : public ImzMLBin
{
  public: 
    ImzMLBinRead(const char* ibd_fname, unsigned int num_of_pixels, imzMLDataType mzType, imzMLDataType intType, bool continuous);
    ~ImzMLBinRead();
    
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
    ImzMLBinWrite(const char* ibd_fname,  unsigned int num_of_pixels, imzMLDataType mzType, imzMLDataType intType, bool continuous);
    ~ImzMLBinWrite();
    
    //Write N elements from the ibd file at the given offset as m/z channels
    //Data is obtained from ptr pointer
    void writeMzData(unsigned long offset, unsigned int N, double* ptr );
    
    //Write N elements from the ibd file at the given offset as spectrum intensities
    //Data is obtained from ptr pointer
    void writeIntData(unsigned long offset, unsigned int N, double* ptr );
};

#endif
