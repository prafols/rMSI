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

enum imzMLDataType
{
  int32, //32 bits integer  (int)
  int64, //64 bits interger (long)
  float32, //32 bits float  (float)
  flaot64  //64 bits float  (double)
};

class ImzMLBin
{
  public:
    ImzMLBin(const char* ibd_fname, imzMLDataType mzType, imzMLDataType intType);
    ~ImzMLBin();
    
  protected:
    std::fstream ibdFile; 
    unsigned int mzDataPointBytes; //Number of bytes used to encode a mass channel 
    unsigned int intDataPointBytes; //Number of bytes used to encode a intensity data point
    imzMLDataType mzDataType, intDataType;
    
    template<typename T> 
    void covertBytes2Double(char* inBytes, double* outPtr, unsigned int N);
};

class ImzMLBinRead : public ImzMLBin
{
  public: 
    ImzMLBinRead(const char* ibd_fname, imzMLDataType mzType, imzMLDataType intType);
    ~ImzMLBinRead();
    
    //Read N elements from the ibd file at the given offset as m/z channels
    //Data is sotred at the ptr pointer
    void readMzData(unsigned long offset, unsigned int N, double* ptr );
    
    //Read N elements from the ibd file at the given offset as spectrum intensities
    //Data is sotred at the ptr pointer
    void readIntData(unsigned long offset, unsigned int N, double* ptr );
};

class ImzMLBinWrite : public ImzMLBin
{
  public: 
    ImzMLBinWrite(const char* ibd_fname, imzMLDataType mzType, imzMLDataType intType);
    ~ImzMLBinWrite();
    
    //Write N elements from the ibd file at the given offset as m/z channels
    //Data is obtained from ptr pointer
    void writeMzData(unsigned long offset, unsigned int N, double* ptr );
    
    //Write N elements from the ibd file at the given offset as spectrum intensities
    //Data is obtained from ptr pointer
    void writeIntData(unsigned long offset, unsigned int N, double* ptr );
};

#endif

