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


#ifndef RMSI_XBIN_H
#define RMSI_XBIN_H

#include <Rcpp.h>
#include <string>
#include <fstream>

class rMSIXBin
{
  public:
    rMSIXBin(Rcpp::List rMSIobject);
    ~rMSIXBin();
    
    //Return a copy of the rMSIObj
    Rcpp::List get_rMSIObj(); 
    
    //Create the ImgStream in the rMSXBin (both XML and binary parts)
    void CreateImgStream(); 
    
    //Delete any previouly created ImgStream in the rMSIXBin (both XML and binary parts)
    void DeleteImgStream();
    
    //Get a ion image in a matrix object by decoding the ImgStream
    Rcpp::NumericMatrix decodeImgStream2IonImage(unsigned int ionIndex);
      
    //TODO add methods for the rest of data encoded in rMSIXBin: normalizations... etc  
      
  private:
    Rcpp::List* rMSIObj; //Pointer to the rMSI object
    
    double* mass; //The mass axis in C format
    unsigned int massLength; //Number of mass channels
    
    unsigned int img_width, img_height; //Image size in pixels
    
    typedef struct
    {
      std::ifstream ibd_file; //the binary part of the imzML data (the XML part has been parsed already)
      bool imzMLContinuous;
      unsigned int *iX;
      unsigned int *iY;
      unsigned long *lmzLength;
      unsigned long *lmzOffset;
      unsigned long *lintLength;
      unsigned long *lintOffset;
    }imzMLHandler;
    
    typedef struct
    {
      std::string XML_file; //rMSIXBin XML file (.XrMSI)
      std::string Bin_file; //rMSIXBin Binary file (.BrMSI)
      float* fScaling; //A vector to store the ImgStream scaling factor for each m/z channel
      unsigned long * iByteLen; //ImgStream byte lengths of each encoded ion image
      unsigned long * iByteOffset; //ImgStrem byte offset of each ion image
    }rMSIXBin_Handler;
    
    rMSIXBin_Handler* _rMSIXBin; 
    
    //Encoide various Ion images to ImgStream in a single buffered execution
    void encodeMultipleIonImage2ImgStream(imzMLHandler* imzMLData, unsigned int ionIndex, unsigned int ionCount);
};

#endif
