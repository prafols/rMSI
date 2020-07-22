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
#include "imzMLBin.h"

//TODO decidir quants bits em quedo... amb 16bits el brain TOF son 524MB amb 8bits son 327MB, no hi ha tanta differencia?
//TODO per decidir quants bits son necessaris fes una scala d color rainbow amb diferents tramps: 256, 2^12, 2^16...etc... l'ull hi veu diferencies?
#define IMG_STREAM_8bits //Comment this line out if you prefere to encode the imgStream as 16bits plus a mask

#define IONIMG_BUFFER_MB 250 //I think 250 MB of RAM is a good balance for fast hdd operation and low memory footprint

#ifdef IMG_STREAM_8bits
typedef  unsigned char imgstreamencoding_type;
#define ENCODER_RANGE 255.0
#define ENCODING_BITS 8
#define ENCODING_BIT_MASK 0xFF //imgStream encoding mask (in 8 bit encoding, no masking)
#else
typedef  unsigned short imgstreamencoding_type;
#define ENCODER_RANGE 65535.0
#define ENCODING_BITS 16
#define ENCODING_BIT_MASK 0xFFC0 //imgStream encoding mask (only used for 16 bit encoding)
#endif

class rMSIXBin
{
  public:
    
    //constructor using a .XrMSI file
    rMSIXBin(Rcpp::String path, Rcpp::String fname);
    
    //Constructor using an already filled rMSIobject
    rMSIXBin(Rcpp::List rMSIobject);
    ~rMSIXBin();
    
    //Return a copy of the rMSIObj
    Rcpp::List get_rMSIObj(); 
    
    //Create the ImgStream in the rMSXBin (both XML and binary parts). Any previois rMSXBin files will be deleted!
    void CreateImgStream(); 
    
    //Get multiple ion image in a matrix object by decoding the ImgStream
    //The MAX operator will be used to merge all ion images in a single image matrix
    Rcpp::NumericMatrix decodeImgStream2IonImages(unsigned int ionIndex, unsigned int ionCount);
      
    //TODO add methods for the rest of data encoded in rMSIXBin: normalizations... etc  puc fer un spl metode per passar tot un imzML i calcular:
    //(Tot aixo hauria d'anar a la classe del imzML ja que u fa amb imzML)
        // - Espectre mig
        // - Espectre max o skyline 
        // - normalitzacions
        
        
      
  private:
    unsigned int irMSIFormatVersion; //An integer to record the rMSI format version
    std::string sImgName; //A string to record the MS image name.
    Rcpp::List rMSIObj; //Internal copy of the rMSI object
    char UUID_imzML[16]; //The linked imzML UUID in raw bytes
    std::string sUUID_imzML; //The linked imzML UUID as a string
    char UUID_rMSIXBin[16]; //The rMSIXBin UUID in raw bytes
    std::string sUUID_rMSIXBin; //The linked imzML UUID as a string
    
    double* mass; //The mass axis in C format
    unsigned int massLength; //Number of mass channels
    double pixel_size_um; //pixel resolution in microns
    unsigned int img_width, img_height; //Image size in pixels
    
    typedef struct
    {
      unsigned int numOfPixels; //Total number of pixel in the image;
      std::string XML_file; //rMSIXBin XML file (.XrMSI)
      std::string Bin_file; //rMSIXBin Binary file (.BrMSI)
      unsigned long* iByteLen; //ImgStream byte lengths of each encoded ion image
      unsigned long* iByteOffset; //ImgStrem byte offset of each ion image
      unsigned int* iX; //Corrected X coordinates (non motor coords)
      unsigned int* iY; //Corrected Y coordinates (non motor coords)
    }rMSIXBin_Handler;
    
    rMSIXBin_Handler* _rMSIXBin; 
    
    //Encoide various Ion images to ImgStream in a single buffered execution
    void encodeMultipleIonImage2ImgStream(ImzMLBinRead* imzMLHandler, unsigned int ionIndex, unsigned int ionCount);
    void encodeMultipleIonImage2ImgStream_continuous(ImzMLBinRead* imzMLHandler, unsigned int ionIndex, unsigned int ionCount);
    void encodeMultipleIonImage2ImgStream_processed(ImzMLBinRead* imzMLHandler, unsigned int ionIndex, unsigned int ionCount);
    
    //Get the byte representation from a 16 bytes uuid string
    void hexstring2byteuuid(std::string hex_str, char* output);
    
    //Display the current encoded persentage on console
    void coutEncodingPersentage(unsigned int ionIndex);
    
    //Write the XML file, any previous .XrMSI file will be deleted
    bool writeXrMSIfile();
    
    //Load a XML file
    void readXrMSIfile();
    
    //Copy imgStream to the rMSIObject
    void copyimgStream2rMSIObj();
    
};

#endif
