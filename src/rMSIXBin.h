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
    rMSIXBin(Rcpp::List rMSIobject, int nThreads);
    ~rMSIXBin();
    
    //Return a copy of the rMSIObj
    Rcpp::List get_rMSIObj(); 
    
    //Get the number number of mass channels in the common mass axis
    unsigned int get_massChannels();
    
    //Get the number of pixels
    unsigned int get_numOfPixels();
    
    //Create the ImgStream in the rMSXBin (both XML and binary parts). Any previois rMSXBin files will be deleted!
    void CreateImgStream(); 
    
    //Get multiple ion image in a matrix object by decoding the ImgStream
    //The MAX operator will be used to merge all ion images in a single image matrix
    Rcpp::NumericMatrix decodeImgStream2IonImages(unsigned int ionIndex, unsigned int ionCount, Rcpp::NumericVector normalization_coefs);
      
    //Calculate average spectrum, base spectrum and normalizations
    //imzMLreader: pointer to an imzMLreader initialized object
    void CalculateAverageBaseNormalizations(ImzMLBinRead *imzMLreader);
    
    //Read a single spectrum from the imzML data
    //If data is in processed mode the spectrum will be interpolated to the common mass axis
    //imzMLreader: pointer to an imzMLreader initialized object
    //pixelID: the pixel ID of the spectrum to read.
    //ionIndex: the ion index at which to start reading the spectrum (0 means reading from the begining).
    //ionCount: the number of mass channels to read (massLength means reading the whole spectrum).
    //out: a pointer where data will be stored (must be allocated before calling this function).
    void ReadSpectrum(ImzMLBinRead *imzMLreader, int pixelID, unsigned int ionIndex, unsigned int ionCount, double *out);
    
  private:
    unsigned int irMSIFormatVersion; //An integer to record the rMSI format version
    std::string sImgName; //A string to record the MS image name.
    Rcpp::List rMSIObj; //Internal copy of the rMSI object
    char UUID_imzML[16]; //The linked imzML UUID in raw bytes
    std::string sUUID_imzML; //The linked imzML UUID as a string
    char UUID_rMSIXBin[16]; //The rMSIXBin UUID in raw bytes
    std::string sUUID_rMSIXBin; //The linked imzML UUID as a string
    
    Rcpp::NumericVector massAxis; //The common mass axis
    double pixel_size_um; //pixel resolution in microns
    unsigned int img_width, img_height; //Image size in pixels
    
    unsigned int number_of_encoding_threads; //Max number of threads to use for the imgStream encoding
    
    typedef struct
    {
      unsigned int numOfPixels; //Total number of pixel in the image;
      std::string XML_file; //rMSIXBin XML file (.XrMSI)
      std::string Bin_file; //rMSIXBin Binary file (.BrMSI)
      unsigned long* iByteLen; //ImgStream byte lengths of each encoded ion image
      unsigned long* iByteOffset; //ImgStrem byte offset of each ion image
      unsigned int* iX; //Corrected X coordinates (non motor coords)
      unsigned int* iY; //Corrected Y coordinates (non motor coords)
      unsigned long* normByteOffsets; //the offset of the start of eahc normalization vector in the binary file.
      std::vector<std::string> normNames; //A vector to store normalization names during XML parsing
    }rMSIXBin_Handler;
    
    rMSIXBin_Handler* _rMSIXBin; 
    
    typedef struct
    {
      unsigned int ionIndex; //Ion index of the current encoded image
      float scaling; //The scaling factor of an ion image
      std::vector<unsigned char> png_stream; //the encoded png stream
    }ImgStreamEncoder_result;
    
    //Threaded encoding model 
    ImgStreamEncoder_result encodeBuffer2SingleImgStream(double *buffer, unsigned int ionIndex, unsigned int bufferIonIndex, unsigned int bufferIonCount); //Threaded method
    void startThreadedEncoding(double *buffer, unsigned int ionIndex, unsigned int ionCount); //Threaded method
    
    //Threaded normalizations and average spectrum
    void startThreadedAverageBaseNormalizations(double *intensity, int pixel, 
                                                Rcpp::NumericVector *averageSpectrum, Rcpp::NumericVector *baseSpectrum,
                                                Rcpp::NumericVector *normTIC, Rcpp::NumericVector *normRMS, Rcpp::NumericVector *normMAX );
    
    //Store normalization vectors
    void storeNormalizations2Binary();
    
    //Load normalization vectors
    void loadNormalizationFromBinary();
    
    //Get the byte representation from a 16 bytes uuid string
    void hexstring2byteuuid(std::string hex_str, char* output);
    
    //Write the XML file, any previous .XrMSI file will be deleted
    bool writeXrMSIfile();
    
    //Load a XML file
    void readXrMSIfile();
    
    //Copy imgStream to the rMSIObject
    void copyimgStream2rMSIObj();
    
    //read the BrMSI uuid, mass axis, average spectrum and base spectrum
    void readBrMSI_header();
    
};

#endif
