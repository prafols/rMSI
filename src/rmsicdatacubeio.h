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
#ifndef RMSI_DATA_CUBE_IO_H
  #define RMSI_DATA_CUBE_IO_H
#include <Rcpp.h>
#include "imzMLBin.h"

/********************************************************************************
 *  CrMSIDataCubeIO: A C++ class to extract datacubes from multiple imzML files
 *  
 *  This class provides a simple interface to load and save data from imzML files
 *  using the rMSIproc multithread processing approach.
 ********************************************************************************/
class CrMSIDataCubeIO
{
  public:
    CrMSIDataCubeIO(Rcpp::NumericVector massAxis, double cubeMemoryLimitMB);
    ~CrMSIDataCubeIO();
    
    //Struct to define a whole data cube in memory
    typedef struct
    {
      int cubeID;
      int ncols;
      int nrows;
      double **data;  
    } DataCube;
    
    //Set de data output path. New imzML files will be created with the processed data. If not called, no data is saved aside of the returning values of the processing function.
    void setDataOutputPath(const char* out_path);
    
    //Appends an image to be processed.
    // The input argument is the rMSI object describing a complete MSI image.
    void appedImageData(Rcpp::List rMSIOoj);
    
    //Loads a data cube specified by iCube into data_ptr
    //WARNING: This is not a thread-safe method, it must be only used on the main thread who'll take care of loading data from HDD and copy it to other thread mem space.
    //It returns a pointer to an allocated structure containing the datacube.
    //It is responisability of user to free memmory of the loaded datacube using freeDataCube() function.
    DataCube *loadDataCube( int iCube);
    void freeDataCube(DataCube *data_ptr);
    
    //Stores a datacube to the path assosiated with its ID
    void storeDataCube(int iCube, DataCube *data_ptr);
    
    //Return the total number of cubes in the ramdisk
    int getNumberOfCubes();
    
    //Return the length of the common mass axis
    int getMassAxisLength();
    
    //Return the total number of pixels in each cube
    int getNumberOfPixelsInCube(int iCube);

  private:
    bool storeData; //A bool indcatinc whether processed data must be stored or not. By default it is false and the function setDataOutputPath() must be called to set it to true.
    unsigned int cubeMaxNumRows; //The maximum rows in a datacube calculated from the maximum memory allowed by each cube and the mass axis length.
    Rcpp::NumericVector mass; //A common mass axis for all images to process
    std::vector<ImzMLBinRead*> imzMLReaders; //Pointers to multiple imzMLReadrs initialized with openIbd = false to avoid exiding the maximum open files.
    
    //Struct to internally handle data cube accessors
    typedef struct
    {
      int pixel_ID; //ID of a pixel in the imzML file
      int imzML_ID; //integer pointing to a imzML Handle (reader and writer) for a pixel ID, -1 value is allowed to indicate to an unallocated imzML
    } PixelDescription;
    
    typedef std::vector<PixelDescription> DataCubeDescription; //Each data cube is defined as standard vector of pixel descriptions
    
    std::vector<DataCubeDescription> dataCubesDesc; //A vector to describe all data cubes in the data set
};

#endif
