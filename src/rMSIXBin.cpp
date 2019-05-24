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

#include <iostream>
#include <fstream>
#include "rMSIXBin.h"
#include "lodepng.h"

#define IONIMG_BUFFER_MB 250 //I think 250 MB of RAM is a good balance for fast hdd operation and low memory footprint
#define IMG_ENCODING_BYTES 2 //I think 16 bits provides with enought dynamic range (remember I also have the scaling factor)
  
rMSIXBin::rMSIXBin(Rcpp::List rMSIobject)
{
  rMSIObj = new Rcpp::List();
  Rcpp::List* data = new Rcpp::List();
  Rcpp::List* rMSIXBin = new Rcpp::List();
  
  *rMSIObj = rMSIobject; 
  *data=(*rMSIObj)["data"]; 
  *rMSIXBin = (*data)["rMSIXBin"];
  
  Rcpp::NumericVector imgSize = (*rMSIObj)["size"]; 
  img_width = (unsigned int) imgSize["x"];
  img_height = (unsigned int) imgSize["y"];
  
  //File handlers to rMSIXBin
  _rMSIXBin = new rMSIXBin_Handler;
  std::string sFnameImgStream = Rcpp::as<std::string>((*rMSIXBin)["file"]);
  _rMSIXBin->XML_file = sFnameImgStream + ".XrMSI";
  _rMSIXBin->Bin_file = sFnameImgStream + ".BrMSI";
  
  //Get the mass axis
  Rcpp::NumericVector massAxis = (*data)["mass"];
  massLength = massAxis.length();
  mass = new double[massLength];
  memcpy(mass, massAxis.begin(), sizeof(double)*massLength);
  
  //Init the ImgStream offset/pointers
  _rMSIXBin->fScaling = new float[massLength];
  _rMSIXBin->iByteLen = new unsigned long[massLength]; //using long instead of int to allow extra room to store all offsets!
  _rMSIXBin->iByteOffset = new unsigned long[massLength]; //using long instead of int to allow extra room to store all offsets!
  
  //TODO: if there is already a ImgStream fill the vectors: fScaling, iByteLen, iByteOffset using XML data...
  //TODO I need to place the XML stuf in a new XML parser method
  
  
  delete rMSIXBin;
  delete data;
}

rMSIXBin::~rMSIXBin()
{
  delete[] mass;
  delete[] _rMSIXBin->fScaling;
  delete[] _rMSIXBin->iByteLen;
  delete[] _rMSIXBin->iByteOffset;
  delete _rMSIXBin;
  delete rMSIObj;
}

Rcpp::List rMSIXBin::get_rMSIObj()
{
  return *rMSIObj;
}


//TODO: What if a imgstream already exists in the binary file?
//TODO: This is just the imgstream, other methods must be created to fill the rest of rMSIXBin data
void rMSIXBin::CreateImgStream()
{
  
  Rcpp::List* data = new Rcpp::List();
  
  Rcpp::List* imzML = new Rcpp::List();
  Rcpp::DataFrame* imgStream = new Rcpp::DataFrame(); 
  Rcpp::DataFrame* imzMLrun = new Rcpp::DataFrame();
  
  Rcpp::NumericVector* Xcoords = new Rcpp::NumericVector();
  Rcpp::NumericVector* Ycoords = new Rcpp::NumericVector();
  Rcpp::NumericVector* imzML_mzLength = new Rcpp::NumericVector();
  Rcpp::NumericVector* imzML_mzOffsets = new Rcpp::NumericVector();
  Rcpp::NumericVector* imzML_intLength = new Rcpp::NumericVector();
  Rcpp::NumericVector* imzML_intOffsets = new Rcpp::NumericVector();
  
  (*data)=(*rMSIObj)["data"]; 
  *imzML = (*data)["imzML"];
  
  imzMLHandler* imzMLData = new imzMLHandler;
  imzMLData->imzMLContinuous = (*imzML)["continuous_mode"];
  std::string sFnameImzML = Rcpp::as<std::string>((*imzML)["file"]);
  sFnameImzML= sFnameImzML + ".ibd";
  
  *imzMLrun = (*imzML)["run"];
  *Xcoords = (*imzMLrun)["x"];
  *Ycoords = (*imzMLrun)["y"];
  *imzML_mzLength = (*imzMLrun)["mzLength"];
  *imzML_mzOffsets = (*imzMLrun)["mzOffset"];
  *imzML_intLength = (*imzMLrun)["intLength"];
  *imzML_intOffsets = (*imzMLrun)["intOffset"];
  
  //Copy the imzMLrun to C data types
  imzMLData->iX = new unsigned int[(*imzMLrun).nrows()];
  imzMLData->iY = new unsigned int[(*imzMLrun).nrows()];
  imzMLData->lmzLength  = new unsigned long[(*imzMLrun).nrows()];
  imzMLData->lmzOffset  = new unsigned long[(*imzMLrun).nrows()];
  imzMLData->lintLength = new unsigned long[(*imzMLrun).nrows()];
  imzMLData->lintOffset = new unsigned long[(*imzMLrun).nrows()];
  
  for(int i = 0; i < (*imzMLrun).nrows(); i++)
  {
    imzMLData->iX[i] = (unsigned int)(*Xcoords)[i];
    imzMLData->iY[i] = (unsigned int)(*Ycoords)[i];
    imzMLData->lmzLength[i] = (unsigned long)(*imzML_mzLength)[i];
    imzMLData->lmzOffset[i] = (unsigned long)(*imzML_mzOffsets)[i];
    imzMLData->lintLength[i] = (unsigned long)(*imzML_intLength)[i];
    imzMLData->lintOffset[i] = (unsigned long)(*imzML_intOffsets)[i];
  }
  
  //All R data related to imzML have been already copied to C data, so delete them
  delete Xcoords;
  delete Ycoords;
  delete imzML_mzLength;
  delete imzML_mzOffsets;
  delete imzML_intLength;
  delete imzML_intOffsets;
  delete imzML;
  delete imzMLrun;
  

  
  /* iIonImgCount calculation
  *  
  *  bytesPerIonImg = img_width * img_height * IMG_ENCODING_BYTES + 4 (32bits float scalingFactor)
  *  iIonImgCount = IONIMG_BUFFER_MB * 1024 * 1024 / bytesPerIonImg
  */
  unsigned int iIonImgCount = (unsigned int)(  ((double)(IONIMG_BUFFER_MB * 1024 * 1024)) / ((double)( img_width *img_height * IMG_ENCODING_BYTES + 4 )) );

  //Prepare file connections
  imzMLData->ibd_file.open(sFnameImzML.c_str(), std::ios::in|std::ios::binary);
  if (imzMLData->ibd_file.is_open())
  {
    //TODO open the rMSIXBin binary file here and check for errors 
    //Loop to create each ion image
    for( int i = 0; i < massLength; i = i + iIonImgCount)
    {
      encodeMultipleIonImage2ImgStream(imzMLData, i, iIonImgCount); 
    }
    //TODO close the rMSIXBin files
    imzMLData->ibd_file.close();
  }
  else
  {
    Rcpp::Rcout << "Error: imzML not available" << std::endl; //TODO improve this error handling in case of not having imzML data
  }
  
  
  //TODO make some Rcout to check some internal data as for example the iIonImgCount
  
  //TODO copy data back to imgStream... 
  
  
  
  //TODO: pensar funcionament de la part XML del rMSIXBin, crec que el millor es el seguent:
   //1. El constructor parseja el fixer XML original (si existeix i guarda les dades)
   //2. Es modifiquen les dades correponents al ImgStream
   //3. Es sobrescriu el XML sencer amb les dades originals+les noves dades del ImgStream
  
  delete imgStream;
  delete data;
  
  delete[] imzMLData->iX;
  delete[] imzMLData->iY;
  delete[] imzMLData->lmzLength;
  delete[] imzMLData->lmzOffset;
  delete[] imzMLData->lintLength;
  delete[] imzMLData->lintOffset;
  delete imzMLData;
}


void rMSIXBin::DeleteImgStream()
{
  //TODO: delete just the binary part correponding to the imgStream
  //TODO: delete just the XML part correponding to the imgStream
  
}

void rMSIXBin::encodeMultipleIonImage2ImgStream(imzMLHandler* imzMLData, unsigned int ionIndex, unsigned int ionCount)
{
  if(imzMLHandler->imzMLContinuous)
  {
    encodeMultipleIonImage2ImgStream_continuous(imzMLData, ionIndex, ionCount);
  }
  else
  {
    encodeMultipleIonImage2ImgStream_processed(imzMLData, ionIndex, ionCount);
  }
}

void rMSIXBin::encodeMultipleIonImage2ImgStream_continuous(imzMLHandler* imzMLData, unsigned int ionIndex, unsigned int ionCount)
{
   //TODO implemnet me!
  //Step1, get the complete buffer from imzML: seek, binread, seek, binread... 
  
  //Step2, call the single image encoding function
}

void rMSIXBin::encodeMultipleIonImage2ImgStream_processed(imzMLHandler* imzMLData, unsigned int ionIndex, unsigned int ionCount)
{
  //TODO implement the imzML processed mode methods
  Rcpp::Stop("TODO: The imzML processed mode is not implemented yet, sorry.");
}

void encodeSingleIonImage2ImgStream()
{
  //TODO implementar, pensar paramteres i possaro dins la classe rMSIXBin
}


//' decodePngStream2IonImage.
//'
//' Obtain a single mass channel ion image be decoding the hdd png stream at a specified ionIndex.
//'
//' @param ionIndex the index of ion to extract from the png stream. //TODO starting by 0 o 1? Im going to call this from R or from C++?
//' 
//' @return A NumerixMatrix containing the ion image. //TODO think if row and cols correspond to witdh and height or visaversa
//' 
Rcpp::NumericMatrix rMSIXBin::decodeImgStream2IonImage(unsigned ionIndex)
{
  //TODO falta normalitzations o scaling factor!!!
  
  //TODO investigate the following lodepng function....
  /*
   unsigned lodepng_decode_memory(unsigned char** out, unsigned* w, unsigned* h,
   const unsigned char* in, size_t insize,
   LodePNGColorType colortype, unsigned bitdepth);
   
   */ 
  
  unsigned width, height;
  
  //TODO arreglar aixo obtenint el tamany del png stream
  width = 55;
  height = 55;
  
  Rcpp::NumericMatrix ionImage(width, height);
  
  return ionImage;
}

//' TODO ideas varies:
//' 
//' - utilitzar un XML de la mateixa forma que fa imzML per tal d'agregar de forma ordenada tot el que necessitis en un binary:
//'   -- Tot el que no pot estar en imzML de espectres: coordenas corregides (actualment pos), normalizacions, etc...
//'   -- Nou ramdisk basat en png: png stream + factor d'escalat de cada m/z channel, pot seguir la seguent trama binaria:
//'         factor_escalat(float 32bit) + png_stream(N bits) [aixo es repetiria per cada pixel]
//'   -- Peak matrix resultant del rMSIproc, seria un camp opcional: m/z vector (centroides), SNR, Intensity i Area (no caldria guardar ni normalizacions i coordenas pq ja hi son!)
//'      
//' - compatibilitat: puc fer que el nou rMSI i rMSIproc detecti automaticament si les dades son nou o antic format i ho carregui?
//' 
//' - comentar coses en el lodepng.h per fer el binary resultant mes petit: el que no facis servir fora! esta documentat en el propi lodepng.h
//' 
//' - imzML processats: no cal crear el imzML continu equivalent, pots mantenir el processat com a "ramdisk" i constriuir els png mitjancant funcio
//'   approx() o interpolacio. Aixi doncs, un paramatre sera l'eix de massa merged previament calculat. 
//'   Es a dir rMSI no convertira imzML processat en continu sino que simplement calculara leix d massa commu i fara servir una interpolacio quan calgui obtenir un espectre concret.
//'   El mateix passara en la part d processat d rMSIproc, en anar carregant espectres s'anira calculant la interpolacio "on the fly"
//'
//' - ajuntar rMSI i rMSIproc: importar tots els metodes de rMSIproc dins del rMSI, es a dir, rMSIproc deixara d existir. Aixi sera mes facil d mantenir.
//'   tb em permetra fer us d funcions d rMSIproc en rMSI (per exemple multi-threading d boost)
//'   oju amb el Namespace d R, crec que rMSI el fa automatic amb roxigen i rMSIproc el fa manual, que vull al final?
//'
//' - creacio d pngstream amb 2 threads. un thread accedeix a disc i l'altre fa encoding. 
//'   Axi sera mes rapid pq mentres un thread fa encoding laltre ja va carregnat el seguent buffer.
//'   
//' - Reconstruccio d'imatges multi-threading. Per construir la imatge dun io i mostrar-la per pantalla cal llegir N png's del pngstream 
//'   (on N correspont al numero d mass chanels o frames de disc), per tan N es correspon amb el valor de tolerancia selecionat per usuari.
//'   La lectura de disc la faria un sol thread k carregaria en RAM N frames. Despres, varis threads es poden repartir la feina de fer el decoding dels N frames.
//'     


/// Pensa en estructura d dades acutal d rMSI i com arribari invocant una sola funcio C++. La estrucura es mante canviant el camp data per les refs a 
/// imzML. Pot ser interessant afegir un camp amb la veriso de format. Axi segons el contingut d versio d format (o si versio no existeix) podre carregar
/// indistintament dades antigue si noves!






//TODO DELETE ME! el k ve tot seguit es un exemple per probar que lodepng funciona i linka be

/*
 3 ways to encode a PNG from RGBA pixel data to a file (and 2 in-memory ways).
 NOTE: this samples overwrite the file or test.png without warning!
 */

//g++ lodepng.cpp examples/example_encode.cpp -I./ -ansi -pedantic -Wall -Wextra -O3

//Example 1
//Encode from raw pixels to disk with a single function call
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  //Encode the image
  unsigned error = lodepng::encode(filename, image, width, height);
  
  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//Example 2
//Encode from raw pixels to an in-memory PNG file first, then write it to disk
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeTwoSteps(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  std::vector<unsigned char> png;
  
  unsigned error = lodepng::encode(png, image, width, height);
  if(!error) lodepng::save_file(png, filename);
  
  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//Example 3
//Save a PNG file to disk using a State, normally needed for more advanced usage.
//The image argument has width * height RGBA pixels or width * height * 4 bytes
void encodeWithState(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
  std::vector<unsigned char> png;
  lodepng::State state; //optionally customize this one
  
  unsigned error = lodepng::encode(png, image, width, height, state);
  if(!error) lodepng::save_file(png, filename);
  
  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
}

//saves image to filename given as argument. Warning, this overwrites the file without warning!
// [[Rcpp::export]]
void testingLodepng() 
{
  //NOTE: this sample will overwrite the file or test.png without warning!
  const char* filename = "test_lodepng.png";
  
  //generate some image
  unsigned width = 512, height = 512;
  std::vector<unsigned char> image;
  image.resize(width * height * 4);
  for(unsigned y = 0; y < height; y++)
    for(unsigned x = 0; x < width; x++) {
      image[4 * width * y + 4 * x + 0] = 255 * !(x & y);
      image[4 * width * y + 4 * x + 1] = x ^ y;
      image[4 * width * y + 4 * x + 2] = x | y;
      image[4 * width * y + 4 * x + 3] = 255;
    }
    
    encodeOneStep(filename, image, width, height);
}


