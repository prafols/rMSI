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
#include <sstream>
#include <stdexcept>
#include <string>
#include <cstdlib>

#include "rMSIXBin.h"
#include "lodepng.h"
#include "pugixml.hpp"
#include "common_methods.h"

using namespace Rcpp;
using namespace pugi;

rMSIXBin::rMSIXBin(String path, String fname)
{
  //Start setting pointers to null to let the destructor to not crash in case of error
  mass = nullptr;
  _rMSIXBin = nullptr;
  
 rMSIObj["data"] = List::create(Named("path") = path, 
                            Named("rMSIXBin") = List::create(Named("file") = fname));
 readXrMSIfile();
}

rMSIXBin::rMSIXBin(List rMSIobject)
{
  rMSIObj = rMSIobject;
  
  irMSIFormatVersion = as<unsigned int>(rMSIObj["rMSI_format_version"]);
  sImgName = as<std::string>(rMSIObj["name"]);
  
  List data = rMSIObj["data"];
  List rMSIXBinData = data["rMSIXBin"];
  
  //Get the UUID's rom the XML part (R  is responsible of verifying the bin part)
  List imzML = data["imzML"];
  sUUID_imzML = as<std::string>(imzML["uuid"]);
  sUUID_rMSIXBin = as<std::string>(rMSIXBinData["uuid"]);
  hexstring2byteuuid(sUUID_imzML, UUID_imzML);
  hexstring2byteuuid(sUUID_rMSIXBin, UUID_rMSIXBin);
  
  NumericVector imgSize = rMSIObj["size"]; 
  img_width = (unsigned int) imgSize["x"];
  img_height = (unsigned int) imgSize["y"];
  
  //File handlers to rMSIXBin
  _rMSIXBin = new rMSIXBin_Handler;
  std::string sFilePath = as<std::string>(data["path"]);
  std::string sFnameImgStream = as<std::string>(rMSIXBinData["file"]);
  _rMSIXBin->XML_file = sFilePath + "/" + sFnameImgStream + ".XrMSI";
  _rMSIXBin->Bin_file = sFilePath + "/" + sFnameImgStream + ".BrMSI";
  
  //Get the mass axis
  NumericVector massAxis =  rMSIObj["mass"];
  massLength = massAxis.length();
  mass = new double[massLength];
  memcpy(mass, massAxis.begin(), sizeof(double)*massLength);
  
  //Get the pixel resolution
  pixel_size_um = as<double>(rMSIObj["pixel_size_um"]);
  
  //Get corrected pixels coords
  NumericMatrix XYCoords = rMSIObj["pos"];
  _rMSIXBin->numOfPixels = XYCoords.nrow();
  _rMSIXBin->iX = new unsigned int[_rMSIXBin->numOfPixels]; 
  _rMSIXBin->iY = new unsigned int[_rMSIXBin->numOfPixels]; 
  for(unsigned int i = 0; i < _rMSIXBin->numOfPixels; i++)
  {
    //Copy the R coords and change the index to make it to start by 0 instead than 1
    _rMSIXBin->iX[i] = XYCoords(i,0) - 1;
    _rMSIXBin->iY[i] = XYCoords(i,1) - 1;
  }
  
  //Get the imgStream from the rMSIObject
  _rMSIXBin->iByteLen = new unsigned long[massLength]; //using long instead of int to allow extra room to store all offsets!
  _rMSIXBin->iByteOffset = new unsigned long[massLength]; //using long instead of int to allow extra room to store all offsets!

  NumericVector RByteLengths = (as<List>((as<List>((as<List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteLength"];
  NumericVector RByteOffsets = (as<List>((as<List>((as<List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteOffset"];
  for(int i=0; i < massLength; i++)
  {
    _rMSIXBin->iByteOffset[i]  = RByteOffsets[i];
    _rMSIXBin->iByteLen[i] = RByteLengths[i];
  }
}

rMSIXBin::~rMSIXBin()
{
  if( mass != nullptr ) //Avoid crashing in case mass ptr was not set
  {
    delete[] mass;
  }
  
  if( _rMSIXBin != nullptr ) //Avoid crashing in case _rMSIXBIN ptr was not set
  {
    delete[] _rMSIXBin->iByteLen;
    delete[] _rMSIXBin->iByteOffset;
    delete[] _rMSIXBin->iX;
    delete[] _rMSIXBin->iY;
  }
  delete _rMSIXBin;
}

List rMSIXBin::get_rMSIObj()
{
  return rMSIObj; 
}

void rMSIXBin::CreateImgStream()
{
  List data;
  List imzML;
  DataFrame imgStream; 
  DataFrame imzMLrun;
  
  NumericVector imzML_mzLength;
  NumericVector imzML_mzOffsets;
  NumericVector imzML_intLength;
  NumericVector imzML_intOffsets;
  
  data=rMSIObj["data"];
  imzML = data["imzML"];
  
  //Get imzML uuid from the XML part (R  is responsible of verifying the bin part)
  hexstring2byteuuid(as<std::string>(imzML["uuid"]), UUID_imzML);
  
  std::string sFilePath = as<std::string>(data["path"]);
  std::string sFnameImzML = as<std::string>(imzML["file"]);
  sFnameImzML= sFilePath + "/" + sFnameImzML + ".ibd";
  
  imzMLrun = imzML["run"];
  imzML_mzLength = imzMLrun["mzLength"];
  imzML_mzOffsets = imzMLrun["mzOffset"];
  imzML_intLength = imzMLrun["intLength"];
  imzML_intOffsets = imzMLrun["intOffset"];

  imzMLDataType mzDataType;
  if( as<std::string>(imzML["mz_dataType"]) == "float" )
  {
    mzDataType = imzMLDataType::float32;
  }
  else if(as<std::string>(imzML["mz_dataType"]) == "double" )
  {
    mzDataType = imzMLDataType::float64;
  }
  else if(as<std::string>(imzML["mz_dataType"]) == "int" )
  {
    mzDataType = imzMLDataType::int32;
  }
  else if(as<std::string>(imzML["mz_dataType"]) == "long" )
  {
    mzDataType = imzMLDataType::int64;
  }
  else
  {
    throw std::runtime_error("ERROR: CreateImgStream() invalid m/z datatype.\n");
  }
  
  imzMLDataType intDataType;
  if( as<std::string>(imzML["int_dataType"]) == "float" )
  {
    intDataType = imzMLDataType::float32;
  }
  else if(as<std::string>(imzML["int_dataType"]) == "double" )
  {
    intDataType = imzMLDataType::float64;
  }
  else if(as<std::string>(imzML["int_dataType"]) == "int" )
  {
    intDataType = imzMLDataType::int32;
  }
  else if(as<std::string>(imzML["int_dataType"]) == "long" )
  {
    intDataType = imzMLDataType::int64;
  }
  else
  {
    throw std::runtime_error("ERROR: CreateImgStream() invalid intensity datatype.\n");
  }
  
  //Create and init the imzML reader
  ImzMLBinRead* imzMLReader;
  try
  {
    imzMLReader = new ImzMLBinRead(sFnameImzML.c_str(), _rMSIXBin->numOfPixels, mzDataType, intDataType, as<bool>(imzML["continuous_mode"])); 
    imzMLReader->set_mzLength(&imzML_mzLength);  
    imzMLReader->set_mzOffset(&imzML_mzOffsets);
    imzMLReader->set_intLength(&imzML_intLength);
    imzMLReader->set_intOffset(&imzML_intOffsets);
  }
  catch(std::runtime_error &e)
  {
    delete imzMLReader;
    stop(e.what());
  }
  
  //Create the binary file (.BrMSI) any previous file will be deleted.
  std::ofstream fBrMSI;
  fBrMSI.open (_rMSIXBin->Bin_file, std::ios::out | std::ios::trunc | std::ios::binary);
  if(!fBrMSI.is_open())
  {
    throw std::runtime_error("Error: rMSIXBin could not open the BrMSI file.\n");
  }
  char UUIDbuffer[16];
  imzMLReader->readUUID(UUIDbuffer);
  fBrMSI.write(UUID_imzML, 16); //Write imzML UUID;
  fBrMSI.write(UUID_rMSIXBin, 16); //Write rMSIXBin UUID;
  
  fBrMSI.write((const char*)mass, sizeof(double) * massLength);
  
  //TODO escriure espectre mig al rMSIXbin (pots omplir amb -1 de moment i ja ho implementaras!)
  fBrMSI.write((const char*)mass, sizeof(double) * massLength); //TODO de momen faig una copia d l'eix de massa i a correr!

  fBrMSI.close();
  if(fBrMSI.fail() || fBrMSI.bad())
  {
    throw std::runtime_error("FATAL ERROR: ImzMLBin got fail or bad bit condition reading the imzML ibd file.\n"); 
  }
  
  //Loop to create each ion image
  /* iIonImgCount calculation
   *  
   *  bytesPerIonImg = img_width * img_height * ENCODING_BITS/8 + 4 (32bits float scalingFactor)
   *  iIonImgCount = IONIMG_BUFFER_MB * 1024 * 1024 / bytesPerIonImg
   */
  unsigned int iIonImgCount = (unsigned int)(  ((double)(IONIMG_BUFFER_MB * 1024 * 1024)) / ((double)( img_width *img_height * ENCODING_BITS/8 + 4 )) );
  unsigned int iRemainingIons = massLength;

  try
  {
    unsigned int iIon = 0;
    Rcout << "Encoding m/z:     ";
    while( true )
    {
      iIonImgCount = iIonImgCount <  iRemainingIons ? iIonImgCount :  iRemainingIons;
      encodeMultipleIonImage2ImgStream(imzMLReader, iIon, iIonImgCount);
      iIon += iIonImgCount;
      iRemainingIons = iRemainingIons - iIonImgCount;
      if(iRemainingIons == 0)
      {
        break;
      }
    }
    Rcout << "\n";
  }
  catch(std::runtime_error &e)
  {
    Rcout << "Encoder Error, stopped\n";
    delete imzMLReader;
    stop(e.what());
  }


  //TODO: pensar funcionament de la part XML del rMSIXBin, crec que el millor es el seguent:
   //1. El constructor parseja el fixer XML original (si existeix i guarda les dades)
   //2. Es modifiquen les dades correponents al ImgStream
   //3. Es sobrescriu el XML sencer amb les dades originals+les noves dades del ImgStream

  delete imzMLReader;
  
  //Copy the _rMSIXBin C style offset to the R rMSIObj
  copyimgStream2rMSIObj(); 
  
  //Write the XML file
  if(!writeXrMSIfile())
  {
    Rcout << ".XrMSI write error, stopped\n";
  }
  
}

void rMSIXBin::encodeMultipleIonImage2ImgStream(ImzMLBinRead* imzMLHandler, unsigned int ionIndex, unsigned int ionCount)
{
  if(imzMLHandler->get_continuous())
  {
    encodeMultipleIonImage2ImgStream_continuous(imzMLHandler, ionIndex, ionCount);
  }
  else
  {
    encodeMultipleIonImage2ImgStream_processed(imzMLHandler, ionIndex, ionCount);
  }
}

//ionIndex: the ion index at which the partial encoding process is started.
//ionCount: the number of ion images to encode at current encoding exectuion.
void rMSIXBin::encodeMultipleIonImage2ImgStream_continuous(ImzMLBinRead* imzMLHandler, unsigned int ionIndex, unsigned int ionCount)
{
  std::vector<double> buffer(ionCount*_rMSIXBin->numOfPixels); 
  for(int i=0; i < _rMSIXBin->numOfPixels; i++)
  {
    imzMLHandler->readIntData(imzMLHandler->get_intOffset(i) + ionIndex*imzMLHandler->get_intEncodingBytes(), ionCount, buffer.data() + (i*ionCount) );
  }
  
  //Prepare the image buffer to encode
  //Init with zeros, observe that non-existing MSI pixels will be zero for all spectra, so there is no need to initialize zeros each time
  std::vector<imgstreamencoding_type> image(img_width * img_height, 0);
  double max;
  float scaling;
  
  std::ofstream fBrMSI;
  fBrMSI.open (_rMSIXBin->Bin_file, std::ios::out | std::ios::app | std::ios::binary);
  if(!fBrMSI.is_open())
  {
    throw std::runtime_error("Error: rMSIXBin could not open the BrMSI file.\n");
  }
  
  for(unsigned int i = 0; i < ionCount; i++)
  { 
    coutEncodingPersentage(ionIndex + i);

    //Calculate the scaling factors
    max = 0.0;
    for(int j=0; j < _rMSIXBin->numOfPixels; j++)
    {
      max = buffer[j*ionCount + i] > max ? buffer[j*ionCount + i] : max;
    }
    
    scaling = (float) max;
    
    for(int j=0; j < _rMSIXBin->numOfPixels; j++)
    {
      image[ _rMSIXBin->iX[j]  + img_width*_rMSIXBin->iY[j] ] = ENCODING_BIT_MASK & (imgstreamencoding_type)(ENCODER_RANGE*((buffer[j*ionCount + i]/max))); //apply scalling and adjust dynamic range
    }
    
    //Encode the current ion image
    std::vector<unsigned char> png_stream;
    unsigned encode_error = lodepng::encode(png_stream, 
                    (const unsigned char*) image.data(), img_width, img_height,
                    LodePNGColorType::LCT_GREY, ENCODING_BITS);
        
    if(encode_error)
    {
      std::stringstream ss; 
      ss << "Error: rMSIXBin png encoding excepion: " << lodepng_error_text(encode_error) << "\n";
      fBrMSI.close();
      throw std::runtime_error(ss.str());
    }
    
    //Save the current image to the imgStream on hdd
    fBrMSI.write((const char*)(&scaling), sizeof(float));  
    fBrMSI.write((const char*)png_stream.data(), png_stream.size());

    //Store offsets info
    _rMSIXBin->iByteLen[i + ionIndex] = sizeof(float) + png_stream.size(); //The encoded bytes are 1) the scaling in a float and 2) the bytes in the png
    if((i + ionIndex) == 0)
    {
      //Special case, the first offset is being writen
      //The first ion image in imgStream will be located at iByteOffset[0] positon of the .BrMSI file.
      //So, 16 bytes for each UUID and massLength bytes for the mass axis and the average spectrum.
      _rMSIXBin->iByteOffset[0] = 16 + 16 + sizeof(double) * massLength + sizeof(double) * massLength;
    }
    else
    {
      
      _rMSIXBin->iByteOffset[i + ionIndex] = _rMSIXBin->iByteOffset[i + ionIndex - 1] + _rMSIXBin->iByteLen[i + ionIndex - 1]; 
    }
    
    if(fBrMSI.fail() || fBrMSI.bad())
    {
      fBrMSI.close();
      throw std::runtime_error("FATAL ERROR: ImzMLBin got fail or bad bit condition reading the imzML ibd file.\n"); 
    }

  }
  fBrMSI.close();
}

void rMSIXBin::encodeMultipleIonImage2ImgStream_processed(ImzMLBinRead* imzMLHandler, unsigned int ionIndex, unsigned int ionCount)
{
  //TODO implement the imzML processed mode methods
  throw std::runtime_error("TODO: The imzML processed mode is not implemented yet, sorry.");
}

//Display the current encoded persentage on console
void rMSIXBin::coutEncodingPersentage(unsigned int ionIndex)
{
  double pp = (double)ionIndex*100.0/((double)massLength - 1.0);
  static double pp_ant = 0.0;
  if( (pp - pp_ant) >= 0.9999999999  || ionIndex == 0  || ionIndex == (massLength-1))
  {
    Rcout << "\b\b\b\b" << std::fixed << std::setprecision(0) << std::setfill(' ') << std::setw(3) << pp << "%";  
    pp_ant = pp;
  }
}

//Write the XML file, any previous .XrMSI file will be deleted
bool rMSIXBin::writeXrMSIfile()
{
  //Reusable pugi variables
  xml_node cvParam;
  
  // empty xml document with custom declaration node
  xml_document doc;
  xml_node decl = doc.prepend_child(node_declaration);
  decl.append_attribute("version") = "1.0";
  decl.append_attribute("encoding") = "UTF-8";
  decl.append_attribute("standalone") = "no";
  
  //mzML top level node
  xml_node node_XrMSI = doc.append_child("XrMSI");
  node_XrMSI.append_attribute("version") =  "1.1";
  node_XrMSI.append_attribute("xmlns") = "http://psi.hupo.org/ms/mzml";
  node_XrMSI.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
  
  //cvList top node
  xml_node node_cvList = node_XrMSI.append_child("cvList");
  node_cvList.append_attribute("count") = "3";
  xml_node node_cv = node_cvList.append_child("cv");
  node_cv.append_attribute("id") = "MS";
  node_cv.append_attribute("fullName") = "Proteomics Standards Initiative Mass Spectrometry Ontology";
  node_cv.append_attribute("version") = "1.3.1";
  node_cv.append_attribute("URI") = "http://psidev.info/ms/mzML/psi-ms.obo";
  node_cv = node_cvList.append_child("cv");
  node_cv.append_attribute("id") = "UO";
  node_cv.append_attribute("fullName") = "Unit Ontology";
  node_cv.append_attribute("version") = "1.15";
  node_cv.append_attribute("URI") = "http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo";
  node_cv = node_cvList.append_child("cv");
  node_cv.append_attribute("id") = "IMS";
  node_cv.append_attribute("fullName") = "Imaging MS Ontology";
  node_cv.append_attribute("version") = "0.9.1";
  node_cv.append_attribute("URI") = "http://www.maldi-msi.org/download/imzml/imagingMS.obo";
  
  node_cv = node_cvList.append_child("cv");
  node_cv.append_attribute("id") = "rMSI";
  node_cv.append_attribute("imgName") = sImgName.c_str();
  node_cv.append_attribute("version") =  irMSIFormatVersion;
  node_cv.append_attribute("URI") = "http://github.com/prafols/rMSI";
  
  //fileDescription node
  xml_node node_fdesc = node_XrMSI.append_child("fileDescription");
  
  //fileContent
  xml_node fileContent = node_fdesc.append_child("fileContent");
  
  //Append imzML file name
  List imzML = (as<List>(rMSIObj["data"]))["imzML"];
  cvParam = fileContent.append_child("cvParam");
  cvParam.append_attribute("accession") = "rMSI:1000000";
  cvParam.append_attribute("cvRef") = "rMSI";
  cvParam.append_attribute("name") = "imzML filename";
  cvParam.append_attribute("value") =  (as<std::string>(imzML["file"])).c_str();
  
  //Append imzML UUID
  cvParam = fileContent.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000080";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "universally unique identifier";
  std::string sUUIDparsed = "{";
  sUUIDparsed += sUUID_imzML.substr(0, 8);
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_imzML.substr(8, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_imzML.substr(12, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_imzML.substr(16, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_imzML.substr(20, 12); 
  sUUIDparsed += "}";
  cvParam.append_attribute("value") = sUUIDparsed.c_str();
  
  //Append rMSI UUID
  cvParam = fileContent.append_child("cvParam");
  cvParam.append_attribute("accession") = "rMSI:1000080";
  cvParam.append_attribute("cvRef") = "rMSI";
  cvParam.append_attribute("name") = "rMSI universally unique identifier";
  sUUIDparsed = "{";
  sUUIDparsed += sUUID_rMSIXBin.substr(0, 8);
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_rMSIXBin.substr(8, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_rMSIXBin.substr(12, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_rMSIXBin.substr(16, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUID_rMSIXBin.substr(20, 12); 
  sUUIDparsed += "}";
  cvParam.append_attribute("value") = sUUIDparsed.c_str();
  
  //contact info.
  xml_node node_contact = node_fdesc.append_child("contact");
  cvParam = node_contact.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000586";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "contact name";
  cvParam.append_attribute("value") = "Pere Rafols";
  
  cvParam = node_contact.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000590";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "contact affiliation";
  cvParam.append_attribute("value") = "Universitat Rovira i Virgili";
  
  cvParam = node_contact.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000589";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "contact email";
  cvParam.append_attribute("value") = "pere.rafols@urv.cat";
  
  //scanSettingsList
  xml_node node_scanSetLst = node_XrMSI.append_child("scanSettingsList");
  node_scanSetLst.append_attribute("count") = "1";
  xml_node node_scanSet = node_scanSetLst.append_child("scanSettings");
  node_scanSet.append_attribute("id") = "scanSettings0";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "rMSI:1000010";
  cvParam.append_attribute("cvRef") = "rMSI";
  cvParam.append_attribute("name") = "max count of m/z channels";
  cvParam.append_attribute("value") = massLength;
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000042";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "max count of pixels x";
  cvParam.append_attribute("value") = img_width;
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000043";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "max count of pixels y";
  cvParam.append_attribute("value") = img_height;
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000046";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "pixel size";
  cvParam.append_attribute("value") = pow(pixel_size_um, 2.0);
  
  //Run data spectra list
  xml_node node_spectrum; //Reusable spectrum node
  xml_node node_run = node_XrMSI.append_child("run");
  xml_node node_spectrumLst = node_run.append_child("spectrumList");
  node_spectrumLst.append_attribute("count") = _rMSIXBin->numOfPixels;

  //Get the motor coordinates
  NumericMatrix XYCoordsMotors = rMSIObj["posMotors"];

  for( int i = 0; i <  _rMSIXBin->numOfPixels; i++)
  {
    //Rcout<<"Parsing pixel "<<i <<" of "<< total_num_pixels << "\n";
    node_spectrum = node_spectrumLst.append_child("spectrum");
    node_spectrum.append_attribute("id") = i;
    
    cvParam = node_spectrum.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000050";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "position x";
    cvParam.append_attribute("value") = (int)XYCoordsMotors(i,0);
    
    cvParam = node_spectrum.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000051";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "position y";
    cvParam.append_attribute("value") = (int)XYCoordsMotors(i,1);;
    
    cvParam = node_spectrum.append_child("cvParam");
    cvParam.append_attribute("accession") = "rMSI:1000050";
    cvParam.append_attribute("cvRef") = "rMSI";
    cvParam.append_attribute("name") = "corrected position x";
    cvParam.append_attribute("value") = _rMSIXBin->iX[i];
    
    cvParam = node_spectrum.append_child("cvParam");
    cvParam.append_attribute("accession") = "rMSI:1000051";
    cvParam.append_attribute("cvRef") = "rMSI";
    cvParam.append_attribute("name") = "corrected position y";
    cvParam.append_attribute("value") = _rMSIXBin->iY[i];
  }
  
  //Run data: imgStream
  xml_node node_ionimg; //Reusable ionimg node
  xml_node node_imgStream = node_run.append_child("imgStreamList");
  node_imgStream.append_attribute("count") = massLength;
  for( int i = 0; i < massLength; i++)
  {
    node_ionimg = node_imgStream.append_child("ionImage");
    node_ionimg.append_attribute("id") = i;
    
    cvParam = node_ionimg.append_child("cvParam");
    cvParam.append_attribute("accession") = "rMSI:1000060";
    cvParam.append_attribute("cvRef") = "rMSI";
    cvParam.append_attribute("name") = "image ion byte count";
    cvParam.append_attribute("value") = _rMSIXBin->iByteLen[i];
    
    cvParam = node_ionimg.append_child("cvParam");
    cvParam.append_attribute("accession") = "rMSI:1000061";
    cvParam.append_attribute("cvRef") = "rMSI";
    cvParam.append_attribute("name") = "image ion byte offset";
    cvParam.append_attribute("value") = _rMSIXBin->iByteOffset[i];
  }
  
  //Run data: Normalizations
  //TODO
    
  // save document to file
  return(doc.save_file(_rMSIXBin->XML_file.c_str(), "\t", format_default, encoding_utf8 ) );
}

//Read the XML file and fill all the rMSIobject data
void rMSIXBin::readXrMSIfile()
{
  String xml_path = as<String>((as<List>(rMSIObj["data"]))["path"]);
    xml_path += "/";
    xml_path += as<String>((as<List>((as<List>(rMSIObj["data"]))["rMSIXBin"]))["file"]);
    xml_path += ".XrMSI";
  
  
  xml_document doc;
  xml_parse_result  result = doc.load_file(xml_path.get_cstring());
  
  if (!result)
  {
    String err = "ERROR: XML [";
    err += xml_path.get_cstring();
    err += "] parsed with errors, attr value: [";
    err += doc.child("node").attribute("attr").value();
    err += "]\nDescription:";
    err += result.description();
    throw std::runtime_error(err.get_cstring());
  }
  
  xml_node XrMSI = doc.child("XrMSI");
  if( XrMSI == NULL )
  {
    throw std::runtime_error("XML parse error: no XrMSI node found");
  }
  
  //Generic fields
  std::string accession; //Accession to XML nodes...
  
  //Fiels to retrive
  std::string strImzML_filename;
  
  //Parse cvList node 
  xml_node cvList = XrMSI.child("cvList");
  if( cvList == NULL )
  {
    throw std::runtime_error("XML parse error: no cvList node found");
  } 
  
  for (xml_node cv = cvList.child("cv"); cv; cv = cv.next_sibling("cv"))
  {
    accession = cv.attribute("id").value();
    if(accession == "rMSI")
    {
      sImgName = cv.attribute("imgName").value();
      irMSIFormatVersion = cv.attribute("version").as_uint();
    }
  }
  
  //Parse fileDescription node
  xml_node fileDesc = XrMSI.child("fileDescription");
  if( fileDesc == NULL )
  {
    throw std::runtime_error("XML parse error: no fileDescription node found");
  } 
  
  xml_node fileContent = fileDesc.child("fileContent");
  if( fileContent == NULL )
  {
    throw std::runtime_error("XML parse error: no fileContent node found");
  } 
  
  for (xml_node cvParam = fileContent.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
  {
    accession = cvParam.attribute("accession").value();
    if(accession == "rMSI:1000000")
    {
      //imzML filename
      strImzML_filename = cvParam.attribute("value").value();
    }
    if(accession == "IMS:1000080")
    {
      //imzML uuid
      sUUID_imzML = cvParam.attribute("value").value();
      sUUID_imzML = parse_xml_uuid(sUUID_imzML);
      hexstring2byteuuid(sUUID_imzML, UUID_imzML);
    }
    if(accession == "rMSI:1000080")
    {
      //XrMSI uuid
      sUUID_rMSIXBin = cvParam.attribute("value").value();
      sUUID_rMSIXBin = parse_xml_uuid(sUUID_rMSIXBin);
      hexstring2byteuuid(sUUID_rMSIXBin, UUID_rMSIXBin);
    }
  }
  
  //Parse scanSettingsList node
  xml_node scanSettingsList = XrMSI.child("scanSettingsList");
  if( scanSettingsList == NULL )
  {
    throw std::runtime_error("XML parse error: no scanSettingsList node found");
  }
  
  xml_node scanSettings = scanSettingsList.child("scanSettings");
  if( scanSettings == NULL )
  {
    throw std::runtime_error("XML parse error: no scanSettings node found");
  } 
  
  for (xml_node cvParam = scanSettings.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
  {
    accession = cvParam.attribute("accession").value();
    if(accession == "rMSI:1000010")
    {
      //mass channels
      massLength = cvParam.attribute("value").as_uint();
    }
    if(accession == "IMS:1000042")
    {
      //pixels in X (image width)
      img_width = cvParam.attribute("value").as_uint();
    }
    if(accession == "IMS:1000043")
    {
      //pixels in Y (image height)
      img_height = cvParam.attribute("value").as_uint();
    }
    if(accession == "IMS:1000046")
    {
      //pixel size
      pixel_size_um = sqrt(cvParam.attribute("value").as_double());
    }
  }
  
  //Parse run node
  xml_node run = XrMSI.child("run");
  if( run == NULL )
  {
    throw std::runtime_error("XML parse error: no run node found");
  }
  
  xml_node spectrumList = run.child("spectrumList");
  if( spectrumList == NULL )
  {
    throw std::runtime_error("XML parse error: no spectrumList node found");
  }
  
  xml_node imgStreamList = run.child("imgStreamList");
  if( imgStreamList == NULL )
  {
    throw std::runtime_error("XML parse error: no imgStreamList node found");
  }
  
  if(imgStreamList.attribute("count").as_uint() != massLength)
  {
    throw std::runtime_error("XML parse error: imgStreamList length is different than mass axis length");
  }
  
  _rMSIXBin = new rMSIXBin_Handler;
  std::string sFilePath = as<String>((as<List>(rMSIObj["data"]))["path"]);
  std::string sFnameImgStream =  as<String>((as<List>((as<List>(rMSIObj["data"]))["rMSIXBin"]))["file"]);
  _rMSIXBin->XML_file = sFilePath + "/" + sFnameImgStream + ".XrMSI";
  _rMSIXBin->Bin_file = sFilePath + "/" + sFnameImgStream + ".BrMSI";
  _rMSIXBin->numOfPixels = spectrumList.attribute("count").as_uint();
  _rMSIXBin->iX = new unsigned int[_rMSIXBin->numOfPixels]; 
  _rMSIXBin->iY = new unsigned int[_rMSIXBin->numOfPixels]; 
  _rMSIXBin->iByteLen = new unsigned long[massLength];
  _rMSIXBin->iByteOffset = new unsigned long[massLength];
  unsigned int  id;
  
  //Read the position matrices
  NumericMatrix pos(_rMSIXBin->numOfPixels, 2);
  colnames(pos) = CharacterVector::create("x", "y");
  NumericMatrix posMotors(_rMSIXBin->numOfPixels, 2);
  colnames(posMotors) = CharacterVector::create("x", "y");
  for (xml_node spectrum = spectrumList.child("spectrum"); spectrum; spectrum = spectrum.next_sibling("spectrum"))
  {
    id = spectrum.attribute("id").as_uint();
    for (xml_node cvParam = spectrum.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
    {
      accession = cvParam.attribute("accession").value();
      if(accession == "IMS:1000050")
      {
        //X Motor pos
        posMotors(id, 0) = cvParam.attribute("value").as_double();
      }
      if(accession == "IMS:1000051")
      {
        //Y Motor pos
        posMotors(id, 1) = cvParam.attribute("value").as_double();
      }
      if(accession == "rMSI:1000050")
      {
        //X pos
        _rMSIXBin->iX[id] = cvParam.attribute("value").as_uint();
        pos(id, 0) = _rMSIXBin->iX[id] + 1;  //+1 to get it in R indexing
      }
      if(accession == "rMSI:1000051")
      {
        //Y pos
        _rMSIXBin->iY[id] = cvParam.attribute("value").as_uint();
        pos(id, 1) = _rMSIXBin->iY[id] + 1; //+1 to get it in R indexing
      }
    }
  }
  
  //Read the imgStream offsets
  for (xml_node ionImage = imgStreamList.child("ionImage"); ionImage; ionImage = ionImage.next_sibling("ionImage"))
  {
    id = ionImage.attribute("id").as_uint();
    for (xml_node cvParam = ionImage.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
    {
      accession = cvParam.attribute("accession").value();
      if(accession == "rMSI:1000060")
      {
        //image ion byte count
        _rMSIXBin->iByteLen[id] = cvParam.attribute("value").as_ullong();
      }
      if(accession == "rMSI:1000061")
      {
        //image ion byte offset
        _rMSIXBin->iByteOffset[id] = cvParam.attribute("value").as_ullong();
      }
    }
  }

  //Fill rMSIObj info
  NumericVector base(massLength); //Empty base spectrum
  rMSIObj.push_front(base, "base");
  
  NumericVector mean(massLength); //Empty mean spectrum
  rMSIObj.push_front(mean, "mean");
  
  rMSIObj.push_front(pixel_size_um, "pixel_size_um"); 
  
  rMSIObj.push_front(posMotors, "posMotors");
  rMSIObj.push_front(pos, "pos"); 
  
  NumericVector size;
  size.push_back(img_width, "x");
  size.push_back(img_height, "y");
  rMSIObj.push_front(size, "size");
  
  NumericVector mass(massLength); //Empty mass
  rMSIObj.push_front(mass, "mass");
  
  rMSIObj.push_front(sImgName, "name");
  
  rMSIObj.push_front(irMSIFormatVersion, "rMSI_format_version");
  
  //Get the data field
  List data_lst = as<List>(rMSIObj["data"]);
  
  //Fill rMSIXBin info
  List rMSIXBIN_lst = as<List>(data_lst["rMSIXBin"]);
  
  //Just start with an empty imgStream, it is copied at the end
  List imgStream_lst = List::create(Named("ByteLength") = NumericVector(),
                                    Named("ByteOffset") = NumericVector());
  imgStream_lst.attr("class") = "imgStream"; //Set class type
  
  rMSIXBIN_lst.push_front(imgStream_lst, "imgStream");
  rMSIXBIN_lst.push_front(sUUID_rMSIXBin, "uuid");
  
  rMSIXBIN_lst.attr("class") = "rMSIXBinData"; //Set class name
  data_lst["rMSIXBin"] = rMSIXBIN_lst; //Overwrite the rMSIXBin
  
  //Fill imzML info
  List imzML_lst = List::create(Named("uuid") = sUUID_imzML,  //TODO not working!
                                Named("file") = strImzML_filename //TODO not working!
                                );
  imzML_lst.attr("class") = "imzMLData"; //Set class name
  
  data_lst.push_back(imzML_lst, "imzML");
  data_lst.attr("class") = "rMSIData"; //Set class name
  rMSIObj["data"] = data_lst; //Overwrite the data
  
  //Copy the _rMSIXBin C style offset to the R rMSIObj
  copyimgStream2rMSIObj(); 
  
  rMSIObj.attr("class") = "rMSIObj"; //Set class name
  
}

//Copy imgStream to the rMSIObject
void rMSIXBin::copyimgStream2rMSIObj()
{
  NumericVector RByteLengths(massLength);
  NumericVector RByteOffsets(massLength);
  for(int i=0; i < massLength; i++)
  {
    RByteOffsets[i] = _rMSIXBin->iByteOffset[i]; 
    RByteLengths[i] = _rMSIXBin->iByteLen[i];
  }
  (as<List>((as<List>((as<List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteLength"] = RByteLengths;
  (as<List>((as<List>((as<List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteOffset"] = RByteOffsets ;
}

//' decodePngStream2IonImages.
//'
//' Obtain a multiple mass channel ion image by decoding the hdd img stream at a specified ionIndex.
//' The MAX operator will be used to merge all ion images in a single image matrix.
//'
//' @param ionIndex the index of ion to extract from the img stream. C style indexing, starting with zero.
//' @param ionCount number of ion image to decode.
//' 
//' @return A NumerixMatrix containing the ion image.
//' 
NumericMatrix rMSIXBin::decodeImgStream2IonImages(unsigned int ionIndex, unsigned int ionCount)
{
  if(ionIndex + ionCount > massLength)
  {
    throw std::runtime_error("ERROR in rMSIXBin::decodeImgStream2IonImages(): ionIndex+ionCount is out of range.\n");
  }
  
  //1- Read the complete stream in a buffer
  //Calculate the total number of bytes to read the length vector
  unsigned long byte_count = 0;
  for(int i=ionIndex; i < (ionIndex+ionCount); i++)
  {
    byte_count += _rMSIXBin->iByteLen[i];
  }
  
  //Too large ion image windows, raise an error
  if( byte_count > IONIMG_BUFFER_MB * 1024 * 1024)
  {
    throw std::runtime_error("ERROR in rMSIXBin::decodeImgStream2IonImages(): number of mass channels too large to load in memory.\n");
  }
  
  //Read tje complete buffer
  char* buffer = new char[byte_count];
  std::ifstream  binFile;
  binFile.open(_rMSIXBin->Bin_file, std::fstream::in | std::ios::binary);
  if(!binFile.is_open())
  {
    delete[] buffer;
    throw std::runtime_error("ERROR: rMSIXBin::decodeImgStream2IonImages could not open the .BrMSI file.\n"); 
  }
  
  binFile.seekg(_rMSIXBin->iByteOffset[ionIndex]);
  if(binFile.eof())
  {
    binFile.close();
    delete[] buffer;
    throw std::runtime_error("ERROR: rMSIXBin::decodeImgStream2IonImages reached EOF seeking the .BrMSI file.\n"); 
  }
  if(binFile.fail() || binFile.bad())
  {
    binFile.close();
    delete[] buffer;
    throw std::runtime_error("FATAL ERROR: rMSIXBin::decodeImgStream2IonImages got fail or bad bit condition seeking the .BrMSI file.\n"); 
  }
  
  binFile.read (buffer, byte_count);
  if(binFile.eof())
  {
    binFile.close();
    delete[] buffer;
    throw std::runtime_error("ERROR: rMSIXBin::decodeImgStream2IonImages reached EOF reading the .BrMSI file.\n"); 
  }
  if(binFile.fail() || binFile.bad())
  {
    binFile.close();
    delete[] buffer;
    throw std::runtime_error("FATAL ERROR:  rMSIXBin::decodeImgStream2IonImages got fail or bad bit condition reading the .BrMSI file.\n"); 
  }
  binFile.close();
  
  //2- Decode the buffer
  NumericMatrix ionImage(img_width, img_height);
  std::vector<unsigned char> raw_image;
  float scaling;
  unsigned int png_width, png_height;
  for(int i=0; i < ionCount; i++)
  {
    //Read the scaling factor
    std::memcpy(&scaling, buffer + _rMSIXBin->iByteOffset[i + ionIndex] - _rMSIXBin->iByteOffset[ionIndex], sizeof(float));
    
    unsigned encode_error = lodepng::decode(raw_image, png_width, png_height,
                    (const unsigned char*)(buffer + _rMSIXBin->iByteOffset[i + ionIndex] - _rMSIXBin->iByteOffset[ionIndex] + sizeof(float)), _rMSIXBin->iByteLen[i + ionIndex] - sizeof(float),
                    LodePNGColorType::LCT_GREY, ENCODING_BITS);
    
    if(encode_error)
    {
      std::stringstream ss; 
      ss << "Error: rMSIXBin png decoding excepion: " << lodepng_error_text(encode_error) << "\n";
      delete[] buffer;
      throw std::runtime_error(ss.str());
    }
    
    if(png_width != img_width || png_height != img_height)
    {
      delete[] buffer;
      throw std::runtime_error("ERROR:  rMSIXBin::decodeImgStream2IonImages decoded image size is invalid, possible data corruption in .BrMSI file.\n");
    }
    
    imgstreamencoding_type pixel_value_raw; //Current pixel value in raw format
    double pixel_value; //Current pixel value in R format
    unsigned int img_offset; //Offset inside the raw image
    unsigned int img_x = 0; //current x coordinate in the image
    unsigned int img_y = 0; //current y coordinate in the image
    for( int iPixel = 0; iPixel <  img_width * img_height; iPixel++)
    {
      img_offset = img_x  + img_width*img_y;
      std::memcpy(&pixel_value_raw, raw_image.data() + img_offset*sizeof(imgstreamencoding_type), sizeof(imgstreamencoding_type));
      pixel_value = (((double)pixel_value_raw)/ENCODER_RANGE) * (double)scaling; 
      ionImage(img_x,img_y) = pixel_value > ionImage(img_x,img_y) ? pixel_value : ionImage(img_x,img_y);
      
      img_x++;
      if(img_x >= img_width)
      {
        img_x = 0;
        img_y++;
      }
    }
  }

  delete[] buffer;
  return ionImage;
}

//Convert a std::string containing a 16 bytes UUID to a big-endian formate byte stream ready to write it to a binary file
void rMSIXBin::hexstring2byteuuid(std::string hex_str, char* output)
{
  if(hex_str.length() != 32)
  {
    throw std::invalid_argument("ERROR: hexstring2byteuuid() invalid input string length");
  }
  
  char byteStr[2];
  unsigned int iout = 0;
  for(int i=0; i < 32; i+=2)
  {
    memcpy(byteStr, hex_str.c_str() + i, 2);
    output[iout] = (char) std::strtol( byteStr, NULL, 16);
    iout++;
  }
}



//R exported methods

//' Ccreate_rMSIXBinData.
//' 
//' creates new rMSIXBin files (.XrMSI and .BrMSI). Previous files will be deleted.
//'
//' @param rMSIobj: an rMSI object prefilled with a parsed imzML.
//' @return the rMSI object with rMSIXBin inforation completed. 
// [[Rcpp::export]]
List Ccreate_rMSIXBinData(List rMSIobj)
{
  try
  {
    rMSIXBin myXBin(rMSIobj); 
    myXBin.CreateImgStream();
    return myXBin.get_rMSIObj();
  }
  catch(std::runtime_error &e)
  {
    stop(e.what());
  }
  return NULL;
}

//' Cload_rMSIXBinData.
//' 
//' Loads the data from the rMSIXBin files (.XrMSI and .BrMSI).
//' This method is used to load a previously stored rMSIXBin file.
//'
//' @param path: full path to the .XrMSI file.
//' @param fname: file name of the .XrMSI file without the extension.
//' @return the rMSI object with rMSIXBin inforation completed. 
// [[Rcpp::export]]
List Cload_rMSIXBinData(String path, String fname)
{
  try
  {
    rMSIXBin myXBin(path, fname); 
    return myXBin.get_rMSIObj();
  }
  catch(std::runtime_error &e)
  {
    stop(e.what());
  }
  return NULL; 
}

//' Cload_rMSIXBinIonImage.
//' 
//' loads a ion image from the .BrNSI img stream.
//' 
//' @param rMSIobj: an rMSI object prefilled with a parsed imzML.
//' @param ionIndex: the first mass channel at which the image starts.
//' @param ionCount: the numer of mass channels used to construct the ion image (a.k.a. image tolerance window).
//' @return the ion image as a NumericMatrix using max operator with all the ion images of the mass channels. 
// [[Rcpp::export]]
NumericMatrix Cload_rMSIXBinIonImage(List rMSIobj, unsigned int ionIndex, unsigned int ionCount)
{
  //Check if ion indeces are valid
  if(ionIndex < 1)
  {
    throw std::runtime_error("ERROR in rMSIXBin::Cload_rMSIXBinIonImage(): ionIndex below zero.\n");
  }
  
  ionIndex = ionIndex - 1; //ionIndex-1 to convert from R to C indexing
  
  try
  {
    rMSIXBin myXBin(rMSIobj); 
    return myXBin.decodeImgStream2IonImages(ionIndex, ionCount); //ionIndex-1 to convert from R to C indexing
  }
  catch(std::runtime_error &e)
  {
    stop(e.what());
  }
  return NumericMatrix(); //Returning empty matrix in cas of error
}


//' TODO ideas varies:
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

