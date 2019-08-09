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

rMSIXBin::rMSIXBin(Rcpp::List rMSIobject)
{
  rMSIObj = rMSIobject;
  
  irMSIFormatVersion = Rcpp::as<unsigned int>(rMSIObj["rMSI_format_version"]);
  sImgName = Rcpp::as<std::string>(rMSIObj["name"]);
  
  Rcpp::List data = rMSIObj["data"];
  Rcpp::List rMSIXBinData = data["rMSIXBin"];
  
  //Get the UUID's rom the XML part (R  is responsible of verifying the bin part)
  Rcpp::List imzML = data["imzML"];
  sUUID_imzML = Rcpp::as<std::string>(imzML["uuid"]);
  sUUID_rMSIXBin = Rcpp::as<std::string>(rMSIXBinData["uuid"]);
  hexstring2byteuuid(sUUID_imzML, UUID_imzML);
  hexstring2byteuuid(sUUID_rMSIXBin, UUID_rMSIXBin);
  
  Rcpp::NumericVector imgSize = rMSIObj["size"]; 
  img_width = (unsigned int) imgSize["x"];
  img_height = (unsigned int) imgSize["y"];
  
  //File handlers to rMSIXBin
  _rMSIXBin = new rMSIXBin_Handler;
  std::string sFilePath = Rcpp::as<std::string>(data["path"]);
  std::string sFnameImgStream = Rcpp::as<std::string>(rMSIXBinData["file"]);
  _rMSIXBin->XML_file = sFilePath + "/" + sFnameImgStream + ".XrMSI";
  _rMSIXBin->Bin_file = sFilePath + "/" + sFnameImgStream + ".BrMSI";
  
  //Get the mass axis
  Rcpp::NumericVector massAxis =  rMSIObj["mass"];
  massLength = massAxis.length();
  mass = new double[massLength];
  memcpy(mass, massAxis.begin(), sizeof(double)*massLength);
  
  //Get corrected pixels coords
  Rcpp::NumericMatrix XYCoords = rMSIObj["pos"];
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
  _rMSIXBin->fScaling = new float[massLength];
  _rMSIXBin->iByteLen = new unsigned long[massLength]; //using long instead of int to allow extra room to store all offsets!
  _rMSIXBin->iByteOffset = new unsigned long[massLength]; //using long instead of int to allow extra room to store all offsets!

  Rcpp::IntegerVector RScalings = (Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["Scaling"] ;
  Rcpp::NumericVector RByteLengths = (Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteLength"];
  Rcpp::NumericVector RByteOffsets = (Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteOffset"];
  for(int i=0; i < massLength; i++)
  {
    _rMSIXBin->fScaling[i] = (float)RScalings[i];
    _rMSIXBin->iByteOffset[i]  = RByteOffsets[i];
    _rMSIXBin->iByteLen[i] = RByteLengths[i];
  }
}

rMSIXBin::~rMSIXBin()
{
  delete[] mass;
  delete[] _rMSIXBin->fScaling;
  delete[] _rMSIXBin->iByteLen;
  delete[] _rMSIXBin->iByteOffset;
  delete[] _rMSIXBin->iX;
  delete[] _rMSIXBin->iY;
  delete _rMSIXBin;
}

Rcpp::List rMSIXBin::get_rMSIObj()
{
  return rMSIObj; 
}

void rMSIXBin::CreateImgStream()
{
  Rcpp::List data;
  Rcpp::List imzML;
  Rcpp::DataFrame imgStream; 
  Rcpp::DataFrame imzMLrun;
  
  Rcpp::NumericVector imzML_mzLength;
  Rcpp::NumericVector imzML_mzOffsets;
  Rcpp::NumericVector imzML_intLength;
  Rcpp::NumericVector imzML_intOffsets;
  
  data=rMSIObj["data"];
  imzML = data["imzML"];
  
  //Get imzML uuid from the XML part (R  is responsible of verifying the bin part)
  hexstring2byteuuid(Rcpp::as<std::string>(imzML["uuid"]), UUID_imzML);
  
  std::string sFilePath = Rcpp::as<std::string>(data["path"]);
  std::string sFnameImzML = Rcpp::as<std::string>(imzML["file"]);
  sFnameImzML= sFilePath + "/" + sFnameImzML + ".ibd";
  
  imzMLrun = imzML["run"];
  imzML_mzLength = imzMLrun["mzLength"];
  imzML_mzOffsets = imzMLrun["mzOffset"];
  imzML_intLength = imzMLrun["intLength"];
  imzML_intOffsets = imzMLrun["intOffset"];

  imzMLDataType mzDataType;
  if( Rcpp::as<std::string>(imzML["mz_dataType"]) == "float" )
  {
    mzDataType = imzMLDataType::float32;
  }
  else if(Rcpp::as<std::string>(imzML["mz_dataType"]) == "double" )
  {
    mzDataType = imzMLDataType::float64;
  }
  else if(Rcpp::as<std::string>(imzML["mz_dataType"]) == "int" )
  {
    mzDataType = imzMLDataType::int32;
  }
  else if(Rcpp::as<std::string>(imzML["mz_dataType"]) == "long" )
  {
    mzDataType = imzMLDataType::int64;
  }
  else
  {
    throw std::runtime_error("ERROR: CreateImgStream() invalid m/z datatype.\n");
  }
  
  imzMLDataType intDataType;
  if( Rcpp::as<std::string>(imzML["int_dataType"]) == "float" )
  {
    intDataType = imzMLDataType::float32;
  }
  else if(Rcpp::as<std::string>(imzML["int_dataType"]) == "double" )
  {
    intDataType = imzMLDataType::float64;
  }
  else if(Rcpp::as<std::string>(imzML["int_dataType"]) == "int" )
  {
    intDataType = imzMLDataType::int32;
  }
  else if(Rcpp::as<std::string>(imzML["int_dataType"]) == "long" )
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
    imzMLReader = new ImzMLBinRead(sFnameImzML.c_str(), _rMSIXBin->numOfPixels, mzDataType, intDataType, Rcpp::as<bool>(imzML["continuous_mode"])); 
    imzMLReader->set_mzLength(&imzML_mzLength);  
    imzMLReader->set_mzOffset(&imzML_mzOffsets);
    imzMLReader->set_intLength(&imzML_intLength);
    imzMLReader->set_intOffset(&imzML_intOffsets);
  }
  catch(std::runtime_error &e)
  {
    delete imzMLReader;
    Rcpp::stop(e.what());
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
    Rcpp::Rcout << "Encoding m/z:     ";
    while( true )
    {
      //Rcpp::Rcout << "Encoding m/z " << iIon << " of " << massLength << "\n"; //TODO millorar aixo... es mooolt lleig... potser posaro dins la funcio d encoding
      iIonImgCount = iIonImgCount <  iRemainingIons ? iIonImgCount :  iRemainingIons;
      encodeMultipleIonImage2ImgStream(imzMLReader, iIon, iIonImgCount);
      iIon += iIonImgCount;
      iRemainingIons = iRemainingIons - iIonImgCount;
      if(iRemainingIons == 0)
      {
        break;
      }
    }
    Rcpp::Rcout << "\n";
  }
  catch(std::runtime_error &e)
  {
    Rcpp::Rcout << "Encoder Error, stopped\n";
    delete imzMLReader;
    Rcpp::stop(e.what());
  }


  //TODO: pensar funcionament de la part XML del rMSIXBin, crec que el millor es el seguent:
   //1. El constructor parseja el fixer XML original (si existeix i guarda les dades)
   //2. Es modifiquen les dades correponents al ImgStream
   //3. Es sobrescriu el XML sencer amb les dades originals+les noves dades del ImgStream

  delete imzMLReader;
  
  //Copy imgStream to the rMSIObject
  Rcpp::IntegerVector RScalings(massLength);
  Rcpp::NumericVector RByteLengths(massLength);
  Rcpp::NumericVector RByteOffsets(massLength);
  for(int i=0; i < massLength; i++)
  {
    RScalings[i] =      (double)_rMSIXBin->fScaling[i]; 
    RByteOffsets[i] = _rMSIXBin->iByteOffset[i]; 
    RByteLengths[i] = _rMSIXBin->iByteLen[i];
  }
  (Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteLength"] = RByteLengths;
  (Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["ByteOffset"] = RByteOffsets ;
  (Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>((Rcpp::as<Rcpp::List>(rMSIObj["data"]))["rMSIXBin"]))["imgStream"]))["Scaling"] = RScalings;
  
  //Write the XML file
  if(!writeXrMSIfile())
  {
    Rcpp::Rcout << ".XrMSI write error, stopped\n";
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
    imzMLHandler->readIntData(imzMLHandler->get_intOffset(i) + ionIndex, ionCount, buffer.data() + (i*ionCount) );
  }
  
  //Prepare the image buffer to encode
  //Init with zeros, observe that non-existing MSI pixels will be zero for all spectra, so there is no need to initialize zeros each time
  std::vector<imgstreamencoding_type> image(img_width * img_height, 0);
  double max;
  
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
    
    _rMSIXBin->fScaling[ionIndex + i] = (float) max; //Store the scaling
    
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
    fBrMSI.write((const char*)((_rMSIXBin->fScaling)+ionIndex + i), sizeof(float));  
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
    Rcpp::Rcout << "\b\b\b\b" << std::fixed << std::setprecision(0) << std::setfill(' ') << std::setw(3) << pp << "%";  
    pp_ant = pp;
  }
}

//Write the XML file, any previous .XrMSI file will be deleted
bool rMSIXBin::writeXrMSIfile()
{
  //Reusable pugi variables
  pugi::xml_node cvParam;
  
  // empty xml document with custom declaration node
  pugi::xml_document doc;
  pugi::xml_node decl = doc.prepend_child(pugi::node_declaration);
  decl.append_attribute("version") = "1.0";
  decl.append_attribute("encoding") = "UTF-8";
  decl.append_attribute("standalone") = "no";
  
  //mzML top level node
  pugi::xml_node node_XrMSI = doc.append_child("XrMSI");
  node_XrMSI.append_attribute("version") =  "1.1";
  node_XrMSI.append_attribute("xmlns") = "http://psi.hupo.org/ms/mzml";
  node_XrMSI.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
  
  //cvList top node
  pugi::xml_node node_cvList = node_XrMSI.append_child("cvList");
  node_cvList.append_attribute("count") = "3";
  pugi::xml_node node_cv = node_cvList.append_child("cv");
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
  pugi::xml_node node_fdesc = node_XrMSI.append_child("fileDescription");
  
  //fileContent
  pugi::xml_node fileContent = node_fdesc.append_child("fileContent");
  
  //Append imzML file name
  Rcpp::List imzML = (Rcpp::as<Rcpp::List>(rMSIObj["data"]))["imzML"];
  cvParam = fileContent.append_child("cvParam");
  cvParam.append_attribute("accession") = "rMSI:1000000";
  cvParam.append_attribute("cvRef") = "rMSI";
  cvParam.append_attribute("name") = "imzML filename";
  cvParam.append_attribute("value") =  (Rcpp::as<std::string>(imzML["file"])).c_str();
  
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
  pugi::xml_node node_contact = node_fdesc.append_child("contact");
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
  pugi::xml_node node_scanSetLst = node_XrMSI.append_child("scanSettingsList");
  node_scanSetLst.append_attribute("count") = "1";
  pugi::xml_node node_scanSet = node_scanSetLst.append_child("scanSettings");
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
  cvParam.append_attribute("value") = pow(Rcpp::as<double>(rMSIObj["pixel_size_um"]), 2.0);
  
  //Run data spectra list
  pugi::xml_node node_spectrum; //Reusable spectrum node
  pugi::xml_node node_run = node_XrMSI.append_child("run");
  pugi::xml_node node_spectrumLst = node_run.append_child("spectrumList");
  node_spectrumLst.append_attribute("count") = _rMSIXBin->numOfPixels;

  //Get the motor coordinates
  Rcpp::NumericMatrix XYCoordsMotors = rMSIObj["posMotors"];

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
  pugi::xml_node node_ionimg; //Reusable ionimg node
  pugi::xml_node node_imgStream = node_run.append_child("imgStreamList");
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
  return(doc.save_file(_rMSIXBin->XML_file.c_str(), "\t", pugi::format_default, pugi::encoding_utf8 ) );
}

//' decodePngStream2IonImages.
//'
//' Obtain a multiple mass channel ion image by decoding the hdd img stream at a specified ionIndex.
//' The MAX operator will be used to merge all ion images in a single image matrix.
//'
//' @param ionIndex the index of ion to extract from the img stream. C style indexing, starting with zero.
//' @param ionCount number of ion image to decode.
//' 
//' @return A NumerixMatrix containing the ion image. //TODO think if row and cols correspond to witdh and height or visaversa
//' 
Rcpp::NumericMatrix rMSIXBin::decodeImgStream2IonImages(unsigned int ionIndex, unsigned int ionCount)
{
  
  //TODO check if the selected   ionCount fits in memory buffer: IONIMG_BUFFER_MB, if not throw exception
  
  //1- Read the complete stream in a buffer
  //Calculate the total number of bytes to read the length vector
  unsigned long byte_count = 0;
  for(int i=ionIndex; i < (ionIndex+ionCount); i++)
  {
    byte_count += _rMSIXBin->iByteLen[i];
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
  Rcpp::NumericMatrix ionImage(img_width, img_height);
  std::vector<unsigned char> raw_image(img_width * img_height, 0);
  double scaling;
  unsigned int png_width, png_height;
  for(int i=0; i < ionCount; i++)
  {
    std::memcpy(&scaling, buffer + _rMSIXBin->iByteOffset[i] - _rMSIXBin->iByteOffset[0], sizeof(float));
  
  unsigned encode_error = lodepng::decode(raw_image, png_width, png_height,
                    (const unsigned char*)(buffer + _rMSIXBin->iByteOffset[i] - _rMSIXBin->iByteOffset[0] + sizeof(float)), _rMSIXBin->iByteLen[i] - sizeof(float),
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
    
    double pixel_value; //Current pixel value
    unsigned int img_offset; //Offset inside the raw image
    unsigned int img_x, img_y; //current coordinates in the image
    for( int iPixel = 0; iPixel < (img_width*img_height); iPixel++)
    {
      //TODO i'm here!!!!! pensa els offsets i coordneas imatge, al final tindras la imatge en un numeric matrix apunt per tornar a R i passar per el raster
      //TODO calculate img_offset!
      //TODO calculate img_x and img_y!!!
      std::memcpy(&pixel_value, raw_image.data() + img_offset, sizeof(imgstreamencoding_type));
      pixel_value *= scaling; 
      ionImage(img_x,img_y) = pixel_value > ionImage(img_x,img_y) ? pixel_value : ionImage(img_x,img_y);
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
Rcpp::List Ccreate_rMSIXBinData(Rcpp::List rMSIobj)
{
  try
  {
    rMSIXBin myXBin(rMSIobj); 
    myXBin.CreateImgStream();
    return myXBin.get_rMSIObj();
  }
  catch(std::runtime_error &e)
  {
    Rcpp::stop(e.what());
  }
  return NULL;
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


///DEBUG METHODS////////////////////////////////////////////////////////////////////////

