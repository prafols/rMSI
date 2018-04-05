/*************************************************************************
 *     rMSI - R package for MSI data processing
 *     Copyright (C) 2018 Pere Rafols Soler
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

#include "pugixml.hpp"
#include <Rcpp.h>
#include <string>
#include <cstring>
#include <algorithm> //To get the transform(), min_element(), max_element() functions
#include <cmath>
#include <libgen.h> //TO get basename() and dirname()
#include <sstream>
using namespace Rcpp;
using namespace pugi;

// imzML resulting data structure:
// $UUID a String object with the UUID
// $continuous_mode a boolean which is true if spectra in continuous mode
// $compression_mz a boolean indicating wheher data is compressed or not
// $compression_int a boolean indicating wheher data is compressed or not
// $mz_dataType a String with the data type: "float", "int"... etc...
// $int_dataType a String with the data type: "float", "int"... etc...
// $pixel_size_um a double with the pixel area? check this...
// $run_data a data.frane with the columns: x, y, mzLength, mzOffset, intLength, intOffset
//
// full path to xml_path must be specified... R function path.expand() can be used 4 this
// [[Rcpp::export]]
List CimzMLParse( String xml_path )
{
  xml_document doc;
  xml_parse_result  result = doc.load_file(xml_path.get_cstring());
  
  if (!result)
  {
    Rcout << "XML [" << xml_path.get_cstring() << "] parsed with errors, attr value: [" << doc.child("node").attribute("attr").value() << "]" << std::endl;
    Rcout << "Error description: " << result.description() << std::endl;
    Rcout << "Error offset: " << result.offset << " (error at [..." << (xml_path.get_cstring() + result.offset) << "]" << std::endl;
    return( List::create(Named("Error") = "XML parse error") );
  }
  
  xml_node mzML = doc.child("mzML");
  if( mzML == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no mzML node found") );
  }
  
  xml_node fileDesc = mzML.child("fileDescription");
  if( fileDesc == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no fileDescription node found") );
  } 
  
  xml_node fileCont = fileDesc.child("fileContent");
  if( fileCont == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no fileContent node found") );
  }
  
  // Data fields
  String sUUID = "";
  String sMD5_Checksum = "";
  String sSHA_Checksum = "";
  bool bContinuous;
  bool bCompressionMz = true; //an error will be raised if this reamains true by the end of xml parsing
  bool bCompressionInt = true; //an error will be raised if this reamains true by the end of xml parsing
  String sMzDataType = "";
  String sIntDataType = "";
  double dPixelSize = 0.0;
  
  // Error control
  bool bUUID_present = false;
  bool bChecksum_present = false;
  bool bDataMode_present = false;
  
  
  std::string accession; //Accession to XML nodes...
  std::string id; //Id accession field used for every xml accession
  std::string value; //Value of item used for every xml accession
  for (xml_node cvParam = fileCont.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
  {
    accession = cvParam.attribute("accession").value();
    value = cvParam.attribute("value").value();
    
    if(accession == "IMS:1000080")
    {
      //Parse the UUID to get just the hex representation in a string (without dashes and {})
      std::size_t ipos = value.find('{');
      if( ipos != std::string::npos )
      {
        value.erase(ipos, 1);
      }
      do
      {
        ipos = value.find('-'); 
        if( ipos != std::string::npos )
        {
          value.erase(ipos, 1);
        } 
      } while ( ipos != std::string::npos);
      ipos = value.find('}');
      if( ipos != std::string::npos )
      {
        value.erase(ipos, 1);
      }
      for( int i=0; i < value.length(); i++)
      {
        value[i] = toupper(value[i]);
      }
      sUUID = value;
      bUUID_present = true;
    }
    if( accession == "IMS:1000091") //SHA Checksum
    {
      transform(value.begin(), value.end(), value.begin(),::toupper);
      sSHA_Checksum = value;
      bChecksum_present = true;
    }
    if( accession == "IMS:1000090") //MD5 Checksum
    {
      transform(value.begin(), value.end(), value.begin(),::toupper);
      sSHA_Checksum = value;
      bChecksum_present = true;
    }
    if( accession == "IMS:1000030")
    {
      bContinuous = true;
      bDataMode_present = true;
    }
    if( accession == "IMS:1000031")
    {
      bContinuous = false;
      bDataMode_present = true;
    }
  }
  
  //Error handling
  if( !bUUID_present )
  {
    return( List::create(Named("Error") = "imzML parse error: UUID not present") );
  }
  if( !bChecksum_present )
  {
    return( List::create(Named("Error") = "imzML parse error: Checksum not present") );
  }
  if( !bDataMode_present )
  {
    return( List::create(Named("Error") = "imzML parse error: Data mode continuous/processed not present") );
  }

  //referenceableParamGroupList
  xml_node refParamGrpLst = mzML.child("referenceableParamGroupList");
  if( refParamGrpLst == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no referenceableParamGroupList node found") );
  }
  for (xml_node refParamGrp = refParamGrpLst.child("referenceableParamGroup"); refParamGrp; refParamGrp = refParamGrp.next_sibling("referenceableParamGroup"))
  {
    id = refParamGrp.attribute("id").value();
    if( id == "mzArray" || id == "intensityArray" || id == "intensities" )
    {
      for (xml_node cvParam = refParamGrp.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
      {
        accession = cvParam.attribute("accession").value();
        
        //No compression
        if( accession == "MS:1000576" &&  id == "mzArray")
        {
          bCompressionMz = false;
        }
        if( accession == "MS:1000576" &&  (id == "intensityArray" || id == "intensities") )
        {
          bCompressionInt = false;
        }
        
        //Check if mass is in valid units
        if( id == "mzArray" && accession == "MS:1000514" )
        {
          if( strcmp(cvParam.attribute("unitAccession").value(), "MS:1000040") != 0 )
          {
            return( List::create(Named("Error") = "m/z Array is not in m/z format, rMSI only supports m/z units") );
          }
        }
       
        if( accession == "IMS:1000101") //External data
        {
          if( !cvParam.attribute("value").as_bool() )
          {
            return( List::create(Named("Error") = "imzML parse error: data must be external (stored in the ibd file)") );
          }
        }
        
        //8 bits integer data array
        //TODO This is not encoded in mzML CV obo file but is specified as valid in imzML spec!
                
        //16 bits integer data array
        //TODO This is not encoded in mzML CV obo file but is specified as valid in imzML spec!
        
        //32 bits integer data array
        if(accession == "IMS:1000141" &&  id == "mzArray")
        {
          sMzDataType = "int";
        }
        if(accession == "IMS:1000141" &&  (id == "intensityArray" || id == "intensities"))
        {
          sIntDataType = "int";
        }
        
        //64 bits integer data array
        if(accession == "IMS:1000142" &&  id == "mzArray")
        {
          sMzDataType = "long";
        }
        if(accession == "IMS:1000142" && (id == "intensityArray" || id == "intensities"))
        {
          sIntDataType = "long";
        }
        
        //32 bits float data array
        if(accession == "MS:1000521" &&  id == "mzArray")
        {
          sMzDataType = "float";
        }
        if(accession == "MS:1000521" &&  (id == "intensityArray" || id == "intensities"))
        {
          sIntDataType = "float";
        }
      
        //64 bits float data array
        if(accession == "MS:1000523" &&  id == "mzArray")
        {
          sMzDataType = "double";
        }
        if(accession == "MS:1000523" &&  (id == "intensityArray" || id == "intensities"))
        {
          sIntDataType = "double";
        }
      }
    }
  }
  
  //Error handling
  if( bCompressionInt || bCompressionMz )
  {
    return( List::create(Named("Error") = "imzML parse error: Data must not be compressed in order to be loaded in rMSI") );
  }
  if( strcmp(sMzDataType.get_cstring(), "") == 0)
  {
    return( List::create(Named("Error") = "imzML parse error: No mass data type found") );
  }
  if( strcmp(sIntDataType.get_cstring(), "") == 0)
  {
    return( List::create(Named("Error") = "imzML parse error: No intensity data type found") );
  }

  //Get the pixel size
  xml_node scanSetLst = mzML.child("scanSettingsList");
  if( scanSetLst == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no scanSettingsList node found") );
  }
  xml_node scanSet = scanSetLst.child("scanSettings");
  if( scanSet == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no scanSettings node found") );
  }
  for (xml_node cvParam = scanSet.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
  {
    accession = cvParam.attribute("accession").value();
    if( accession == "IMS:1000046")
    {
      dPixelSize = sqrt( cvParam.attribute("value").as_double());
    }
  }
  
  //Error handling for pixel size
  if( dPixelSize == 0.0)
  {
    return( List::create(Named("Error") = "imzML parse error: No pixel size found") );
  }

  //Parse binary offsets array  
  xml_node runNode = mzML.child("run");
  if( runNode == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no run node found") );
  }
  xml_node specLst = runNode.child("spectrumList");
  if( specLst == NULL )
  {
    return( List::create(Named("Error") = "XML parse error: no spectrumList node found") );
  }
  
  //Get the total number of pixels as the numbe rof items in spcLst
  const int num_of_pixels = specLst.select_nodes("spectrum").size();
  NumericVector iPos_x(num_of_pixels);
  NumericVector iPos_y(num_of_pixels);
  NumericVector imzLength(num_of_pixels);
  NumericVector imzOffset(num_of_pixels);
  NumericVector iintLength(num_of_pixels);
  NumericVector iintOffset(num_of_pixels);
  
  //Pointers to fast acces to imzLength, imzOffset, iintLength and iintOffset
  NumericVector *ptr_iLength;
  NumericVector *ptr_iOffsets;
  
  int iPixel = 0; //pixel iterator
  bool byMissing, bxMissing;
  bool bArrayPtrNotSet, bLengthMissing, bOffsetMissing;
  for (xml_node spectrum = specLst.child("spectrum"); spectrum; spectrum = spectrum.next_sibling("spectrum")) //For each pixel in the dataset
  {
    //Error handling for iPixel overflow
    if( iPixel >= num_of_pixels)
    {
      return( List::create(Named("Error") = "XML parse error: Pixel counter overflow") );
    }
    
    //Parse pixel coordinates
    bxMissing = true;
    byMissing = true;
    xml_node scan = spectrum.child("scanList").child("scan");
    for (xml_node cvParam = scan.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
    {
      accession = cvParam.attribute("accession").value();
      if( accession == "IMS:1000050" )
      {
        iPos_x[iPixel] = cvParam.attribute("value").as_double(); 
        bxMissing = false;
      }
      if( accession == "IMS:1000051" )
      {
        iPos_y[iPixel] = cvParam.attribute("value").as_double();
        byMissing = false;
      }
    }
    
    //Error handling for missing XY values
    if( bxMissing || byMissing )
    {
      return( List::create(Named("Error") = "XML parse error: Some pixel is missing coordinates") );
    }
    
    //Parse arrays offsets
    bArrayPtrNotSet = true;
    xml_node binArrayLst = spectrum.child("binaryDataArrayList");
    for (xml_node binArray = binArrayLst.child("binaryDataArray"); binArray; binArray = binArray.next_sibling("binaryDataArray"))
    {
      std::string sArrayType = binArray.child("referenceableParamGroupRef").attribute("ref").value();
      if( sArrayType == "mzArray" )
      {
        ptr_iLength = &imzLength;
        ptr_iOffsets = &imzOffset;
        bArrayPtrNotSet = false;
      }
      if( sArrayType == "intensityArray" || sArrayType == "intensities")
      {
        ptr_iLength = &iintLength;
        ptr_iOffsets = &iintOffset;
        bArrayPtrNotSet = false;
      }
      
      //Error handling no ptr set so no data available
      if( bArrayPtrNotSet )
      {
        return( List::create(Named("Error") = "XML parse error: Missing binary data reference") );
      }
      
      bLengthMissing = true;
      bOffsetMissing = true;
      for (xml_node cvParam = binArray.child("cvParam"); cvParam; cvParam = cvParam.next_sibling("cvParam"))
      {
        accession = cvParam.attribute("accession").value();
        if( accession == "IMS:1000103")
        {
          (*ptr_iLength)[iPixel] = cvParam.attribute("value").as_double();
          bLengthMissing = false;
        }
        if( accession == "IMS:1000102")
        {
          (*ptr_iOffsets)[iPixel] = cvParam.attribute("value").as_double();
          bOffsetMissing = false;
        }
      }
      
      //Error handling no ptr set so no data available
      if( bLengthMissing ||  bOffsetMissing)
      {
        return( List::create(Named("Error") = "XML parse error: Missing binary data assignment") );
      }
      
    }
    
    iPixel++;
  }
  
  return List::create(Named("UUID") = sUUID , 
                      Named("SHA") = sSHA_Checksum,
                      Named("MD5") = sMD5_Checksum,
                      Named("continuous_mode")= bContinuous,
                      Named("compression_mz")= bCompressionMz,
                      Named("compression_int")= bCompressionInt,
                      Named("mz_dataType") = sMzDataType,
                      Named("int_dataType") = sIntDataType,
                      Named("pixel_size_um") = dPixelSize,
                      Named("run_data") = DataFrame::create( Named("x") = iPos_x,
                                                             Named("y") = iPos_y, 
                                                             Named("mzLength") = imzLength, 
                                                             Named("mzOffset") = imzOffset, 
                                                             Named("intLength") = iintLength, 
                                                             Named("intOffset") = iintOffset)
                      );
}


//Append the cvParam node according the data type to a current xml_node specified by a pointer.
//Return the number of bytes used to encode each datapoint or -1 if error.
int AppendimzMLDataTypeNode(const std::string sDataType, pugi::xml_node *imzMLnode)
{
  pugi::xml_node cvParam;
  int encodingByteSize = -1;
  
  if( sDataType.compare("int") == 0 )
  {
    cvParam = imzMLnode->append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000141";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "32-bit integer";
    encodingByteSize = 4;
  }
  else if( sDataType.compare("long") == 0 )
  {
    cvParam = imzMLnode->append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000142";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "64-bit integer";
    encodingByteSize = 8;
  }
  else if( sDataType.compare("float") == 0 )
  {
    cvParam = imzMLnode->append_child("cvParam");
    cvParam.append_attribute("accession") = "MS:1000521";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "32-bit float";
    encodingByteSize = 4;
  }
  else if( sDataType.compare("double") == 0 )
  {
    cvParam = imzMLnode->append_child("cvParam");
    cvParam.append_attribute("accession") = "MS:1000523";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "64-bit float";
    encodingByteSize = 8;
  }
  else
  {
    stop("Error: No valid imzML data type\n");
  }
  
  return encodingByteSize;
}

// imzML imgInfo data structure:
// $UUID a String object with the UUID
// $continuous_mode a boolean which is true if spectra in continuous mode
// $compression_mz a boolean indicating wheher data is compressed or not
// $compression_int a boolean indicating wheher data is compressed or not
// $mz_dataType a String with the data type: "float", "int"... etc...
// $int_dataType a String with the data type: "float", "int"... etc...
// $pixel_size_um a double with the pixel area? check this...
// $run_data a data.frane with the columns: x, y, mzLength, mzOffset, intLength, intOffset
//
// full path to xml_path must be specified... R function path.expand() can be used 4 this
// [[Rcpp::export]]
bool CimzMLStore( String fname, List imgInfo )
{
  //Reusable pugi variables
  pugi::xml_node cvParam;
 
  //Extract the RunData Matrix for fast access
  DataFrame runMat = as<DataFrame>(imgInfo["run_data"]);
  
  // empty xml document with custom declaration node
  pugi::xml_document doc;
  pugi::xml_node decl = doc.prepend_child(pugi::node_declaration);
  decl.append_attribute("version") = "1.0";
  decl.append_attribute("encoding") = "UTF-8";
  decl.append_attribute("standalone") = "no";
  
  //mzML top level node
  pugi::xml_node node_mzML = doc.append_child("mzML");
  node_mzML.append_attribute("version") =  "1.1";
  node_mzML.append_attribute("xmlns") = "http://psi.hupo.org/ms/mzml";
  node_mzML.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
  node_mzML.append_attribute("xsi:schemaLocation") = "http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0_idx.xsd";

  //cvList top node
  pugi::xml_node node_cvList = node_mzML.append_child("cvList");
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
  
  //fileDescription node
  pugi::xml_node node_fdesc = node_mzML.append_child("fileDescription");
  
  //fileContent
  pugi::xml_node fileContent = node_fdesc.append_child("fileContent");
  cvParam = fileContent.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000080";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "universally unique identifier";
  
  std::string sUUIDorig = (std::string)as<String>(imgInfo["UUID"]);
  std::string sUUIDparsed = "{";
  sUUIDparsed += sUUIDorig.substr(0, 8);
  sUUIDparsed += "-";
  sUUIDparsed += sUUIDorig.substr(8, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUIDorig.substr(12, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUIDorig.substr(16, 4); 
  sUUIDparsed += "-";
  sUUIDparsed += sUUIDorig.substr(20, 12); 
  sUUIDparsed += "}";
  cvParam.append_attribute("value") = sUUIDparsed.c_str();
  cvParam = fileContent.append_child("cvParam");
  if( as<bool>(imgInfo["continuous_mode"]) )
  {
    cvParam.append_attribute("accession") = "IMS:1000030";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "continuous";
  }
  else
  {
    cvParam.append_attribute("accession") = "IMS:1000031";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "processed";
  }

  if( ((std::string)as<String>(imgInfo["MD5"])).length() > 0 )
  {
    cvParam = fileContent.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000090";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "ibd MD5";
    cvParam.append_attribute("value") = as<String>(imgInfo["MD5"]).get_cstring();
  }
  
  if( ((std::string)as<String>(imgInfo["SHA"])).length() > 0)
  {
    cvParam = fileContent.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000091";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "ibd SHA-1";
    cvParam.append_attribute("value") = as<String>(imgInfo["SHA"]).get_cstring();
  }
  
  cvParam = fileContent.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000294";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "mass spectrum";
  
  //sourceFileList
  pugi::xml_node srcFileLst = node_fdesc.append_child("sourceFileList");
  srcFileLst.append_attribute("count") = "1";
  pugi::xml_node srcFile = srcFileLst.append_child("sourceFile");
  srcFile.append_attribute("id") = "sourceFile0";
  
  //Dirname and Basenme must be copied to strings since the POSIX fucntions may modify the content of the char pointer.
  std::string sDirname = fname.get_cstring();
  sDirname = dirname((char*)sDirname.c_str());
  std::string sBasename = fname.get_cstring();
  sBasename = basename((char*)sBasename.c_str());
  srcFile.append_attribute("location") = sDirname.c_str();
  srcFile.append_attribute("name") = sBasename.c_str();
  
  cvParam = srcFile.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000560";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "mass spectrometer file format";
  cvParam.append_attribute("value") = "rMSI exported imzML";
  
  cvParam = srcFile.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1002333";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "conversion software";

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
  
  //referenceableParamGroupList
  pugi::xml_node refParamGrpLst = node_mzML.append_child("referenceableParamGroupList");
  refParamGrpLst.append_attribute("count") = "2";
  
  //The MZ Array
  pugi::xml_node refParamGrpMzArray = refParamGrpLst.append_child("referenceableParamGroup");
  refParamGrpMzArray.append_attribute("id") = "mzArray";
  cvParam = refParamGrpMzArray.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000514";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "m/z array";
  cvParam.append_attribute("unitAccession") = "MS:1000040";
  cvParam.append_attribute("unitCvRef") = "MS";
  cvParam.append_attribute("unitName") = "m/z";
  int mzEncodingBytes = AppendimzMLDataTypeNode(as<String>(imgInfo["mz_dataType"]).get_cstring(), &refParamGrpMzArray);
  if( as<bool>(imgInfo["compression_mz"]) )
  {
    stop("Error: binary data compression is not supported\n");
  }
  else
  {
    cvParam = refParamGrpMzArray.append_child("cvParam");
    cvParam.append_attribute("accession") = "MS:1000576";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "no compression";
  }
  cvParam = refParamGrpMzArray.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000101";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "external data";
  cvParam.append_attribute("value") = "true";
  
  //The Intensities Array
  pugi::xml_node refParamGrpIntArray = refParamGrpLst.append_child("referenceableParamGroup");
  refParamGrpIntArray.append_attribute("id") = "intensityArray";
  cvParam = refParamGrpIntArray.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000515";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "intensity array";
  cvParam.append_attribute("unitAccession") = "MS:1000131";
  cvParam.append_attribute("unitCvRef") = "MS";
  cvParam.append_attribute("unitName") = "number of detector counts";
  int intEncodingBytes = AppendimzMLDataTypeNode(as<String>(imgInfo["int_dataType"]).get_cstring(), &refParamGrpIntArray);
  if( as<bool>(imgInfo["compression_int"]) )
  {
    stop("Error: binary data compression is not supported\n");
  }
  else
  {
    cvParam = refParamGrpIntArray.append_child("cvParam");
    cvParam.append_attribute("accession") = "MS:1000576";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "no compression";
  }
  cvParam = refParamGrpIntArray.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000101";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "external data";
  cvParam.append_attribute("value") = "true";
  
  //sampleList
  pugi::xml_node node_sampleLst = node_mzML.append_child("sampleList");
  node_sampleLst.append_attribute("count") = "1";
  pugi::xml_node node_sample = node_sampleLst.append_child("sample");
  node_sample.append_attribute("id") = "sample1";
  node_sample.append_attribute("name") = "Sample1";
  cvParam = node_sample.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000001";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "sample number";
  cvParam.append_attribute("value") = "1";
  
  //softwareList
  pugi::xml_node node_softLst = node_mzML.append_child("softwareList");
  node_softLst.append_attribute("count") = "1";
  pugi::xml_node node_soft = node_softLst.append_child("software");
  node_soft.append_attribute("id") = "rMSI";
  node_soft.append_attribute("version") = "0.8";
  cvParam = node_soft.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000799";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "custom unreleased software tool";
  cvParam.append_attribute("value") = "rMSI";

  
  //scanSettingsList
  pugi::xml_node node_scanSetLst = node_mzML.append_child("scanSettingsList");
  node_scanSetLst.append_attribute("count") = "1";
  pugi::xml_node node_scanSet = node_scanSetLst.append_child("scanSettings");
  node_scanSet.append_attribute("id") = "scanSettings0";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000401";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "top down";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000413";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "flyback";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000491";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "linescan left right";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000480";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "horizontal line scan";
  
  IntegerVector xCoords = runMat["x"];
  int maxXCoord = *std::max_element( xCoords.begin(), xCoords.end()); //Using the max in coords as Bruker is using it this way... no matter if there is an offset. 
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000042";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "max count of pixels x";
  cvParam.append_attribute("value") = maxXCoord;
  
  IntegerVector yCoords = runMat["y"];
  int maxYCoord = *std::max_element( yCoords.begin(), yCoords.end());  //Using the max in coords as Bruker is using it this way... no matter if there is an offset.
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000043";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "max count of pixels y";
  cvParam.append_attribute("value") = maxYCoord;
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000044";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "max dimension x";
  cvParam.append_attribute("unitAccession") = "UO:0000017";
  cvParam.append_attribute("unitCvRef") = "UO";
  cvParam.append_attribute("unitName") = "micrometer";
  cvParam.append_attribute("value") = round((double)maxXCoord * as<double>(imgInfo["pixel_size_um"]));
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000045";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "max dimension y";
  cvParam.append_attribute("unitAccession") = "UO:0000017";
  cvParam.append_attribute("unitCvRef") = "UO";
  cvParam.append_attribute("unitName") = "micrometer";
  cvParam.append_attribute("value") = round((double)maxYCoord * as<double>(imgInfo["pixel_size_um"]));
    
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000053";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "absolute position offset x";
  cvParam.append_attribute("unitAccession") = "UO:0000017";
  cvParam.append_attribute("unitCvRef") = "UO";
  cvParam.append_attribute("unitName") = "micrometer";
  cvParam.append_attribute("value") = "0";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000054";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "absolute position offset y";
  cvParam.append_attribute("unitAccession") = "UO:0000017";
  cvParam.append_attribute("unitCvRef") = "UO";
  cvParam.append_attribute("unitName") = "micrometer";
  cvParam.append_attribute("value") = "0";
  
  cvParam = node_scanSet.append_child("cvParam");
  cvParam.append_attribute("accession") = "IMS:1000046";
  cvParam.append_attribute("cvRef") = "IMS";
  cvParam.append_attribute("name") = "pixel size";
  cvParam.append_attribute("value") = pow(as<double>(imgInfo["pixel_size_um"]), 2.0);
  
  //instrumentConfigurationList, just declaring it was a file exported from rMSI
  pugi::xml_node node_instrtLst = node_mzML.append_child("instrumentConfigurationList");
  node_instrtLst.append_attribute("count") = "1";
  pugi::xml_node node_instrCfg = node_instrtLst.append_child("instrumentConfiguration");
  node_instrCfg.append_attribute("id") = "instrumentConfiguration0";
  
  //dataProcessingList
  pugi::xml_node node_dataProcLst = node_mzML.append_child("dataProcessingList");
  node_dataProcLst.append_attribute("count") = "1";
  pugi::xml_node node_dataProc = node_dataProcLst.append_child("dataProcessing");
  node_dataProc.append_attribute("id") = "dataProcessing0";
  pugi::xml_node node_procMethod = node_dataProc.append_child("processingMethod");
  node_procMethod.append_attribute("order") = "0";
  node_procMethod.append_attribute("softwareRef") = "software0";
  cvParam = node_procMethod.append_child("cvParam");
  cvParam.append_attribute("accession") = "MS:1000544";
  cvParam.append_attribute("cvRef") = "MS";
  cvParam.append_attribute("name") = "Conversion to mzML";
  
  //Run data
  int total_num_pixels = (int) runMat.nrows();
  pugi::xml_node node_run = node_mzML.append_child("run");
  node_run.append_attribute("defaultInstrumentConfigurationRef") = "instrumentConfiguration0";
  node_run.append_attribute("id") = "run0";
  node_run.append_attribute("sampleRef") = "sample1";
  pugi::xml_node node_spectrumLst = node_run.append_child("spectrumList");
  node_spectrumLst.append_attribute("count") = total_num_pixels;
  node_spectrumLst.append_attribute("defaultDataProcessingRef") = "dataProcessing0";
  
  std::stringstream ssParser; //Auxiliar string to parse scans
  //Copy to numeric to access it using a long long range 
  NumericVector v_mzLength = runMat["mzLength"];
  NumericVector v_mzOffset = runMat["mzOffset"];
  NumericVector v_intLength = runMat["intLength"];
  NumericVector v_intOffset = runMat["intOffset"];
  
  pugi::xml_node node_spectrum; //Reusable spectrum node
  pugi::xml_node node_refParamGrpRef; //Reusable referenceableParamGroupRef node
  pugi::xml_node node_scanLst; //Reusable scanList node
  pugi::xml_node node_scan; //Reusable scan node
  pugi::xml_node node_binDataLst; //Reusable binaryDataArrayList node
  pugi::xml_node node_binData; //Reusable binaryDataArray node
  pugi::xml_node  node_binary; //Reusable binary node
  for( int i = 0; i < total_num_pixels; i++)
  {
    //Rcout<<"Parsing pixel "<<i <<" of "<< total_num_pixels << "\n";
    ssParser.str("");
    ssParser << "Scan=" <<i;
    node_spectrum = node_spectrumLst.append_child("spectrum");
    node_spectrum.append_attribute("defaultArrayLength") = "0";
    node_spectrum.append_attribute("id") = ssParser.str().c_str();
    node_spectrum.append_attribute("index") = i;
    
    node_refParamGrpRef = node_spectrum.append_child("referenceableParamGroupRef");
    node_refParamGrpRef.append_attribute("ref") = "spectrum";
    
    node_scanLst = node_spectrum.append_child("scanList");
    node_scanLst.append_attribute("count") = "1";
    cvParam = node_scanLst.append_child("cvParam");
    cvParam.append_attribute("accession") = "MS:1000795";
    cvParam.append_attribute("cvRef") = "MS";
    cvParam.append_attribute("name") = "no combination";
    
    node_scan =  node_scanLst.append_child("scan");
    node_scan.append_attribute("instrumentConfigurationRef") = "instrumentConfiguration0";
    cvParam = node_scan.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000050";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "position x";
    cvParam.append_attribute("value") = xCoords[i];
    cvParam = node_scan.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000051";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "position y";
    cvParam.append_attribute("value") = yCoords[i];
    
    node_binDataLst = node_spectrum.append_child("binaryDataArrayList");
    node_binDataLst.append_attribute("count") = "2";
    
    //MZ array
    node_binData = node_binDataLst.append_child("binaryDataArray");
    node_binData.append_attribute("encodedLength") = "0";
    node_refParamGrpRef = node_binData.append_child("referenceableParamGroupRef");
    node_refParamGrpRef.append_attribute("ref") = "mzArray";
    
    cvParam = node_binData.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000103";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "external array length";
    cvParam.append_attribute("value") = (unsigned long long)v_mzLength[i];
    
    cvParam = node_binData.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000104";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "external encoded length";
    cvParam.append_attribute("value") = (unsigned long long) (mzEncodingBytes*v_mzLength[i]);
    
    cvParam = node_binData.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000102";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "external offset";
    cvParam.append_attribute("value") = (unsigned long long)v_mzOffset[i];
    
    node_binary = node_binData.append_child("binary");
    
    //INT array
    node_binData = node_binDataLst.append_child("binaryDataArray");
    node_binData.append_attribute("encodedLength") = "0";
    node_refParamGrpRef = node_binData.append_child("referenceableParamGroupRef");
    node_refParamGrpRef.append_attribute("ref") = "intensityArray";
    
    cvParam = node_binData.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000103";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "external array length";
    cvParam.append_attribute("value") = (unsigned long long)v_intLength[i];
    
    cvParam = node_binData.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000104";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "external encoded length";
    cvParam.append_attribute("value") = (unsigned long long)(intEncodingBytes*v_intLength[i]);
    
    cvParam = node_binData.append_child("cvParam");
    cvParam.append_attribute("accession") = "IMS:1000102";
    cvParam.append_attribute("cvRef") = "IMS";
    cvParam.append_attribute("name") = "external offset";
    cvParam.append_attribute("value") = (unsigned long long)v_intOffset[i];
    
    node_binary = node_binData.append_child("binary");
  }
  
 // save document to file
  return(doc.save_file(fname.get_cstring(), "\t", pugi::format_default, pugi::encoding_utf8 ) );
}
